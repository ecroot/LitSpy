import requests
import re
import itertools
import threading
import multiprocessing
import copy

from collections import Counter
from bs4 import BeautifulSoup
from rtgo import ReadyThready
from urllib import parse

from litspy.get_synonyms import ExtractOLSSynonyms, OLSRequests, GetUniprotSynonyms, ArgumentCleaner
from litspy.alternative_characters import hyphens, greek_dict, numerals
from litspy.noisy_phrases import common_gene_noise, stop_words


class Query:
    def __init__(self, df, logger, disease=None, tissue=None, kwds=None, others=None, expand=None,
                 n_threads=multiprocessing.cpu_count()):
        """
        initialise object with relevant values

        :param pandas.DataFrame df: dataframe containing gene name, gene id type, taxonomy id, keywords
        :param logging.Logger logger: project logger
        :param list disease: a single-element list containing the specified disease
        :param list tissue: a single-element list containing the specified tissue
        :param str kwds: a comma-separated string of specified keywords
        :param dict others: a dictionary of EPMC field: value
        :param bool expand: whether to get synonyms for keywords
        """
        self.logger = logger
        self.df = df
        self.disease = disease
        self.tissue = tissue
        self.kwds = kwds
        self.other_args = others
        self.expand = expand
        self.n_threads = n_threads
        self.greek_chars = []
        for chr_list in greek_dict.values():
            self.greek_chars.extend(chr_list)
        self.disease_syns = []
        self.tissue_syns = []
        self.kwd_syns = []
        self.kwd_syn_lists = []
        self.gene_syns = []
        self.gene_syn_roots = []

    @staticmethod
    def add_abstract_title_and_join_synonyms(syns, search_in_kwds=True, join_on_and=False):
        """
        add parameters for EMPC search in abstract and title regions

        :param list syns: list of synonyms
        :param bool search_in_kwds: whether to search in the keywords list for the term (false for gene roots)
        :param bool join_on_and: whether to join terms on AND (default false, to join on OR)
        :return: string of joined synonyms with relevant title and abstract parameters
        """
        # initialise
        exp_syns = []

        # create e.g. ABSTRACT:"synonym" OR TITLE:"synonym" for every synonym, and add to a list
        for syn in syns:
            # if the synonym already has the fields added (possible for keyword synonyms), add it to the list
            if syn.startswith("(TITLE:"):
                exp_syns.append(syn.strip("()"))  # remove brackets to avoid extra brackets getting added later
            # else, create the field query and add it to the list
            elif search_in_kwds:
                exp_syns.append(f'TITLE:"{syn}" OR KW:"{syn}" OR ABSTRACT:"{syn}"')
            else:
                exp_syns.append(f'TITLE:"{syn}" OR ABSTRACT:"{syn}"')

        # join the items in the list to create e.g. (ABSTRACT:"1" OR TITLE:"1" OR ABSTRACT:"2" OR TITLE:"2")
        if join_on_and:
            query_string = '(' + ") & (".join(exp_syns) + ')'
        else:
            query_string = '(' + " OR ".join(exp_syns) + ')'

        return query_string

    def clean_orig_term(self, term):
        """
        ensure the term is a string and return it

        :param list or str term: the term supplied at the command line
        :return: str term
        :rtype: str
        """
        # make sure the term is a string
        cleaner = ArgumentCleaner(term=term, logger=self.logger)
        clean_term = cleaner.get_clean_arg()
        return clean_term

    def get_syns_in_parallel(self, term, req, iris, is_tissue=False):
        """
        for each iri in the list, parse the json and extract the synonyms

        :param str term: cleaned input term
        :param req: OLSRequests object initialised with an input term and logger
        :param list iris: list of iris derived from querying for the term in a relevant ontology
        :param bool is_tissue: whether to collect anatomy qualifiers
        :return: list of synonyms
        """
        # initialise job handler and synonyms objects, and list of synonyms
        ReadyThready.set_logger(self.logger)
        synonyms = ExtractOLSSynonyms(term, self.logger, n_threads=self.n_threads)
        syns = [term]

        self.logger.info(f"{threading.current_thread().name}: Collecting synonyms for '{term}'")

        # for each iri, retrieve the information from the page as a dict and extract synonyms from it
        if len(iris) > 1:
            parsed_json = ReadyThready.go(req.get_json_for_iri, [iris], 0, n_threads=self.n_threads)
            syn_results = ReadyThready.go(synonyms.get_syns, [parsed_json, is_tissue], 0, n_threads=self.n_threads)

            # add collected syns to the syn list
            for syn_list in syn_results:
                syns.extend(syn_list)
        else:
            for iri in iris:
                parsed_json = req.get_json_for_iri(iri)
                syns.extend(synonyms.get_syns(parsed_json, is_tissue))

        if is_tissue:
            # check the anatomy qualifiers dictionary for further synonyms and add them to the synonyms list
            syns = synonyms.get_anatomy_qualifiers(syns)

        self.logger.info(f"{threading.current_thread().name}: Done collecting synonyms for '{term}'")

        # clean the list of synonyms
        self.logger.info(f"{threading.current_thread().name}: Cleaning collected synonyms for '{term}'")
        if is_tissue:
            clean_syns = synonyms.clean_syn_list(syns_list=syns, syn_type="tissue")
        else:
            clean_syns = synonyms.clean_syn_list(syns_list=syns)
        self.logger.info(f"{threading.current_thread().name}: Done: {len(clean_syns)} cleaned synonyms for '{term}'")

        return clean_syns

    def get_disease_syns(self):
        """
        collect synonyms for the disease term

        :return: list of synonyms for the supplied disease term
        :rtype: list
        """
        # ensure the disease term is a string
        clean_disease = self.clean_orig_term(self.disease)
        self.disease = clean_disease

        # initialise requests object with the disease string
        arg_requests = OLSRequests(clean_disease, self.logger)

        # get EBI OLS IRIs for the disease
        self.logger.info(f"{threading.current_thread().name}: Getting synonyms for '{clean_disease}'")
        iris = arg_requests.get_iris('&ontology=mondo')

        # get the synonyms
        syns = self.get_syns_in_parallel(clean_disease, arg_requests, iris)

        self.disease_syns = syns

    def get_tissue_syns(self):
        """
        collect synonyms for the tissue term

        :return: non-redundant synonyms for the supplied tissue term
        :rtype: list
        """
        # ensure the tissue term is a string
        clean_tissue = self.clean_orig_term(self.tissue)
        self.tissue = clean_tissue

        # initialise requests object with the disease string
        arg_requests = OLSRequests(clean_tissue, self.logger)

        # get EBI OLS IRIs for UBERON ontology nodes that exactly match the tissue
        self.logger.info(f"{threading.current_thread().name}: Getting synonyms for '{clean_tissue}'")
        iris = arg_requests.get_iris("&exact=on&ontology=uberon")

        # get the synonyms
        syns = self.get_syns_in_parallel(clean_tissue, arg_requests, iris, is_tissue=True)

        self.tissue_syns = syns

    def get_kwd_synonyms(self, kwd_string):
        """
        collect synonyms for a keyword

        :param str kwd_string: comma-separated list of keywords (can be a phrases rather than a single words)
        :return: list of lists of synonyms for the keyword
        """
        kwd_syns = []
        raw_kwds = kwd_string.split(',')
        kwds = [x.strip() for x in raw_kwds]
        kwds = list(set(kwds))

        # get synonyms for the keywords if specified, add keywords and their synonyms to the list of kwd syns
        if self.expand:
            for kwd in kwds:
                # initialise requests object with the keyword
                arg_requests = OLSRequests(kwd, self.logger)

                # get EBI OLS IRIs for the kwd
                self.logger.info(f"{threading.current_thread().name}: Getting synonyms for keyword '{kwd}'")
                iris = arg_requests.get_iris('&exact=on')

                # get the synonyms
                syns = self.get_syns_in_parallel(kwd, arg_requests, iris)

                if len(syns) > 100:
                    self.logger.warning(f"Noise suspected: Too many synonyms ({len(syns)}) identified for keyword "
                                        f"'{kwd}'. Only the supplied keyword will be used to search Europe "
                                        f"PMC, not its synonyms")
                    kwd_syns.append([kwd])
                else:
                    kwd_syns.append(syns)

        else:
            for kwd in kwds:
                kwd_syns.append([kwd])

        self.kwd_syn_lists = kwd_syns

        return

    def get_constant_query_strings(self):
        """
        for the universal parameters (disease, tissue, other), obtain synonyms where relevant, and make query strings

        :return: query strings for the supplied disease, tissue and other parameters
        :rtype: tuple[str, str, str, list, list]
        """
        # initialise to empty strings
        disease_query_string = ""
        tissue_query_string = ""
        others_query_string = ""
        kwd_query_strings = []
        kwds = []

        # if EPMC search fields were supplied, format them in to a sting and check that the string is not too long
        if self.other_args:
            list_of_other_args = []
            # create a list of strings in the format "parameter:value" from the keys and values in the other_args dict
            for k, v in self.other_args.items():
                list_of_other_args.append(f"{k}:{v}")
            # join the list on ampersands to create a query string
            others_query_string = ' & '.join(list_of_other_args)
            # this query string contains settings, which can not be split in to multiple queries if too long. Therefore,
            # log an informative error if the query string is too long, but do not exit; continue with empty string
            if len(others_query_string) > 4500:
                self.logger.error(f"The string of settings was too long to produce a query. No settings will be "
                                  f"included in the queries. Please consider removing some of the settings from the "
                                  f"following in order to allow the query to run with settings: {others_query_string}")
                others_query_string = ""

        # if a tissue was supplied, perform synonym expansion and return a string for querying
        if self.tissue:
            self.get_tissue_syns()
            tissue_query_string = self.add_abstract_title_and_join_synonyms(self.tissue_syns)

        # if a disease was supplied, perform synonym expansion and return a string for querying
        if self.disease:
            self.get_disease_syns()
            disease_query_string = self.add_abstract_title_and_join_synonyms(self.disease_syns)

        if self.kwds:
            # store lists of keywords for each keyword
            self.get_kwd_synonyms(self.kwds)

            for each in self.kwd_syn_lists:
                # create a single list combining all keyword syns
                self.kwd_syns.extend(each)
                # construct the query string for the syns of the keyword
                kwd_query_strings.append(self.add_abstract_title_and_join_synonyms(each))
            raw_kwds = self.kwds.split(',')
            kwds = [x.strip() for x in raw_kwds]

        return disease_query_string, tissue_query_string, others_query_string, kwd_query_strings, kwds

    def get_keyword_query_string_and_words(self, row):
        """
        get synonyms for the keywords (if specified) and make query strings from the keywords

        :param list[str] row: row from a df of query inputs represented as a list
        :return: query parts for the keywords, and a list of keywords for the output file
        :rtype: tuple[list, list, list]
        """
        kwds = []
        kwd_query_strings = []
        kwd_flatlist = []

        # split keywords in to a list
        if row[4]:
            self.logger.info(f"{threading.current_thread().name}: Creating query parts for keyword(s) '{row[4]}'")
            raw_kwds = row[4].split(',')
            kwds = [x.strip() for x in raw_kwds]
            self.get_kwd_synonyms(row[4])
            for each_list in self.kwd_syn_lists:
                kwd_query_string = self.add_abstract_title_and_join_synonyms(each_list)
                kwd_query_strings.append(kwd_query_string)

            self.logger.info(f"{threading.current_thread().name}: Done: successfully created query parts for the "
                             f"keyword(s) '{row[4]}'")

            kwd_flatlist = [item for sublist in self.kwd_syn_lists for item in sublist]

        return kwd_query_strings, kwds, kwd_flatlist

    def get_search_terms_string(self, kwds):
        """
        join the keywords, disease term, tissue term and other search settings in to a single string for each search,
        to be printed in to the output files

        :param list kwds: list of keywords
        :return: string of search terms for output files
        """
        clean_terms = []
        # join the list of keywords on commas
        kwds = ', '.join(kwds)

        for each in [self.disease, self.tissue, self.other_args, kwds]:
            # if the item is a dict, then create a string from each key and value, add these to a list, join the
            # list on commas to get a string of search terms and re-assign this string to the variable
            if isinstance(each, dict):
                new_terms = []
                for k, v in each.items():
                    new_terms.append(f"{k} {v}")
                each = ', '.join(new_terms)

            if each:
                each = each.replace("'", "")  # remove single quote marks
                each = re.sub(r"\s{2,}", " ", each)  # multiple spaces with a single space
                clean = each.strip(", ")  # remove trailing spaces and commas
                # add the search term string to the list of clean search terms
                clean_terms.append(clean)

        # create and return the full search term string from the list of clean search terms
        kwds_disease_tissue = f"{', '.join(clean_terms)}"
        return kwds_disease_tissue

    def final_gene_clean(self, gene, all_syns):
        """
        perform the standard cleaning and expansion on the synonym list, and also remove any common gene noise and from
        the list (common biomedical abbreviations, short words etc.)

        :param str gene: gene name
        :param list all_syns: list of synonyms
        :return: list: cleaned synonym list
        """
        # initialise
        filtered_syns = []
        cleaner = ExtractOLSSynonyms(gene, self.logger)

        # perform the normal synonym cleaning
        first_clean = cleaner.clean_syn_list(syns_list=all_syns)

        for syn in first_clean:
            # filter common noise and two-character synonyms
            if len(syn) > 2 and not \
                    re.fullmatch(r"([A-z]|CI|CD|CT|CRP|PP|LAG|PER|period|TC|UP) \d+\s?\d*", syn, re.IGNORECASE) \
                    and not re.fullmatch(r"[vVlL]\d+\s?\d*", syn) and syn.upper() not in common_gene_noise:
                # filter phrases that start or end with stop words
                split_syn = re.split(f"{'|'.join(hyphens)}| ", syn.lower())
                last_part = split_syn[-1]
                first_part = split_syn[0]
                if last_part not in stop_words and first_part not in stop_words:
                    filtered_syns.append(syn)

        # ensure the original term is still in the list following cleaning
        filtered_syns.insert(0, gene)
        final_syns = list(set(filtered_syns))

        return final_syns

    def get_recurrent_root_phrases(self, syns):
        """
        find which root phrases recur among synonyms for wildcard genes, and return the recurring phrases

        :param list syns: list of synonyms for a wildcard gene
        :return: list of root phrases from supplied synonyms
        """
        # de-duplicate
        syns = list(set(syns))
        root_phrases = []
        recurrent_root_phrases = []
        for syn in syns:
            root_phrase = self.get_root_name_from_gene_syn(syn)
            if root_phrase:
                root_phrases.append(root_phrase)
        # count occurrences of number-stripped phrases, and if any are present more than once then return them
        root_counts = Counter(root_phrases)
        for k, v in root_counts.items():
            if v > 1:
                recurrent_root_phrases.append(k)
        return recurrent_root_phrases

    @staticmethod
    def get_root_name_from_gene_syn(syn):
        """
        strip the gene synonym down to its root

        :param str syn: synonym for a gene
        :rtype: str
        :return: root phrase of the synonym
        """
        root_phrase = ""
        pattern = re.compile(fr".+?[{''.join(hyphens)}\s,]*(\d| [{''.join(numerals)}]+)+[A-z]?(\d| [{''.join(numerals)}]+)*$")

        if re.search(pattern, syn) and not re.match(r"[Cc]\d+orf\d+", syn) and not re.match(r"UNQ\d+/PRO\d+", syn) \
                and not re.match(r"KIAA\d+", syn) and not syn.startswith('type'):
            match = re.match(
                fr"(.+?)(TYPE|[Tt]ype)?[{''.join(hyphens)}\s,]*(\d| [{''.join(numerals)}]+)+([{''.join(hyphens)}\s,]|[A-z]|MOTIF|[Mm]otif|PROTEIN|[Pp]rotein|DOMAIN|[Dd]omain|PSEUDOGENE|[Pp]seudogene|CONTAINING|[Cc]ontaining)*(\d| [{''.join(numerals)}]+)*$",
                syn)
            root_phrase = match.group(1)
            root_phrase.strip(',- ')

        return root_phrase

    @staticmethod
    def systematic_family_naming_test(gene):
        """
        determine whether a gene's name suggests it is part of a similarly named family/group, and should therefore
        be searched as a partial name in lists (e.g. ADAMTS5 in the phrase "ADAMTS4 and 5")

        :param str gene: a gene name (not fully written out)
        :return: bool: true if gene ends with digits/similar
        """
        # if the synonym is not an ORF, UNQ or KIAA type name (these are not families), return true if the synonym ends
        # with numbers (and optionally capital letters)
        if not re.match(r"[Cc]\d+orf\d+", gene) and not re.match(r"UNQ\d+/PRO\d+", gene) and not \
                re.match(r"KIAA\d+", gene):
            if re.match(r".*\d+[A-z]?\d*$", gene, re.IGNORECASE):
                return True
        else:
            return False

    def get_gene_syn_iris_from_ebi_ols(self, syn, search_human_only):
        """
        get a list of iris for the gene synonym, using exact settings if not wildcard and searching under the human
        node if no other organism was specified

        :param str syn: gene synonym
        :param bool search_human_only: whether to only consider human genes
        :rtype: list[str]
        :return: list of iri strings
        """
        req = OLSRequests(syn, self.logger)
        if '*' in syn:
            # if the supplied synonym contains a wildcard character, do not query OLS with exact mode
            if search_human_only:
                iris = req.get_iris(search_settings="ontology=ogg&rows=2000&allChildrenOf="
                                                    "http://purl.obolibrary.org/obo/OGG_2000009606")
            else:
                iris = req.get_iris(search_settings="&ontology=ogg&rows=2000")
        else:
            if search_human_only:
                iris = req.get_iris(search_settings="&exact=on&ontology=ogg&rows=2000&"
                                                    "allChildrenOf=http://purl.obolibrary.org/obo/OGG_2000009606")
            else:
                iris = req.get_iris(search_settings="&exact=on&ontology=ogg&rows=2000")

        return iris

    def get_gene_synonyms(self, row):
        """
        collect and clean gene synonyms

        :param list row: row from a df of query input values, containing gene id, gene type, organism id and keywords
        :return: list of gene family roots
        :rtype: list[str]
        """
        search_human_only = False
        fam_roots = []
        all_iris_for_row = []

        # initialise uniprot syns object with DF row containing relevant values, and get the values
        get_uniprot_syns = GetUniprotSynonyms(row, self.logger)
        uniprot_gene_syns, gene_found = get_uniprot_syns.get_gene_syns_from_uniprot()

        if '*' in get_uniprot_syns.gene_id:
            self.logger.warning(f"'{get_uniprot_syns.gene_id}' is a wildcard: searching for a wildcard term may take "
                                f"longer than usual")

        # determine whether to only search under the human gene branch
        if str(get_uniprot_syns.tax_id) == '9606':  # could be str or int, so make str
            search_human_only = True

        # get the ontlogy of genes and genomes syns for the uniprot syns
        for syn in uniprot_gene_syns:
            # get iris
            iris = self.get_gene_syn_iris_from_ebi_ols(syn, search_human_only)

            if iris:
                all_iris_for_row.extend(iris)
            else:
                if not gene_found:
                    self.logger.warning(f"No synonyms were found for '{syn}' in UniProt {row[2]}s or the Ontology of "
                                        f"Genes and Genomes. It is recommended to check for mistakes in this entry")

            # if the gene syn looks like it's from a family, get its root and append to the roots list
            if self.systematic_family_naming_test(syn):
                fam_root = self.get_root_name_from_gene_syn(syn)
                if fam_root:
                    fam_roots.append(fam_root)

        # get the gene and use it to initialise a request object
        orig_gene = get_uniprot_syns.gene_id
        req = OLSRequests(orig_gene, self.logger)

        # get the synonyms for the gene
        syns = self.get_syns_in_parallel(orig_gene, req, all_iris_for_row)

        if len(syns) > 30:
            self.logger.warning(f"{len(syns)} synonyms found for {orig_gene}; you may want to check for noise in the "
                                f"results of the associated search")

        # if the entered gene is a wildcard, add the root phrases for the wildcard's synonyms to the synonym list
        if '*' in orig_gene:
            self.logger.info("Adding synonym roots to the list of gene synonyms")
            roots = self.get_recurrent_root_phrases(syns)
            syns.extend(roots)

        # clean the final lists
        self.logger.info(f"{threading.current_thread().name}: Performing additional cleaning steps specific to genes "
                         f"for '{orig_gene}'")
        final_syns = self.final_gene_clean(orig_gene, syns)
        self.logger.info(f"{threading.current_thread().name}: Done: {len(final_syns)} cleaned synonyms for "
                         f"'{orig_gene}' after additional cleaning")

        if fam_roots:
            self.logger.info(f"{threading.current_thread().name}: One or more synonyms of '{orig_gene}' are predicted "
                             f"to be in a systematically named family of genes")
            root_clean = ExtractOLSSynonyms(orig_gene, self.logger)
            self.logger.info(f"{threading.current_thread().name}: Cleaning collected root phrases for synonyms of "
                             f"'{orig_gene}'")
            fam_roots = root_clean.clean_syn_list(syns_list=fam_roots)
            if orig_gene in fam_roots:
                fam_roots.remove(orig_gene)
            self.logger.info(f"{threading.current_thread().name}: Done: the additional root synonyms {fam_roots} will "
                             f"be used to search for lists containing '{orig_gene}' and its synonyms indirectly")
        return final_syns, fam_roots

    def get_gene_query_string_and_roots(self, row):
        """
        get gene synonyms and make a query string from them. Also get family roots for the gene

        :param list[str] row: row from df of query inputs (index, gene, id type, species, optionally keyword, gene key)
        :return: gene query string, family roots, list of all synonyms searched
        """
        self.logger.info(f"{threading.current_thread().name}: Creating query part for gene '{row[1]}'")

        # get clean lists of synonyms for the supplied gene, and its family root phrases (e.g. ABC is the root for ABC1)
        syn_list, fam_roots = self.get_gene_synonyms(row)
        self.gene_syns = syn_list
        self.gene_syn_roots = fam_roots

        # create the gene query string
        gene_query_string = self.add_abstract_title_and_join_synonyms(syn_list)
        self.logger.info(f"{threading.current_thread().name}: Done: successfully created query part for '{row[1]}'")

        return gene_query_string

    def estimate_query_string_length(self, q_string):
        """
        estimate the length of a query string after URI encoding, with some additional characters for bringing the
         query pieces together
        :param str q_string: the query string
        :return: int
        """
        len_est = 0

        # if the query part exists, then estimate its length
        if q_string:
            encoded_query_string = parse.quote_plus(q_string)
            len_est = len(encoded_query_string) + 12   # + 12 for the surrounding brackets, spaces, ampersand

        return len_est

    def make_len_dict(self, diseases, tissues, others, genes, kwds):
        """
        determine the length of each section of the query, and record it in a dict

        :param str diseases: query string for disease and its synonyms
        :param str tissues: query string for tissue and its synonyms
        :param str others: query string for other settings (e.g. publication year, author etc.)
        :param str genes: query string for genes and their synonyms
        :param list kwds: list of query string for other key words
        :return: dictionary of
                 category: {'string': string, 'len' string length, 'gr_chrs': number of greek characters in the string
        :rtype: dict
        """
        len_dict = {"diseases": {"string": diseases, "len": self.estimate_query_string_length(diseases)},
                    "tissues": {"string": tissues, "len": self.estimate_query_string_length(tissues)},
                    "others": {"string": others, "len": self.estimate_query_string_length(others)},
                    "genes": {"string": genes, "len": self.estimate_query_string_length(genes)}}
        i = 1
        for kwd in kwds:
            len_dict[f"kwd{i}"] = {"string": kwd, "len": self.estimate_query_string_length(kwd)}
            i += 1

        return len_dict

    @staticmethod
    def create_query_length_list(query_len_dict):
        """
        sum the lengths of the splittable strings in the dictionary (splittable strings are any strings except 'others',
        because 'others' can include settings that must be applied to every query)

        :param query_len_dict: {arg: {'string': 'query items', 'len': 11, 'gr_chrs': 5}}
        :return: an estimate for the total length of splittable query strings in the dict
        """
        lens = []
        for key in query_len_dict.keys():
            if key != 'others':
                lens.append(query_len_dict[key]['len'])
        return lens

    def identify_long_queries(self, len_dict, max_length):
        """
        identify the longest query strings that are causing the total query to go over the length limit

        :param dict len_dict: dict containing query strings, their lengths and the number of greek characters in them
        :param int max_length:
        :return: len_dict with longer query strings removed, long_strs dict containing the longer queries to be split
        :rtype: tuple[dict, dict]
        """
        # initialise
        long_strs = {}
        # sum the lengths of the query elements from the dictionary
        len_list = self.create_query_length_list(len_dict)

        while sum(len_list) > (max_length - 500):  # ensures that there are always at least 500 characters for the query
            for key in len_dict.copy().keys():
                # the 'others' query string shouldn't be split up; it contains search settings relevant to every
                # query (its length was checked when created, so should usually be sufficiently short)
                if key != "others":  # the length of 'others' isn't included in the len list
                    # if the query string is the longest, add its relevant info to the long_strs dict
                    if len_dict[key]['len'] == max(len_list):
                        long_strs[key] = {'string': len_dict[key]['string']}
                        # delete the query string's entry in the dictionary
                        del len_dict[key]
            # recreate the len_list from the lengths of the remaining queries in the dictionary
            len_list = self.create_query_length_list(len_dict)

        return len_dict, long_strs

    def split_long_queries(self, len_dict, max_length):
        """
        create multiple smaller queries that cover the full combinations of terms present in supplied query strings

        :param dict len_dict: dict containing query strings, their lengths and the number of greek characters in them
        :param int max_length: 8000 (URI limit) - the length of the constant parts of the URIs (base url, settings)
        :return: len_dict with longer query strings removed, list of lists of combinations of split longer query strings
        :rtype: tuple[dict, list]
        """
        # initialise list of lists
        query_string_chunks = []

        # identify the long query string(s) to be split
        shorter_len_dict, long_strs = self.identify_long_queries(len_dict, max_length)

        # determine the size to split the longer query elements in to:
        # # take away the estimated total length of the shorter queries from max_length,
        #   including an estimate for the number of brackets and ampersands used
        # # divide (to the nearest whole number) by the number of long strings that need to be split
        # # minus 100 more characters to approximately account for the final 'ABSTRACT:"synonym"' at the end of each
        #   chunk (particularly long synonyms may exceed this, but this will be rare)

        est_total_len = sum(self.create_query_length_list(shorter_len_dict)) + (len(shorter_len_dict) * 10)
        chunk_size = ((max_length - est_total_len) // len(long_strs)) - 100

        # for each of the long query strings (often only one)
        for k, v in long_strs.items():
            # initialise variables
            syns = []
            querystring = ""
            querystrings = []
            i = 0
            finished = False

            # get the cleaned synonyms
            if k == 'genes':
                syns = copy.deepcopy(self.gene_syns)
            elif k == 'diseases':
                syns = copy.deepcopy(self.disease_syns)
            elif k == 'tissues':
                syns = copy.deepcopy(self.tissue_syns)
            elif k.startswith("kwd"):
                syns = copy.deepcopy(self.kwd_syn_lists[int(k.strip("kwd")) - 1])
            else:
                self.logger.error(f"Unrecognised query element type '{k}'. Cannot split query.")
                exit() # note: exit commands do not work when running within a thread - which this command often is

            # build the query strings to be no longer than the chunk size
            while syns:
                # calculate the length the query string will be upon the addition of the next synonym, and if the
                # total is lower than the chunk size then incorporate the next synonym
                qstring_len = len(parse.quote_plus(querystring)) + (len(parse.quote_plus(syns[i])) * 3) + 62
                if i == 0 and qstring_len > chunk_size:
                    self.logger.warning(f"The {k} synonym '{syns[i]}' is too long, and will be excluded from the "
                                        f"queries")
                else:
                    while qstring_len <= chunk_size and not finished:
                        i += 1
                        querystring = self.add_abstract_title_and_join_synonyms(syns=syns[:i])
                        try:
                            qstring_len = len(parse.quote_plus(querystring)) + (len(parse.quote_plus(syns[i])) * 3) + 62
                        except IndexError:
                            finished = True  # when i is no longer accessible, every syn is accounted for
                    querystrings.append(querystring)
                    if len(syns) >= i:
                        del syns[:i]
                    else:
                        syns = []
                    querystring = ""
                    i = 0
            # add the query strings for the key to the list of chunks of query strings
            query_string_chunks.append(querystrings)
        # create all combinations of items in all lists in the list of lists query_string_chunks
        chunk_combinations = list(itertools.product(*query_string_chunks))

        return shorter_len_dict, chunk_combinations

    def create_epmc_queries(self, search_description, diseases="", tissues="", others="", genes="", kwds=[]):
        """
        create the EPMC query/queries from the supplied arguments (and their synonyms where relevant)

        :param str search_description: string of gene and search terms
        :param str diseases: string in the format (ABSTRACT:"syn" OR TITLE:"syn" OR ABSTRACT:"syn2"...) for disease syns
        :param str tissues: string in the format (ABSTRACT:"syn" OR TITLE:"syn" OR ABSTRACT:"syn2"...) for tissue syns
        :param str others: string in the format PARAM1:value & PARAM2:value
        :param str genes: string in the format (ABSTRACT:"syn" OR TITLE:"syn" OR ABSTRACT:"syn2"...) for gene syns
        :param list kwds: list of strings in the format (ABSTRACT:"kwd" OR TITLE:"kwd"...) for key words
        :return: list of strings of EPMC queries
        """
        # initialise, calculate the sum of the lengths of query strings for keywords
        queries = []
        self.logger.info(f"{threading.current_thread().name}: Building query")

        len_dict = self.make_len_dict(diseases=diseases, tissues=tissues, genes=genes, others=others, kwds=kwds)

        # the character limit for EPMC API URIs (after encoding) is around 8000, so accounting for constants
        # e.g. some brackets and ampersands, the base url 'https://www.ebi.ac.uk/europepmc/webservices/rest/search'
        # and settings e.g. 'format=xml', use 7500 as the max length to be safe
        leftover_amount = 7500 - len_dict['others']['len']

        query_lengths_list = self.create_query_length_list(len_dict)

        # if the query could be too long, split it in to multiple smaller queries, else make a single query
        # + 20 to approximately account for encoded brackets and ampersands added around/between queries
        if sum(query_lengths_list) + 20 > leftover_amount:
            # log and initialise
            shorter_parts = []
            self.logger.info(f"{threading.current_thread().name}: Query would be too long! Creating multiple shorter "
                             f"queries for the search '{search_description}'")

            # identify and split query parts that are too long
            len_dict, chunk_combinations = self.split_long_queries(len_dict, leftover_amount)
            self.logger.info(f"{threading.current_thread().name}: {len(chunk_combinations)} shorter queries will be "
                             f"built for '{search_description}'")

            # create a partial query from the shorter parts, with other settings last (if present)
            for each in len_dict.keys():
                # don't add empty queries to the final query
                if len_dict[each]["string"] != "" and len_dict[each]["string"] != "()":
                    shorter_parts.append(len_dict[each]["string"])
            consistent_query_string = " & ".join(shorter_parts)

            for combo in chunk_combinations:  # for list in list
                # strip the query strings add brackets
                combo = ["(" + each.strip("(,) ") + ")" for each in combo]
                # join the query string combinations on &
                combo_string = " & ".join(combo)
                # create the full query string and add to the queries list
                queries.append(f"{combo_string} & {consistent_query_string}")

            self.logger.info(f"{threading.current_thread().name}: Done: successfully created {len(queries)} shorter "
                             f"queries to submit for the search '{search_description}'")
        else:
            non_empty = []
            # add non-empty query parts to the list of non-empty query parts
            for each in [genes, diseases, tissues, others]:
                if each != "" and each != "()":
                    non_empty.append(each)
            for each in kwds:
                if each != "" and each != "()":
                    non_empty.append(each)
            # join non-empty query parts on & to form the final query, and add the final query to the queries list
            queries.append(" & ".join(non_empty))

            self.logger.info(f"{threading.current_thread().name}: Done: successfully built Europe PMC query")
        return queries

    def get_roots_and_parts_from_family_synonyms(self):
        """
        create a dictionary of roots of gene names and endings of gene names (e.g. ABC1 --> {ABC: [1]})

        :return: dict of roots and lists of specific parts
        """
        roots_and_parts = {}
        for root in self.gene_syn_roots:
            syns_containing_root = []
            for syn in self.gene_syns:
                # if the root phrase is in the synonym, then remove it from the synonym to leave only the
                # non-root part of the syn
                if root in syn:
                    syns_containing_root.append(syn.replace(root, "").strip())  # non-root part of the gene
            # add the root and the de-duplicated list of syn parts it precedes to the dict of roots and parts
            if syns_containing_root:
                roots_and_parts[root] = list(set(syns_containing_root))

        return roots_and_parts

    def get_query_and_kwd_strings(self, row, diseases=None, tissues=None, others=None, kwds_strings=None, kwds=None):
        """
        create the Europe PMC query/queries for the gene in the row

        :param list[str] row: list of elements of a df row containing gene information and optionally kwd information
        :param str diseases: query string for disease term synonyms
        :param str tissues: query string for tissue term synonyms
        :param str others: query string for other settings
        :param list[str] kwds_strings: list of query strings for keywords
        :param list[str] kwds: keywords to be used in the output
        :return: dict of gene key, search terms, queries and gene family synonyms
        """
        roots_and_parts = {}
        root_queries = []
        all_kwds = []

        # get keywords from the dataframe row if they were not supplied globally at the command line
        if not self.kwds:
            kwds_strings, kwds, all_kwds = self.get_keyword_query_string_and_words(row)

        # create the string of search terms for the relevant output column
        search_term_string = self.get_search_terms_string(kwds)
        search_description = f"{row[1]} and {search_term_string}"

        # create the query string for the gene
        gene_query_string = self.get_gene_query_string_and_roots(row)

        # create the query/queries for the gene and relevant other search terms
        queries = self.create_epmc_queries(search_description, diseases=diseases, tissues=tissues, others=others,
                                           genes=gene_query_string, kwds=kwds_strings)

        # if the gene is in a family, then also query for the root phrases of the family and add these results to a
        # separate dict to be checked for indirect references of the gene in a list (e.g. ABC2 in the list ABC1, 2)
        if self.gene_syn_roots:
            # get a dict of roots of gene synonyms and their counterpart non-root parts
            roots_and_parts = self.get_roots_and_parts_from_family_synonyms()
            if roots_and_parts:
                roots = [f"{x}*" for x in roots_and_parts.keys()]

                gene_fam_query_string = self.add_abstract_title_and_join_synonyms(roots, search_in_kwds=False)
                search_description = f"roots of {search_description}"

                root_queries = self.create_epmc_queries(search_description, diseases=diseases, tissues=tissues,
                                                        others=others, genes=gene_fam_query_string, kwds=kwds_strings)

        specifics = self.gene_syns + self.gene_syn_roots + all_kwds
        # assemble results in to a dictionary
        queries_dict = {'gene key': row[5], 'gene name': row[1], 'search terms': search_term_string, 'queries': queries,
                        'root_syns': roots_and_parts, 'root queries': root_queries, 'query-specific syns': specifics}

        return queries_dict

    @staticmethod
    def extract_results_if_available(result, tag_name):
        """
        tags are not always present and/or populated, so extract their contents conditionally to prevent errors

        :param BeautifulSoup result: soup of a result
        :param str tag_name: name of the tag to find in the soup
        :return: content of the relevant tag, or 'unavailable'
        :rtype: bs4.element.NavigableString or string
        """
        value = "unavailable"
        # if the tag exists
        if result.find(tag_name):
            # try to extract its content
            try:
                value = result.find(tag_name).contents[0]
            # if the tag has no contents, then ignore it
            except IndexError:
                pass
        return value

    def get_relevant_info_from_results(self, result_soup):
        """
        add relevant details about each result to the dictionary of result details

        :param BeautifulSoup result_soup: soup of xml for a result (publication)
        :return: dict of relevant info for the doc
        """
        # initialise
        keywords = []
        preprint_of = "nothing"

        # there should always be an ID and source tag, as these are used by EPMC to create a URL for the document
        pub_id = result_soup.find('id').contents[0]
        url = f"https://europepmc.org/abstract/{result_soup.find('source').contents[0]}/{pub_id}"

        # if a standard tag is present, assign its value to the variable. else, assign "unavailable"
        title = self.extract_results_if_available(result_soup, 'title')
        year = self.extract_results_if_available(result_soup, 'pubyear')
        authorstring = self.extract_results_if_available(result_soup, 'authorstring')
        abstract = self.extract_results_if_available(result_soup, 'abstracttext')
        p_types = self.extract_results_if_available(result_soup, 'pubtype')

        # add each keyword to the list of keywords (or leave the list empty if there are no keyword tags in the doc)
        for res in result_soup.find_all('keyword'):
            try:
                keywords.append(res.contents[0])
            # handle empty keyword tags (i.e. ignore them)
            except IndexError:
                pass

        if pub_id.startswith("PPR"):
            comment_list = result_soup.find_all('commentcorrection')
            for comment in comment_list:
                comment_type = comment.find('type')
                if comment_type.contents[0] == "Preprint of":
                    preprint_of = comment.find('id').contents[0]

        return pub_id, {"ID": pub_id, "Title": title, "Year": year, "Author": authorstring, "Publication type": p_types,
                        "Abstract": abstract, "Keywords": keywords, "url": url, "prep_of": preprint_of}

    def parse_potential_result_content_for_relevant_info(self, res, ignore_ids):
        """
        result contents can be large, so perform initial parsing to collect just the relevant information

        :param str res: string of XML result page
        :param list ignore_ids: list of ids to ignore
        :return: list of dicts of relevant info for each document,
        """
        doc_info_list = []
        cursor_mark = "*"

        # get soup
        whole_results_soup = BeautifulSoup(res, features='lxml')

        # get the next page cursor code from the soup
        # EDIT: now using a try, because the tag is no longer present in the XML if there is not a next page
        try:
            cursor_mark = whole_results_soup.find("nextcursormark").contents[0]
        except AttributeError:
            self.logger.info(f"No next cursor mark element in results")

        # get the hit count
        hit_count = whole_results_soup.find("hitcount").contents[0]

        # for each result listed in the page, get its important information
        result_tags = whole_results_soup.find_all("result")
        if result_tags:
            for result_tag in result_tags:
                if result_tag.id.contents[0] not in ignore_ids:
                    doc_id, doc_info_dict = self.get_relevant_info_from_results(result_tag)
                    doc_info_list.append(doc_info_dict)

        unique_docs = list({v['ID']: v for v in doc_info_list}.values())

        return unique_docs, cursor_mark, hit_count

    def parse_actual_result_content_for_relevant_info(self, result_text_list):
        """
        result contents can be large, so perform initial parsing to collect just the relevant information

        :param list result_text_list: list of strings of XML result pages
        :return: string of query strings used, list of dicts of relevant info for each document, list of doc IDs found
        :rtype: tuple(str, list, list)
        """
        string_number = 0
        query_strings = []
        doc_info_list = []
        id_list = []

        for res in result_text_list:
            # get soup
            whole_results_soup = BeautifulSoup(res, features='lxml')

            # get the number of results and query string for the results page, log them, and add the string and its
            # index number to the list of query strings to be printed in the output files
            result_count = whole_results_soup.hitcount.contents[0]
            query_string = whole_results_soup.querystring.contents[0]
            self.logger.info(f"{threading.current_thread().name}: {result_count} results found for the query string "
                             f"'{query_string}'")
            string_number += 1
            query_strings.append(f"{string_number}: {query_string}")

            # for each result listed in the page, get its important information
            result_tags = whole_results_soup.find_all("result")
            if result_tags:
                for result_tag in result_tags:
                    doc_id, doc_info_dict = self.get_relevant_info_from_results(result_tag)
                    doc_info_list.append(doc_info_dict)
                    id_list.append(str(doc_id))
            else:
                # create an empty dict of relevant headings if there were no results returned
                doc_info_list.append({"ID": "none", "Title": "", "Year": "", "Author": "", "Publication type": "",
                                      "Abstract": "", "Keywords": "", "url": "", "prep_of": ""})

        # make unique
        unique_docs = list({v['ID']: v for v in doc_info_list}.values())
        unique_ids = list(set(id_list))

        # make single string
        string_of_queries = " split_here ".join(query_strings)

        return string_of_queries, unique_docs, unique_ids

    def run_epmc_query(self, query):
        """
        run the supplied EPMC query and return the first page of results (up to 1000 results; no need to collect more)

        :param str query: the query string to search for
        :return: EPMC query response
        :raises HTTPError: if a HTTP error occurs
        :raises ConnectionError: if a connection error occurs
        """
        res_text = "No result text"

        # run the query and collect its content, or return appropriate errors
        try:
            self.logger.info(f"{threading.current_thread().name} Querying Europe PMC for {query}")
            res = requests.get('https://www.ebi.ac.uk/europepmc/webservices/rest/search',
                               headers={'Content-type': "application/x-www-form-urlencoded"},
                               params={'query': query,
                                       'resultType': 'core',
                                       'pageSize': 1000,
                                       'format': 'xml'}
                               )
            res.raise_for_status()
            res_text = res.content  # use content rather than text to ensure proper encoding of greek characters etc
            # log status code of response and close connection
            self.logger.info(f"{threading.current_thread().name}: {res}")
            res.close()

        # raise relevant exceptions if there are errors for the request
        except (requests.exceptions.HTTPError, requests.exceptions.ConnectionError) as err:
            self.logger.error(err)
            exit() # note: exit commands do not work when running within a thread - which this command often is

        return res_text

    def analyse_text_for_gene_in_list(self, abstract, title, root_and_syn_part, doc_id):
        """
        check the abstract and title for lists that indirectly contain the gene of interest (e.g. the list
        "ABC1, 2 and 3" contains ABC2 and ABC3 indirectly). return true if either contains a relevant list

        :param bs4.element.NavigableString abstract: abstract text (extracted from soup)
        :param bs4.element.NavigableString title: title text (extracted from soup)
        :param dict root_and_syn_part: root phrases are keys, lists of synonyms with root phrases removed are values
        :param str doc_id: document ID
        :return: bool: True if the abstract contains the gene of interest
        """
        text = ". ".join([str(title), str(abstract)]).lower()

        # if the text contains a root and any associated synonym part, then check if the root and synonym are in a
        # list-type structure
        for root, syn_part_list in root_and_syn_part.items():
            if root.lower() in text:
                for syn_part in syn_part_list:
                    if syn_part.lower() in text:
                        # if the text contains the structure "root phrase optionally followed by ', s and/or
                        # spaces, optionally followed by any number of comma-separated numbers with optional
                        # letters (e.g. '1, 1A1'), non-optionally followed by [and , or] and the relevant part
                        # of the synonym", then return True
                        if re.match(fr".*{root.lower()}'?({'|'.join(hyphens)}|\s)*\d+[a-z]?\d*({'|'.join(hyphens)}|\s)*"
                                    fr"(,({'|'.join(hyphens)}|\s)*\d+[a-z]?\d*)*"
                                    fr"(and|or|,)\s({'|'.join(hyphens)}|\s)*{syn_part.lower()}.*", text):
                            self.logger.info(
                                f"Document {doc_id} contains the synonym '{root}{syn_part}' indirectly in a list")
                            return True
        return False

    def run_epmc_query_all_results(self, query, res_count, id_list, root_and_syn_parts, gene_name):
        """
        run the supplied EPMC query and return all results, even if split over multiple pages

        :param str query: the query string to search for
        :param int res_count: the number of actual results
        :param list id_list: the list of IDs of actual results:
        :param dict root_and_syn_parts: dict containing syn roots and their final parts (e.g. {'ABC': ['3']}
        :param str gene_name: gene name
        :return: list of dicts of relevant information from relevant documents
        :raises HTTPError: if a HTTP error occurs
        :raises ConnectionError: if a connection error occurs
        """
        res_info = []
        c_mark = '*'
        run_number = 0
        continue_querying = True

        while continue_querying and res_count <= 1000:
            run_number += 1
            res_text = "No result text"

            # run the query and collect its content, or return appropriate errors
            try:
                self.logger.info(f"{threading.current_thread().name} Querying Europe PMC with root query {run_number}: "
                                 f"{query}")
                res = requests.get('https://www.ebi.ac.uk/europepmc/webservices/rest/search',
                                   headers={'Content-type': "application/x-www-form-urlencoded"},
                                   params={'query': query,
                                           'resultType': 'core',
                                           'pageSize': 1000,
                                           'format': 'xml',
                                           'cursorMark': c_mark}
                                   )
                res.raise_for_status()
                res_text = res.content  # use content rather than text to ensure proper encoding of greek characters etc
                # log status code of response and close connection
                self.logger.info(f"{threading.current_thread().name}: {res}")
                res.close()

            # raise relevant exceptions if there are errors for the request
            except (requests.exceptions.HTTPError, requests.exceptions.ConnectionError) as err:
                self.logger.error(err)
                exit() # note: exit commands do not work when running within a thread - which this command often is

            potential_doc_info, c_mark, hit_count = \
                self.parse_potential_result_content_for_relevant_info(res_text, id_list)

            if int(hit_count) > 5000 and run_number == 1:
                self.logger.warning(f"{hit_count} results were identified when searching for a root of a synonym of "
                                    f"{gene_name}. It may take several minutes to parse this many documents")

            if root_and_syn_parts:
                for doc in potential_doc_info:
                    abstract = doc['Abstract']
                    title = doc['Title']
                    doc_id = doc['ID']
                    correct = self.analyse_text_for_gene_in_list(abstract, title, root_and_syn_parts, doc_id)
                    # if the root is followed by the value in a list, then add its info to the list of actual results
                    if correct:
                        res_info.append(doc)
                        res_count += 1

            # if total hits - max hits per page processed <= 0 (i.e. if all hits processed), stop querying
            if int(hit_count) - (run_number * 1000) <= 0:
                continue_querying = False
        return res_info

    def run_queries_and_get_results(self, queries_dict):
        """
        run the generated Europe PMC queries and add their results to the dictionary

        :param queries_dict: dict of gene key, search terms, queries and family synonyms
        :return: dictionary of gene keys and relevant information about query results for the output
        """
        # initialise job handler and potential results list
        ReadyThready.set_logger(self.logger)
        root_doc_info = []

        # extract the lists of queries from the dict
        queries = queries_dict['queries']
        root_qs = queries_dict['root queries']
        root_and_syn_parts = queries_dict['root_syns']

        self.logger.info(f"Running {len(queries)} queries in Europe PMC for '{queries_dict['gene name']}'")
        result_text = ReadyThready.go(self.run_epmc_query, [queries], 0, self.n_threads)
        self.logger.info(f"Done running '{queries_dict['gene name']}' queries in Europe PMC")
        query_strings, doc_info, id_list = self.parse_actual_result_content_for_relevant_info(result_text)

        res_count = len(id_list)

        if root_qs and res_count < 1000:
            self.logger.info(f"Running {len(root_qs)} queries in Europe PMC for "
                             f"roots of '{queries_dict['gene name']}'")
            root_results = ReadyThready.go(self.run_epmc_query_all_results,
                                           [root_qs, res_count, id_list, root_and_syn_parts, queries_dict['gene name']],
                                           0, self.n_threads)
            # flatten list (list of items expected, not list of lists)
            self.logger.info(f"Done running '{queries_dict['gene name']}' root "
                             f"queries in Europe PMC")

            root_doc_info = [item for sublist in root_results for item in sublist]
            # print(len(root_doc_info))
            doc_info.extend(root_doc_info)

        # add the list of results text strings to the dictionary
        queries_dict["query_strings"] = query_strings
        queries_dict["true_result_doc_info"] = doc_info
        queries_dict["res_count"] = res_count + len(root_doc_info)

        # remove the queries from the dict as they are no longer needed
        del queries_dict['root queries']
        del queries_dict['queries']
        del queries_dict['root_syns']

        return queries_dict

    def get_results_from_epmc(self):
        """
        create and run queries using the supplied terms and their synonyms, and parse the relevant information from the
        returned results

        :return: list of dicts of information about the queries and results for each input
        :rtype: list[dict]
        """
        ReadyThready.set_logger(self.logger)
        self.logger.info(f"{threading.current_thread().name}: Creating and running Europe PMC queries")

        # create query strings for universal parameters
        dis_q_string, tis_q_string, others_q_string, kwd_q_strings, kwds = self.get_constant_query_strings()

        # create list of lists of input information for multithreading
        row_list = [[row[0], row[1], row[2], row[3], row[4], row[5]] for row in self.df.itertuples()]

        # create queries using multiple threads
        queries = ReadyThready.go(func=self.get_query_and_kwd_strings, arg_data_index=0, n_threads=self.n_threads,
                                  args=[row_list, dis_q_string, tis_q_string, others_q_string, kwd_q_strings, kwds])
        self.logger.info(f"Done: finished creating queries")

        # run the created queries using multiple threads
        self.logger.info(f"Running queries")
        results = ReadyThready.go(func=self.run_queries_and_get_results, args=[queries], arg_data_index=0,
                                  n_threads=self.n_threads)
        self.logger.info(f"Done: successfully finished running all queries")

        return results
