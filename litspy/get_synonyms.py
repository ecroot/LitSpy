import requests
import json
import re
from time import sleep
import threading
from rtgo import ReadyThready

import litspy.alternative_characters as chars
import litspy.noisy_phrases as noise
from litspy.anatomy_qualifiers import anatomy_qualifiers as anatomy_qualifiers


class OLSRequests:
    """class for making requests to the EBI OLS"""

    def __init__(self, original_term, logger):
        """
        initialise with logger and the argument supplied at the command line

        :param str original_term: a supplied gene, disease or tissue/organ
        :param logging.Logger logger: project logger
        """
        self.logger = logger
        self.original_term = original_term

    def get_request_parse_result(self, url):
        """
        get request for the url, or raise relevant exceptions and parse the result in to a dict

        :param str url: url
        :return: information parsed from the request result
        :rtype: dict
        :raises HTTPError: if a HTTP error occurs
        :raises ConnectionError: if a connection error occurs
        """
        self.logger.info(f"{threading.current_thread().name}: Requesting {url}")
        # make the request
        try:
            res = requests.get(url)
            res.raise_for_status()
        # raise relevant exceptions if there are errors for the request
        except (requests.exceptions.HTTPError, requests.exceptions.ConnectionError) as err:
            self.logger.error(err)
            exit() # note: exit commands do not work when running within a thread - which this command often is
        # log status code of response
        self.logger.info(f"{threading.current_thread().name}: {res} for {url}")
        # parse the test of the response in to a dict
        parsed_json = json.loads(res.text)
        # close the connection
        res.close()
        # return parsed JSON
        return parsed_json

    def get_iris(self, search_settings=None):
        """
        search for the term in EBI OLS, return a IRI (Internationalized Resource Identifier) for each search result

        :param str search_settings: optional string of search settings, e.g. ontology=ontology_name
        :return: list of unique IRIs
        """
        # initialise the iri list, and prevent performing too many requests if this method is called in a loop
        iris = []
        sleep(0.1)

        # replace any spaces in the term with +, for the url
        term = self.original_term.replace(" ", "+")

        # if exact mode is not on, then wrap term in quote marks to ensure the whole phrase is searched as one
        if 'exact=on' not in search_settings:
            term = f'"{term}"'

        # query EBI OLS for the term, with settings if supplied
        if search_settings:
            url = f"http://www.ebi.ac.uk/ols/api/search?q={term}{search_settings}"
        else:
            url = f"http://www.ebi.ac.uk/ols/api/search?q={term}"
        parsed_json = self.get_request_parse_result(url)

        # extract IRIs and append them to the iri list
        for result in parsed_json["response"]["docs"]:
            iri = result["iri"]
            if "/obo/BFO_0" in iri:
                pass
            else:
                iris.append(iri)

        # de-duplicate list, log and return
        iris = list(set(iris))
        self.logger.info(f"{threading.current_thread().name}: Found {len(iris)} IRIs for ontology nodes relevant to "
                         f"'{term}'")
        # return unique IRIs
        return iris

    def get_json_for_iri(self, iri):
        """
        query EBI OLS for the IRI, return parsed response text

        :param str iri: IRI for a relevant ontology node
        :return: dict response text parsed in to a dict
        :rtype: dict
        """
        # prevent performing too many requests if this method is called in a loop
        sleep(0.1)
        self.logger.info(f"{threading.current_thread().name}: Querying {iri} for synonyms of '{self.original_term}'")
        # get response for the IRI
        url = f"http://www.ebi.ac.uk/ols/api/terms?iri={iri}&size=1000"
        parsed_json = self.get_request_parse_result(url)
        # return parsed JSON
        return parsed_json

    def get_json_for_full_url(self, url, description):
        """
        get a full the URL, parse the response in to a dict

        :param str url: a full URL for a relevant page, e.g. hierarchical descendants for a term
        :param str description: description of the expected URL contents
        :return: parsed JSON of the response text
        :rtype: dict
        """
        self.logger.info(f"{threading.current_thread().name}: Querying {url} for {description}")
        parsed_json = self.get_request_parse_result(url)
        return parsed_json


class ExtractOLSSynonyms:
    """class for extracting and cleaning synonyms from OLS's json-derived dicts"""
    def __init__(self, original_term, logger, min_syn_len, n_threads=0):
        """
        initialise synonyms object with logger, original term, synonym list, headers, noise greek letters

        :param str original_term: a supplied gene, disease or tissue/organ
        :param logging.Logger logger: project logger
        :param int min_syn_len: synonyms shorter than this value will be discarded
        :param int n_threads: number of threads to use in any multithreading tasks
        """
        self.logger = logger
        self.original_term = original_term
        self.min_syn_len = min_syn_len
        self.n_threads = n_threads
        self.syns = []

        # EBI OLS annotation headers that can contain synonyms
        self.relevant_keys = ["has_related_synonym", "alternative term", "comment", "description",
                              "symbol from nomenclature authority", "hasExactSynonym"]

        # create a list of all the greek characters that are handled in the alternative characters module
        greek_char_names = list(chars.greek_dict.keys())
        [greek_char_names.extend(each) for each in list(chars.greek_dict.values())]
        self.greek_char_list = greek_char_names

    def extract_syns(self, term):
        """
        extract synonyms from relevant headers

        :param dict term: the 'term' part of a full page's parsed JSON
        :return: append synonyms to the synonym list attribute
        :rtype: None
        """
        self.logger.info(f"{threading.current_thread().name}: extracting synonyms from element")

        # if exists and not null, capture synonyms and add them to the synonym list
        if "synonyms" in term and term["synonyms"]:  # if exists and not null
            for syn in term["synonyms"]:
                self.syns.append(syn)

        # if exists and not null, capture label and add to the synonym list
        if "label" in term and term["label"]:  # if exists and not null
            self.syns.append(term["label"])

        # if the term has annotations and there are any relevant headers within the annotations, then capture
        # synonyms from the relevant header(s) and add them to the synonym list
        if "annotation" in term:
            for relevant_key in self.relevant_keys:
                if relevant_key in term["annotation"] and term["annotation"][relevant_key]:  # if exists and not null
                    for syn in term["annotation"][relevant_key]:  # values of children of annotations are lists
                        # get synonyms from a pipe-separated list of "other designations"
                        if syn.startswith("Other designations:"):
                            syn = syn.replace("Other designations:", "")
                            other_designations = syn.split("|")
                            for o_des in other_designations:
                                self.syns.append(o_des.strip())
                        else:
                            self.syns.append(syn)

        # de-duplicate as we go to prevent the list getting too long if terms have many syns/many ontology nodes
        self.syns = list(set(self.syns))

    def get_elems_and_log(self, parsed_json):
        """
        get a count of elements from the json object and log it appropriately

        :param dict parsed_json: a parsed JSON object for an EBI OLS page
        :return: elems
        :rtype: int or None
        """
        elems = None
        # if the expected standard keys exist, obtain the count of elements
        if "page" in parsed_json and 'totalElements' in parsed_json['page']:
            elems = parsed_json['page']['totalElements']
            # log relevant messages about the number of elements found, only once (for the first page)
            if 'totalElements' in parsed_json['page'] and int(parsed_json['page']['number']) == 0:
                if elems > 50:
                    self.logger.warning(f"{elems} elements were found within one of the synonym search results for "
                                        f"'{self.original_term}'. This may cause a longer than usual running time")
                else:
                    self.logger.info(f"{elems} elements were found within one of the synonym search results for "
                                     f"'{self.original_term}'")
        return elems

    def get_syns_from_next_page(self, parsed_json):
        """
        get the url for the next page and get synonyms for it

        :param dict parsed_json: a parsed JSON object for an EBI OLS page
        :return:
        """
        # short sleep to prevent sending too many requests
        sleep(0.05)

        # extract URL string for next page
        next_page_url = parsed_json["_links"]["next"]["href"]

        # initialise request object with term and logger
        next_page_request = OLSRequests(original_term=self.original_term, logger=self.logger)

        # get and parse the next page
        parsed_next_page = next_page_request.get_json_for_full_url(
            url=next_page_url, description=f"{self.original_term} [next page]"
        )

        # get synonyms from the parsed next page
        self.get_syns_from_json(parsed_next_page)

    def get_syns_from_json(self, parsed_json, descendants=False):
        """
        access relevant values within the dict of the parsed JSON, extract synonyms from the page, and move on to the
        next page to do the same if one exists

        :param dict parsed_json: a parsed JSON object for an EBI OLS page
        :param bool descendants: whether getting descendants
        :return: add synonyms to the synonym list attribute
        :rtype: None
        """
        desc = ""
        # initialise job handler object
        ReadyThready.set_logger(self.logger)
        if descendants:
            desc = " descendants of"
        self.logger.info(f"{threading.current_thread().name}: Collecting synonyms from page for{desc} '"
                         f"{self.original_term}'")

        elems = self.get_elems_and_log(parsed_json)

        # if expected keys exist and are not null, extract synonyms from relevant fields
        if "_embedded" in parsed_json and "terms" in parsed_json["_embedded"]:
            if elems > 1:
                ReadyThready.go(self.extract_syns, [parsed_json["_embedded"]["terms"]], 0, self.n_threads)
            else:
                for term in parsed_json["_embedded"]["terms"]:
                    self.extract_syns(term)
            self.logger.info(f"{threading.current_thread().name}: Done: successfully collected synonyms of{desc} "
                             f"'{self.original_term}' from page")
        else:
            # if there are no '_embedded' and 'terms' keys because there are no results, then log this to info
            if elems == 0:
                self.logger.info(f"{threading.current_thread().name}: No synonyms found for '{self.original_term}'")
            # if there are different unexpected keys/headers, log a warning containing the unexpected json as a dict
            else:
                self.logger.warning(f"{threading.current_thread().name}: Unexpected values in EBI OLS json object when "
                                    f"searching for '{self.original_term}'. No synonyms obtained from the following: "
                                    f"{parsed_json}")

        # if the results page is part of a set and there is a 'next' page, then also get the synonyms from the next page
        if "_links" in parsed_json and "next" in parsed_json["_links"]:
            self.get_syns_from_next_page(parsed_json)

    def get_syns_of_descendants(self, parsed_json):
        """
        get the descendants for all search results using the specialised EBI OLS API URL

        :param dict parsed_json: request response's text information parsed in to a dict
        :return: append synonyms to the synonym list
        """
        # initialise OLS request object
        descendant_request = OLSRequests(original_term=self.original_term, logger=self.logger)
        ids = []
        description = f"descendants of '{self.original_term}'"
        self.logger.info(f"{threading.current_thread().name}: Collecting synonyms for {description}")
        # get the ID for each of the query result terms
        if "_embedded" in parsed_json and "terms" in parsed_json["_embedded"]:
            for term in parsed_json["_embedded"]["terms"]:
                ids.append(term["obo_id"])
        else:
            self.logger.warning(f"{threading.current_thread().name}: Unexpected values in EBI OLS json object for a "
                                f"descendant of search term {self.original_term}. No synonyms will be derived from: "
                                f"'{parsed_json}'")
        # get unique IDs
        ids = set(ids)
        # for each id, get synonyms of its hierarchical descendants
        for obo_id in ids:
            # get parsed JSON of descendants
            parsed_descendants = descendant_request.get_json_for_full_url(
                f'http://www.ebi.ac.uk/ols/api/ontologies/uberon/hierarchicalDescendants?id={obo_id}&size=1000',
                description
            )
            # get synonyms from parsed JSON of descendants
            self.get_syns_from_json(parsed_descendants, descendants=True)

    @staticmethod
    def word_search(word):
        """
        determine whether a word/phrase is a boundaried substring of another phrase

        :param str word: a word or phrase that is a synonym of the original term, and is not known to be redundant
        :return: compiled regex with escaped special characters and boundaries placed around the word/phrase
        """
        # escape special characters
        word = re.escape(word)
        return re.compile(fr'\b({word})\b', flags=re.IGNORECASE).search

    def remove_noise_and_punctuation(self, syn_list, syn_type):
        """
        initial cleaning of the synonym list

        :param list syn_list: list of synonyms
        :param str syn_type: the type of synonym being cleaned
        :return: partially cleaned list of synonyms
        :rtype: list
        """
        # initialise clean list, de-duplicate
        filtered_terms = []
        unique = set(syn_list)

        # clean each term
        for term in unique:
            # ignore terms that contain the original term within them
            if self.word_search(self.original_term)(term):
                pass
            # ignore terms that contain noise indicators, e.g. "Editor note"
            elif not term.startswith("GO:") and \
                    any(noisy_term.upper() in term.upper() for noisy_term in noise.synonym_noise_indicators):
                pass
            # ignore terms that contain . unless there are numbers in the term
            elif "." in term and not re.search(r"\d+", term):
                pass
            else:
                # replace terms and punctuation with spaces (some syns use punctuation instead of spaces)
                term = term.replace("EXACT", " ")  # remove "EXACT" (commonly added to the end of synonyms)
                term = term.replace("susceptibility to", " ")  # remove this common phrase
                term = term.replace("working designation", " ")  # remove this common phrase
                term = term.replace("_", " ")  # replace underscores with spaces
                term = term.replace(",", " ")  # EPMC handles commas/no commas
                term = term.replace("?", " ")  # common in some ontologies, possibly in place of hyphen characters
                term = term.replace("\"", " ")  # remove quote marks
                # remove fake quote marks
                term = term.replace("“", " ")
                term = term.replace("”", " ")
                for char in chars.hyphens:
                    term = term.replace(char, " ")  # EPMC treats spaces/hyphens/dashes the same, so normalise to spaces
                # remove brackets and their contents, unless the contents are digits or numerals
                if not re.findall(fr"\(.*[{''.join(chars.numerals)}\d]+.*\)", term):
                    term = re.sub(r"[\(\[]+.*?[\)\]]+", "", term)
                # EPMC handles brackets/no brackets, so remove to reduce character count
                term = term.replace("(", " ")
                term = term.replace(")", " ")
                term = re.sub(r"\n", " ", term)  # remove any newline characters
                term = re.sub(r"\s{2,}", " ", term)  # de-duplicate spaces
                term = term.strip()  # strip trailing spaces

                # keep cleaned terms that are at least the minimum syn length
                if len(term) >= self.min_syn_len:
                    # for tissues, remove syns that are a single letter followed by digits
                    # (e.g. A10 is very noisy, but only remove these syns from tissues because e.g. p53 is a valid gene)
                    if syn_type == "tissue":
                        if not re.match(r"[A-Z]\d+$", term):
                            filtered_terms.append(term)
                    else:
                        filtered_terms.append(term)

        # add the original term to the cleaned syn list, and return the list
        filtered_terms.insert(0, self.original_term)
        return filtered_terms

    def remove_redundant_syns(self, syns):
        """
        identify and remove synonyms that contain a shorter synonym within them
        A synonym is considered redundant if another synonym in the list is a substring of it, because querying for
        the shorter substring will also return all instances of the longer substring

        :param list syns: list of synonyms with noise filtered out
        :return: list of non-redundant synonyms
        """
        redundant = []

        # sort synonym list on ascending number of words in synonym
        sorted_on_phrase_length = [syn.split() for syn in syns]
        sorted_on_phrase_length.sort(key=len)
        sorted_on_phrase_length = [" ".join(syn) for syn in sorted_on_phrase_length]

        # determine redundant synonyms:
        # starting from the phrase with the smallest number of words, if the phrase is not known to be redundant,
        # then check whether the phrase is a substring of any of the other phrases
        for term in sorted_on_phrase_length:
            if term not in redundant:
                for syn in syns:
                    if syn == term:
                        pass  # ignore self
                    # if the term is a substring of the syn, add the syn to the list of redundant phrases
                    elif self.word_search(term)(syn):
                        redundant.append(syn)

        # remove redundant synonyms from the list of synonyms
        for redundant_syn in redundant:
            if redundant_syn in syns:
                syns.remove(redundant_syn)
        return syns

    @staticmethod
    def add_space_before_number(syns):
        """
        EPMC treats spaces/hyphens/dashes the same, but does not treat spaces & no spaces the same:
        ADAMTS-5 == ADAMTS 5, but ADAMTS 5 =/= ADAMTS5
        So, for any instances of synonyms without spaces preceding numbers (except C11orf22 type synonyms), also add
        a version with a space (e.g. for "ADAMTS5", also add "ADAMTS 5" to the list)

        :param list syns: list of synonyms
        :return: list of synonyms
        """
        diff_spacing = []
        
        # get numbers (and optionally capital letters) that are not preceded by hyphens, numbers or spaces
        number_in_gene = re.compile(fr".*[^\s\d]+\d+\s*[A-z]?[\s\d]*\b")

        for syn in syns:
            # replace any hyphens with spaces for easier handling
            # (in other cleaning steps, hyphens are removed from synonyms but not from original terms)
            for hyp in chars.hyphens:
                syn = re.sub(hyp, " ", syn)
            # for synonyms except those in formats such as orf, KIAA, UNQ codes, p53/A4, or full phrases
            if not re.match(r"[Cc]\d+orf\d+", syn) and not re.match(r"UNQ\d+/PRO\d+", syn) and not \
                    re.match(r"KIAA\d+", syn) and not re.fullmatch(r"[A-z]\d{1,3}", syn) \
                    and not re.fullmatch(r"(CI|RR|AIM|FBS|IOP)\d+", syn) and not len(syn.split(" ")) > 2:
                for syn_match in re.finditer(number_in_gene, syn):
                    res = syn_match.group()
                    n = 0
                    for num_match in re.finditer(r"\d+", res):
                        n += 1
                        match = num_match.group()
                        space_variant = (re.sub(match, " " + match, syn, n))
                        diff_spacing.append(space_variant)
                        position = n
                        while re.search(number_in_gene, space_variant):
                            # create other spacing variants for the term (e.g. if the term contains >1 number)
                            position += 1
                            space_variant = (re.sub(r"(\d+)", r" \1", space_variant, position))
                            diff_spacing.append(re.sub(r"\s+", " ", space_variant))

        # add the spacing variant synonyms to the synonym list, de-duplicate and return
        syns.extend(diff_spacing)
        return list(set(syns))

    @staticmethod
    def remove_space_hyphen_before_number(syns):
        """
        EPMC treats spaces/hyphens/dashes the same, but does not treat spaces & no spaces the same:
        ADAMTS-5 == ADAMTS 5, but ADAMTS 5 =/= ADAMTS5
        so, for any instances of a number preceded by a space/hyphen/dash, also add a version without the space/hyphen/
        dash to the list - e.g. for "ADAMTS-5" also add "ADAMTS5" to the list

        :param list syns: list of synonyms
        :return: list of synonyms
        """
        diff_spacing = []
        # regex for one or more numbers, plus optional capital letters, preceded by a space
        number_in_gene = re.compile(fr"\s+\d+\s*[A-z]?[\s\d]*\b")

        for syn in syns:
            n = 0
            # replace any hyphens with spaces for easier handling
            # (in other cleaning steps, hyphens are removed from synonyms but not from original terms)
            for hyp in chars.hyphens:
                syn = re.sub(hyp, " ", syn)
            for match in re.finditer(number_in_gene, syn):
                n += 1
                res = match.group(0)
                diff_spacing.append(re.sub(res, res.lstrip(), syn, n))

        # add the spacing variants to the list of synonyms, de-duplicate and return
        syns.extend(diff_spacing)
        return list(set(syns))

    @staticmethod
    def remove_multiple_spaces(syns):
        """
        following creation of space variants, remove multiple spaces from terms

        :param list syns: list of synonyms
        :return: list of synonyms
        """
        cleaned_syns = []
        for syn in syns:
            cleaned_syns.append(re.sub(r"\s{2,}", " ", syn))  # replace all instances of 2 or more spaces with 1 space
        return cleaned_syns

    def greek_char_and_spacing_expansion(self, syns):
        """
        for more complete literature search, account for variations in e.g. hyphenation, numerals in synonyms

        :param list syns: list of synonyms
        :return: list of synonyms with expanded variations of punctuation
        """
        # expand greek letters to include variants with words, chars and chars commonly used in place of true chars
        syns = self.expand_greek_letters(syns)

        # add and remove spaces: in EPMC searching, 'ABC-1' and 'ABC 1' are equivalent, but 'ABC1' and 'ABC 1' are not
        self.logger.debug("Creating spacing variants between letters and numbers in synonyms")
        syns = self.add_space_before_number(syns)
        syns = self.remove_space_hyphen_before_number(syns)
        syns = self.remove_multiple_spaces(syns)
        return syns

    def create_type_variations(self, type_phrase, syn):
        """
        create variations for synonyms containing types (e.g. "type 1 collagen" and "collagen type 1").
        Add the variations to the list of synonyms (EPMC handles variations in hyphenation, commas and brackets, so
        there is no need to account for these in the generated variations)
        :param type_phrase: string - the phrase concerning the type (e.g. "a-type", "type 1a")
        :param syn: string - the full synonym
        :return: list of synonyms
        """
        syn_list = []
        type_phrase = type_phrase.strip()
        intermediate = re.sub(type_phrase, "", syn, flags=re.IGNORECASE).strip()
        end = intermediate + " " + type_phrase
        start = type_phrase + " " + intermediate
        for each in end, start:
            syn_list.append(each)
        syn_list = self.remove_multiple_spaces(syn_list)
        return syn_list

    def expand_types(self, syns):
        """
        if the word 'type' is found in a synonym followed by a number/numeral/greek letter, then ensure that other
        variants of the phrase are added (e.g. "collagen (type 1)" should be expanded to "collagen 1", "type 1 collagen"

        :param list syns: list of synonyms
        :return: list of synonyms
        """
        # log and initialise
        self.logger.debug("Creating variants of synonyms that contain 'type' to include multiple phrase orders")
        additional_syns = []

        for syn in syns:
            if "type" in syn.lower():  # if the synonym contains 'type'
                # create compiled regexes for finding phrases containing types
                x_hyphen_type = re.compile(fr"(.+\s*[{''.join(chars.hyphens)}]\s*[Tt][Yy][Pp][Ee])")
                type_x = re.compile(fr"^.*([Tt][Yy][Pp][Ee][\s\d{''.join(chars.numerals)}]*\b[a-zA-Z]?\b[\s\d"
                                    fr"{''.join(chars.numerals)}]*({'|'.join(self.greek_char_list)})*[\s\d]*)")

                hyphen_result = re.search(x_hyphen_type, syn.lower())
                if hyphen_result:  # unlikely; most hyphens have been removed by this stage of processing
                    # create variations of the synonym and add them to additional_syns
                    additional_syns.extend(self.create_type_variations(hyphen_result.group(1), syn))

                type_x_result = re.search(type_x, syn)
                # if type is followed by expected characters such as numbers, numerals, greek characters, get variations
                if type_x_result and type_x_result.group(1).strip().lower() != "type":
                    # create variations of the synonym and add them to additional_syns
                    additional_syns.extend(self.create_type_variations(type_x_result.group(1), syn))

        # add the additional synonyms to the list of synonyms
        syns.extend(additional_syns)
        return list(set(syns))

    def remove_chain_from_end(self, syns):
        """
        for syns that end in "chain", remove "chain"

        :param list syns: list of synonyms
        :return: list of synonyms
        """
        # log and initialise
        self.logger.debug("Stripping 'chain' from the end of synonyms")
        new_syns = []

        # for every synonym in the list, remove trailing spaces, remove 'chain'/'chains', remove trailing spaces again
        for syn in syns:
            syn = syn.strip()
            if syn.endswith("chain"):
                syn = re.sub("chain", "", syn, flags=re.IGNORECASE).strip()
            elif syn.endswith("chains"):
                syn = re.sub("chains", "", syn, flags=re.IGNORECASE).strip()
            new_syns.append(syn)

        # return de-duplicated synonyms with chain/chains removed
        return list(set(new_syns))

    def expand_greek_letters(self, syns):
        """
        if a version of a greek letter is found in a syn, create and add versions of the syn with other formats of the
        character
        :param syns: list of synonyms
        :return: list of synonyms
        """
        # log and initialise
        self.logger.debug("Expanding greek letters within synonyms to include words and greek characters")
        sub_terms = []
        syns_to_remove = []
        for syn in syns:
            # keys are words such as 'alpha', value_lists are characters such as ["α", "∝", "𝛂", "𝛼"]
            for key, value_list in chars.greek_dict.items():
                if re.match(fr"{key}\s?\d+", syn, re.IGNORECASE):  # prevent adding noisy syns like 'gamma 2'
                    syns_to_remove.append(syn)
                else:
                    if (re.match(fr".*\b{key}\b.*", syn, re.IGNORECASE) or
                            re.match(fr".*[^a-zA-Z]{key}[^a-zA-Z].*", syn, re.IGNORECASE)):
                        for character in value_list:
                            sub_terms.append(re.sub(key, character, syn, flags=re.IGNORECASE))
                for char in value_list:
                    if not (re.match(fr"{char}\s?\d+", syn, re.IGNORECASE)):  # prevent adding noisy syns like 'gamma 2'
                        if char.upper() in syn.upper():
                            for other_char in value_list:
                                sub_terms.append(syn.replace(char, other_char))
                            sub_terms.append(syn.replace(char, key))
        cleaner_syns = [i for i in syns if i not in syns_to_remove]
        cleaner_syns.extend(sub_terms)

        return list(set(cleaner_syns))

    def get_anatomy_qualifiers(self, syns_list):
        """
        searches the dictionary of anatomy qualifiers and adds relevant qualifiers to the synonym list

        :param list syns_list: list of synonyms
        :return: list of synonyms with added anatomy qualifiers
        :rtype: list
        """
        # initialise
        new_syns = []

        self.logger.info("Adding anatomy qualifiers to synonyms")

        # remove the worst noise and de-duplicate the current syns list, and add these to the new syns list
        syns = self.remove_noise_and_punctuation(syns_list, "tissue")
        syns.insert(0, self.original_term)
        syns = list(set(syns))
        new_syns.extend(syns)

        # if any of the synonyms in the list are present as keys in the anatomy qualifiers dictionary, then add the
        # relevant anatomy qualifiers values to the synonym list
        for syn in syns:
            if syn.lower() in anatomy_qualifiers.keys():
                new_syns.extend(anatomy_qualifiers[syn.lower()])

        self.logger.info("Done adding anatomy qualifiers to synonyms")

        return new_syns

    def clean_syn_list(self, syns_list=None, syn_type="other"):
        """
        remove noisy and redundant synonyms from a list of synonyms, remove unnecessary punctuation and terms from
        synonyms and expand character variants such as hyphens, presence/lack of spaces and greek characters in synonyms

        :param list syns_list: list of synonyms
        :param str syn_type: type of synonyms (tissue requires an extra cleaning step compared to others)
        :return: non-redundant list of synonyms
        :rtype: list or None
        """

        # if a list of synonyms is supplied, use it. Otherwise, use self.syns
        if syns_list:
            syns = syns_list
        else:
            syns = self.syns

        # remove non-required punctuation and syns that contain the original term
        # (log here so prevent logging multiple times; this function is also called before adding anatomy qualifiers)
        self.logger.debug("Removing noise and unnecessary punctuation characters from collected synonyms")
        filtered_terms = self.remove_noise_and_punctuation(syns, syn_type)

        # remove redundant synonyms from the list of filtered synonyms
        syns = self.remove_redundant_syns(filtered_terms)

        # expand synonyms to include multiple orderings of phrases including "type"s
        type_expanded_terms = self.expand_types(syns)

        # remove "chain" from the end of synonyms
        removed_chain = self.remove_chain_from_end(type_expanded_terms)

        # expand synonyms to include variants of spacing, hyphenation and greek characters (logging is in the function)
        exp_syns = self.greek_char_and_spacing_expansion(removed_chain)

        # de-duplicate
        self.logger.debug("Removing duplicate and redundant synonyms")  # log here so that this is only logged once
        syns = self.remove_redundant_syns(exp_syns)

        # return the synonyms or set them as self.syns
        if syns_list:
            return syns
        else:
            self.syns = syns

    def get_syns(self, parsed_json, descendants=False):
        """
        get synonyms for the provided term and, if descendants = True, its descendants

        :param dict parsed_json: information from the EBI OLS page for an ontology node
        :param bool descendants: whether to get the descendants of the node
        :return: non-redundant synonyms
        :rtype: list
        """
        self.get_syns_from_json(parsed_json)
        if descendants:
            self.get_syns_of_descendants(parsed_json)
        return self.syns


class GetUniprotSynonyms:
    """Get gene synonyms from uniprot"""
    def __init__(self, row, logger):
        """
        init
        """
        self.logger = logger
        self.id_type = row[2]
        self.gene_id = row[1]
        self.tax_id = row[3]

    def build_uniprot_query(self):
        """
        create the URL with the query terms
        :return: url (str)
        """
        # build the query string with the id type, gene id and tax id from a single data frame row
        query = f"https://rest.uniprot.org/uniprotkb/search?query={self.id_type}:{self.gene_id}+organism_id:{self.tax_id}+" \
                f"reviewed:true&fields=gene_names&format=tsv"
        return query

    def run_uniprot_query(self, query):
        """
        run the query and either retrieve the gene synonyms page contents or return an appropriate error
        :return: result or error
        """
        # run the query and either retrieve the gene synonyms or return an appropriate error
        try:
            res = requests.get(query)
            res.raise_for_status()
            # log status code of response
            self.logger.info(f"{threading.current_thread().name}: {res} for {query}")
            res.close()
            return res

        except (requests.exceptions.HTTPError, requests.exceptions.ConnectionError,
                requests.exceptions.SSLError) as err:
            self.logger.error(f"{threading.current_thread().name}: Unable to collect gene synonyms for "
                              f"'{self.gene_id}' (type: {self.id_type}, taxonomy id: {self.tax_id}) from uniprot.org")
            self.logger.error(err)
            # exit() # note: exit commands do not work when running within a thread - which this command often is
            return "uniprot request failed"

    def extract_uniprot_syns(self, res):
        """
        process the query result to obtain a list of synonyms
        :param res: result object from querying Uniprot
        :return: list of syns
        """
        gene_syns_found = False

        # process the retrieved text in to a list of gene synonyms
        gene_syns = res.text.replace("Gene names", "")
        gene_syns = gene_syns.strip()
        if gene_syns == "":
            gene_syns = [self.gene_id]
        else:
            gene_syns = re.split(r"[\n ]", gene_syns)
            self.logger.info(f"{threading.current_thread().name}: Uniprot synonyms {gene_syns} "
                             f"found for {self.gene_id}")
            gene_syns_found = True
        return gene_syns, gene_syns_found

    def get_gene_syns_from_uniprot(self):
        """
        create and run a query to capture uniprot synonyms of the gene, and return them
        :return: list of gene synonyms
        """
        query = self.build_uniprot_query()
        self.logger.info(f"{threading.current_thread().name}: Using the following URL to query uniprot.org: {query}")
        res = self.run_uniprot_query(query)

        if res == "uniprot request failed":
            return [self.gene_id], False

        else:
            syns, syns_found = self.extract_uniprot_syns(res)
            # add the original term back in to the list of synonyms, and de-duplicate
            syns.insert(0, self.gene_id)
        return list(set(syns)), syns_found


class ArgumentCleaner:
    """
    return the argument as a string
    """
    def __init__(self, term, logger):
        """
        init
        """
        self.term = term
        self.logger = logger

    def clean_original_arg(self):
        """
        sets self.term as a string of the supplied argument (converts single-item lists or ints to strings)

        :return: None
        """
        new_term = ""
        if type(self.term) is list:
            # if the list has more than one element, log a warning and skip
            if len(self.term) > 1:
                self.logger.warning(f"The supplied term '{self.term}' should not be a list. Skipping term")
            else:
                # if single-item list, return the element as a string
                new_term = str(self.term[0])
        elif type(self.term) is str:
            new_term = self.term
        elif type(self.term) is int:
            new_term = str(self.term)
        self.term = new_term

    def get_clean_arg(self):
        """
        clean and return the argument as a string

        :return: the original argument as a string
        :rtype: str
        """
        self.clean_original_arg()
        return self.term
