import pandas
import datetime
import wordcloud
import textwrap
import os
import webbrowser
import re
import multiprocessing

from matplotlib import pyplot as plt, ticker
from collections import Counter
from rtgo import ReadyThready

from litspy.input_args import Arguments
from litspy.logger import Logger
from litspy.epmc_query import Query
from litspy.create_html import HtmlResults
from litspy.noisy_phrases import stop_words, not_top_ten

import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)


def create_keyword_frequency_chart(name_result_tuple):
    """
    create a bar chart showing the frequency of publication year for the results in the supplied df

    :param tuple name_result_tuple: tuple of gene key (str), keywords (list)
    :return: None (chart png file created)
    """
    chart_name, kwd_list = name_result_tuple
    kwd_fig = plt.figure()

    if kwd_list:
        xlabel_list = []

        # count keywords and find the 20 most common
        raw_kwd_count = Counter(kwd_list)
        top_20 = raw_kwd_count.most_common(20)
        kwd_count = Counter(dict(top_20))

        # get the axes to format the yticks labels
        ax = kwd_fig.add_subplot(111)

        # plot chart of keywords & their frequency
        plt.bar(range(len(kwd_count)), kwd_count.values())

        # wrap the top 20 keywords to max 21 characters per line
        for label in kwd_count.keys():
            wrapped_label = textwrap.wrap(label, width=21)
            x_label = "\n".join(wrapped_label)
            xlabel_list.append(x_label)

        # add the keywords and frequencies to the axes, and optimise formats
        plt.xticks(range(len(kwd_count)), xlabel_list, rotation=90, ma='right', ha='center', va='top', size=8)
        ax.yaxis.set_major_locator(ticker.MaxNLocator(integer=True))

    # add title and axes labels
    plt.title("Top 20 keywords present in results")
    plt.ylabel("Frequency")
    plt.xlabel("Keyword")

    # reduce margins around bars and ensure that the bottom doesn't get cut off
    plt.margins(x=0.01)
    plt.tight_layout()

    # save and close
    kwd_fig.savefig(f"html_results/{chart_name}_top20_kwds.png")
    plt.close()

    return


def create_wordcloud(name_result_tuple):
    """
    create a wordcloud from the abstracts of results

    :param tuple name_result_tuple: tuple of chart name with gene key (str), abstracts (list)
    :return: None (png file created)
    """
    _, abstract_text, chart_name = name_result_tuple
    plt.figure()

    # turn the abstracts in to a long string, or assign the string "none" to the variable if no abstracts
    if abstract_text:
        # create a wordcloud object from the long string of abstracts
        wc = wordcloud.WordCloud(max_font_size=60,
                                 min_font_size=8,
                                 height=400,
                                 width=600,
                                 prefer_horizontal=0.6,
                                 relative_scaling=0.35,
                                 background_color="white").generate_from_frequencies(abstract_text)
        # create the wordcloud plot
        plt.imshow(wc, interpolation="bilinear")

    # remove axes
    plt.axis("off")

    # ensure no edges are cut off, save and close
    plt.tight_layout()
    plt.savefig(f"html_results/{chart_name}_wordcloud.png")
    plt.close()

    return


class LitSpy:
    def __init__(self):
        """
        initialise with command-line arguments, date stamp of start time, and logger
        """
        # get the date and time stamp to use in file names
        self.date_stamp = datetime.datetime.now().strftime("%Y%m%d-%H%M%S")

        # get arguments and relevant logging messages
        args = Arguments()
        args.get_args()
        self.args = args

        # create logger, initialise with level, whether to log to file and the date stamp
        lgr = Logger()
        self.logger = lgr.initialise_logger(level=self.args.logging, logfile=self.args.log_file,
                                            datestamp=self.date_stamp)
        # log any messages created in arg parsing
        self.log_arg_messages(self.args.log_msgs)

    def log_arg_messages(self, msg_dict):
        """
        perform logging for the messages collected when parsing input arguments

        :param dict msg_dict: dictionary of logging levels (str) and associated messages (list)
        :return: None
        """
        # if logging level key has messages, log each message at the relevant level
        if msg_dict['info']:
            for msg in msg_dict['info']:
                self.logger.info(msg)
        if msg_dict['warning']:
            for msg in msg_dict['warning']:
                self.logger.warning(msg)
        if msg_dict['error']:
            # if there are error messages, log them all then exit
            for msg in msg_dict['error']:
                self.logger.error(msg)
            exit("Error in supplied arguments. See log for more details")

    # function for determining if doc type is review. potentially use later to help determine score
    # def determine_if_review(self, res_soup):
    #     """
    #     determine whether the result is labelled as a review
    #
    #     :param BeautifulSoup res_soup: result soup object
    #     :return: string (y, n, n/a as appropriate)
    #     """
    #     # initialise
    #     review = "N"
    #
    #     # determine if the publication is listed as a review
    #     p_types = self.extract_results_if_available(res_soup, 'pubtype')
    #     if p_types == "unavailable":
    #         review = "n/a"
    #     else:
    #         # epmc splits pubtypes on both ',' and ';', so normalise to a single split character and split
    #         p_types = p_types.replace(",", ";").split(";")
    #         for pt in p_types:
    #             if pt.lower().strip() == "review":
    #                 review = "Y"
    #     return review

    def create_html_output_directories(self):
        """
        creates the directory structure for the charts and html results files

        :return: None
        :raises OSError: if unable to create the output directories, e.g. due to a permissions issue
        """
        self.logger.info(f"Creating output directories for html results files, under 'html_results/{self.date_stamp}'")
        try:
            os.makedirs(f'html_results/{self.date_stamp}')
            os.makedirs(f'html_results/{self.date_stamp}/charts')
            os.makedirs(f'html_results/{self.date_stamp}/results')
            self.logger.info("Done: successfully created html output directories")
        except OSError as err:
            self.logger.error('Could not create directories for HTML output files')
            raise err

    def add_unique_key_to_df(self, my_df):
        """
        create and add unique keys based on gene names to each row

        :param pandas.DataFrame my_df: data frame containing gene name, id type, tax id and (optionally) keywords
        :return: data frame with added unique keys to be used in output file names
        :rtype: pandas.DataFrame
        """
        self.logger.info("Creating unique keys from input genes")

        # if there are any duplicated genes, then use the cumulative count of each gene name to create unique keys
        if my_df.duplicated('UniProtID').any():
            unique_key_list = []
            gb = my_df.groupby('UniProtID').cumcount()
            for i, row in my_df.iterrows():
                if gb[i] == 0:
                    unique_key_list.append(row[0])
                else:
                    unique_key_list.append(f"{row[0]}_{gb[i]+1}")
            my_df['UniqueKey'] = unique_key_list
        else:
            # if there are no duplicates, then the gene names can be used as unique keys
            my_df['UniqueKey'] = my_df['UniProtID']

        self.logger.info("Done creating unique keys from input genes")

        return my_df

    def get_input_data_as_data_frame(self):
        """
        whether the data has been supplied within an excel file or commandline arguments, produce a data frame
        containing the data

        :return: data frame containing columns of gene, id type, tax id, keywords and a unique key based on gene name
        :rtype: pandas.DataFrame
        """
        # create df of input values: from an input file if supplied, else from relevant args
        if self.args.infile:
            df = pandas.read_excel(self.args.infile, usecols=[0, 1, 2, 3])
            # ensure columns are named as expected
            df.columns = ['UniProtID', 'IDType', 'TaxID', 'Keyword']
            df = df.replace("Gene", "gene_exact")
            # determine whether to make charts
            if self.args.charts:
                num_entries = len(df)
                if num_entries <= 30:
                    if self.args.charts == "auto":
                        self.args.charts = True
                else:
                    if self.args.charts == "auto":
                        self.args.charts = False
                    elif self.args.charts:
                        self.logger.warning(f"You have specified to create charts for {num_entries} searches, which "
                                            f"will result in creation of up to {(3 * num_entries)} images. If this is "
                                            f"undesirable, please cancel and re-run without the -c option")

        else:
            self.args.kwds = ",".join(self.args.kwds) if self.args.kwds else None
            # initialise a data frame with one of the supplied gene names in each row
            df = pandas.DataFrame(self.args.gene_ids, columns=['UniProtID'])
            # add columns containing the other supplied values
            df['IDType'] = self.args.id_type
            df['TaxID'] = self.args.tax_id
            df['Keyword'] = self.args.kwds

        # clean the data
        self.logger.info("Processing submitted genes to remove clone-based gene names and gene map locations")
        df.drop_duplicates(inplace=True)  # remove duplicate rows

        clone_name_df = df[df['UniProtID'].str.match(r"[A-Z]{2,}\d{6}\.\d")]
        map_name_df = df[df['UniProtID'].str.match(r"\d{1,2}[pq]\d+\.?\d*")]
        dirty_df = pandas.concat([clone_name_df, map_name_df])

        df.drop(dirty_df.index, inplace=True)

        self.logger.info(f"Done processing submitted genes. {len(dirty_df)} removed.")

        if not dirty_df.empty:
            removed_genes = ', '.join(dirty_df['UniProtID'].tolist())
            self.logger.warning(f"The supplied gene(s) '{removed_genes}' will not be searched, as they were identified "
                                f"as clone-based genes or gene map locations")

        if df.empty:
            self.logger.warning("No genes remain after cleaning")
            exit("Exit: nothing to query")

        # create and add unique keys for each row (search)
        final_df = self.add_unique_key_to_df(df)

        return final_df

    @staticmethod
    def remove_redundant_preprints(results_dicts):
        """
        identifies and removes preprint documents that have a corresponding publication

        :param list results_dicts: list of dicts of results and information about each search
        :return: results dicts with redundant preprints (preprints of papers that have been published) removed
        :rtype: list[dict]
        """
        for res_dict in results_dicts:
            # collect list of all doc IDs
            all_ids = [x['ID'] for x in res_dict['true_result_doc_info']]

            # only keep docs that are not a preprint of an item that is within the list of IDs
            res_dict['true_result_doc_info'] = [x for x in res_dict['true_result_doc_info'] if x['prep_of']
                                                not in all_ids]

            # remove the preprint_of information to reduce df size later on
            for each in res_dict['true_result_doc_info']:
                each.pop('prep_of')

        return results_dicts

    @staticmethod
    def get_summary_dict(results_dicts):
        """
        create a summary of results information to be used in the summary results pages

        :param list results_dicts: list of dicts of results information for each search
        :return: dict of summary information
        """
        summary_results = {}
        for each in results_dicts:
            if each['res_count'] < 1000:
                # count the unique, non-empty results
                res_list = [d for d in each['true_result_doc_info'] if d.get('ID') != 'none']
                num_results = len(res_list)
            else:
                num_results = "Over 1000"

            summary_results[each['gene key']] = {'gene name': each['gene name'], 'search terms': each['search terms'],
                                                 'results': num_results, 'result page': each['gene key'],
                                                 'query strings': each['query_strings']
                                                 }
        return summary_results

    @staticmethod
    def get_output_data_as_dfs_and_dicts(summary_results, results_dicts):
        """
        process the parsed results dicts in to a df and a dict of dfs

        :param dict summary_results: dict of relevant info from summary of results
        :param list results_dicts: lists of dicts of relevant info from each query result
        :return: summary df, and dict of dfs of result details
        :rtype: tuple[pandas.DataFrame, dict]
        """
        # initialise the dict of details dfs
        details_dfs = {}
        search_terms_dict = {}

        # create df of summary results, and start indexing at 1 rather than 0
        res_summary_df = pandas.DataFrame(summary_results.values())
        res_summary_df.index += 1

        # create dict of dfs of detailed results
        for each in results_dicts:
            # create a df of the relevant info about the result results dict
            result_df = pandas.DataFrame(each['true_result_doc_info'])
            # remove rows where the ID is "none"
            result_df = result_df[result_df.ID != "none"]
            # prevent skipping indexes of removed rows, and start indexing at 1
            result_df.index = range(len(result_df.index))
            result_df.index += 1
            # add the df to the dict of dfs under the relevant gene key
            details_dfs[each['gene key']] = result_df
            search_terms_dict[each['gene key']] = each['query-specific syns']

        return res_summary_df, details_dfs, search_terms_dict

    def create_outfile(self, summary_df, details_dfs):
        """
        create an output excel file containing tables of results similar to OmixLitMiner output

        :param pandas.DataFrame summary_df: df containing summary of results for each query
        :param dict details_dfs: dict of dfs containing details of results for each query
        :return: populated output file
        :rtype: None
        """
        self.logger.info(f"Printing results to output file '{self.args.outfile}'")
        with pandas.ExcelWriter(self.args.outfile) as writer:
            summary_df.to_excel(writer, sheet_name="summary", index=False)
            for k, v in details_dfs.items():
                if not v.empty:
                    v.to_excel(writer, sheet_name=k)
        self.logger.info(f"Done: successfully created output file {self.args.outfile}")

    @staticmethod
    def create_year_frequency_chart(chart_name, year_counts):
        """
        create a bar chart showing the frequency of publication year for the results in the supplied df

        :param chart_name: string containing date stamp and gene key for output file creation
        :param year_counts: counter objects for number of docs/year
        :return: save a png file of the chart
        :rtype: None (png file created)
        """
        if year_counts.empty:
            plt.title("Number of publications per year")
        else:
            # plot the years and their frequencies with relevant title
            ax = year_counts.plot.bar(title="Number of publications per year")

            # ensure the y axis uses integer ticks
            ax.yaxis.set_major_locator(ticker.MaxNLocator(integer=True))

        # add axes labels
        plt.ylabel("Number of publications")
        plt.xlabel("Year")

        # ensure no edges are cut off, save and close
        plt.tight_layout()
        plt.savefig(f"html_results/{chart_name}_by_year.png")
        plt.close()

    @staticmethod
    def split_syn_phrases_and_singles(syns):
        """
        splits the synonyms in to two lists: phrases and single words

        :param list syns: list of synonyms
        :return: pipe-separated string of synonym phrases, list of single-word synonyms
        :rtype: tuple[str, list]
        """
        syn_phrases = []
        single_syns = []

        for syn in syns:
            if " " in syn:
                syn_phrases.append(syn)
            else:
                single_syns.append(syn.upper())

        return syn_phrases, single_syns

    def get_counted_relevant_words(self, g_key, dfs_dict, syn_key_dict, universal_syn_phrases, universal_syn_words):
        """
        return a count of relevant (not queried or stopword) words from the abstracts to be used to determine the top
        ten and/or create the wordcloud

        :param str g_key: key for dicts
        :param dict dfs_dict: dict of dataframes containing results for each gene queried
        :param syn_key_dict: dict of gene key and all the syns specifically used in that query (genes and kwds)
        :param universal_syn_phrases: synonyms phrases used in all queries (e.g. multiword disease syns)
        :param universal_syn_words: synonym words used in all queries (e.g. single word disease syns)
        :return: tuple[gene key, Counter, filepath part]
        """
        cntr = Counter()
        df = dfs_dict[g_key]
        abs_list = df['Abstract'].tolist()
        cleaned_abstracts = []
        file_path_name_part = f"{self.date_stamp}/charts/{g_key.replace('*', '')}"

        phrase_syns, single_syns = self.split_syn_phrases_and_singles(syn_key_dict[g_key])
        all_phrase_syns = universal_syn_phrases + phrase_syns
        all_single_syns = universal_syn_words + single_syns
        or_syns = '|'.join(all_phrase_syns)

        for a in abs_list:
            if not a == "unavailable":
                a.replace("-", " ")  # syns don't contain hyphenation, so remove hyphens for better matching
                a = re.sub(or_syns, "", a, flags=re.IGNORECASE)  # remove all syn phrases from the text
                a = re.sub(r"\s{2,}", " ", a)  # remove multiple spaces following syn removal
                cleaned_abstracts.append(a)

        text = " ".join(cleaned_abstracts)

        # if a string of abstracts exists, then turn it in to a counter (excluding stopwords and other noise)
        if text:
            keys_to_exclude = []
            noise = not_top_ten + stop_words

            # remove the worst of the noise
            # (stopword removal only works for lowercase words, but don't .lower() to maintain important capitalisation)
            counts = wordcloud.WordCloud(stopwords=set(noise)).process_text(text)

            # remove case variants of noise whilst maintaining capitalisation of counted words
            for key in counts.keys():
                if key.upper() in all_single_syns:
                    keys_to_exclude.append(key)
                elif key.lower() in noise:
                    keys_to_exclude.append(key)
                elif any(s in key.upper() for s in all_single_syns) and any(n in key.lower() for n in noise):
                    keys_to_exclude.append(key)

            for each in keys_to_exclude:
                counts.pop(each)

            cntr.update(counts)

        return g_key, cntr, file_path_name_part

    def create_html_output_page(self, g_key, counts, dfs, results_html):
        """
        create a html page for the results of the search

        :param str g_key: gene-based key
        :param dict counts: dictionary of gene key: count of each relevant word in the abstracts
        :param dict dfs: dictionary of gene key: dataframe containing results
        :param HTMLResults results_html: initialised html results object
        :return: saves a html output page
        """
        top_ten = []
        df = dfs[g_key]

        if counts:
            count = counts[g_key]
            if count:
                top_ten = count.most_common(10)

        results_html.make_html_page_for_results(g_key, df, self.args.charts, top_ten)

    def make_outputs_and_get_location(self, summary_df, details_dfs, summary_dict, kwd_and_gene_dict, tissue_syns,
                                      disease_syns, other_syns):
        """
        create summary and individual HTML results pages, and return the location of the summary page

        :param summary_df: dataframe containing summary information for results
        :param details_dfs: dict of dataframes containing results for each gene queried
        :param summary_dict: dictionary of summary information
        :param kwd_and_gene_dict: dict of gene key and all the syns specifically used in that query (genes and kwds)
        :param tissue_syns: list of tissue synonyms queried
        :param disease_syns: list of disease synonyms queried
        :param other_syns: list of other command line-derived synonyms queried
        :return: path to created results summary output file
        """
        ReadyThready.set_logger(self.logger)
        syns_queried = tissue_syns + disease_syns + other_syns
        warning_required = False
        abs_counts = []
        abs_count_dict = {}
        kwds = []

        # create the excel output file if requested
        if self.args.outfile:
            self.create_outfile(summary_df, details_dfs)

        # initialise the html creator
        self.logger.info("Creating HTML results files")
        results_html = HtmlResults(summary_dict, self.date_stamp)

        # get counts of words in abstracts (except stopwords and queried terms)
        if self.args.charts or self.args.top_ten:
            universal_syn_phrases, universal_syn_words = self.split_syn_phrases_and_singles(syns_queried)
            abs_counts = ReadyThready.go(func=self.get_counted_relevant_words, arg_data_index=0,
                                         n_threads=self.args.n_threads,
                                         args=[list(details_dfs.keys()), details_dfs, kwd_and_gene_dict,
                                               universal_syn_phrases, universal_syn_words])
            if self.args.top_ten:
                abs_count_dict = {g_key: abs_count for g_key, abs_count, _ in abs_counts}

        if self.args.charts:
            for k, v in details_dfs.items():
                # get list of tuples of file name part and lowercase article keywords
                name = f"{self.date_stamp}/charts/{k.replace('*', '')}"
                kwds.append((name, [item.lower() for each in v['Keywords'].tolist() for item in each]))
                # workaround: create year charts here to prevent windows permission issues in multiprocess
                year_counts = v.Year.value_counts()
                self.create_year_frequency_chart(name, year_counts.sort_index(ascending=True))

            # create charts via multiprocess
            pool = multiprocessing.Pool()
            # try to create the charts using multiprocess, or warn and continue if any charts are not created
            try:
                pool.map(create_keyword_frequency_chart, kwds)
                pool.map(create_wordcloud, abs_counts)
            except RecursionError as e:
                # windows access errors can occur, possibly https://bugs.python.org/issue38188
                warning_required = e
            pool.close()
            pool.join()
            if warning_required:
                self.logger.warning(f"Unable to create some output charts. All other output information is "
                                    f"unaffected. Error information: {warning_required}")

        ReadyThready.go(func=self.create_html_output_page,
                        args=[list(details_dfs.keys()), abs_count_dict, details_dfs, results_html],
                        arg_data_index=0, n_threads=self.args.n_threads)

        # create the summary html page
        results_html.make_html_page_for_summary(summary_df)

        self.logger.info("Done: successfully created HTML output files")
        return results_html.get_summary_path()

    def main(self):
        """
        parse input data, collect synonyms for genes, diseases and tissues supplied, use these to build queries for
        Europe PMC, run the Europe PMC queries and parse the results, process and print the results to HTML and
        (optionally) excel output files
        """
        # create output directories
        self.create_html_output_directories()

        # convert input data to a data frame
        df = self.get_input_data_as_data_frame()

        # initialise the query object with the data frame of input data, disease, tissue and other EPMC args
        my_query = Query(df=df, logger=self.logger, disease=self.args.disease, tissue=self.args.tissue,
                         kwds=self.args.kwds, others=self.args.other_args, expand=self.args.expand,
                         n_threads=self.args.n_threads)

        # run the queries and get the results
        results_dicts = my_query.get_results_from_epmc()

        # check potential results and add them to the list of actual results, then process in to a summary dict
        self.logger.info("Parsing results from Europe PMC queries")
        no_redundant_preprints = self.remove_redundant_preprints(results_dicts)
        summary_dict = self.get_summary_dict(no_redundant_preprints)
        self.logger.info("Done: successfully parsed results from Europe PMC queries")

        # create dfs of the processed results
        summary_df, details_dfs, kwd_gene_syn_dict = self.get_output_data_as_dfs_and_dicts(summary_dict, results_dicts)

        # make output files and return location of html output file summary
        html_output_location = self.make_outputs_and_get_location(summary_df, details_dfs, summary_dict,
                                                                  kwd_gene_syn_dict, my_query.tissue_syns,
                                                                  my_query.disease_syns, my_query.kwd_syns)
        if not self.args.quiet_results:
            # open the summary results page in the default browser
            try:
                webbrowser.open_new(os.path.realpath(html_output_location))
                self.logger.info("Opening the summary of the results in the default browser")
            except OSError:
                self.logger.warning(f"Unable to automatically open results in browser")

        self.logger.info(f"Done. Results files can be found at {html_output_location}")


if __name__ == '__main__':
    # initialise and run the lit miner
    my_search = LitSpy()
    my_search.main()
