import argparse
import sys
import os
import pandas as pd
import requests
import datetime
from platform import uname
from bs4 import BeautifulSoup


class Arguments:
    """class for setting and getting user input arguments"""
    def __init__(self):
        """
        initialise dictionary of messages to be logged and class attributes
        """
        self.log_msgs = {'error': [], 'warning': [], 'info': []}
        self.gene_ids = None
        self.id_type = None
        self.tax_id = None
        self.kwds = None
        self.infile = None
        self.outfile = None
        self.logging = None
        self.disease = None
        self.tissue = None
        self.other_args = None
        self.log_file = None
        self.expand = None
        self.n_threads = None
        self.charts = False
        self.top_ten = False
        self.quiet_results = False

    @staticmethod
    def determine_required(arg_flag):
        """
        if a list of genes has been supplied, ensure that at least one of disease/tissue/keyword has been supplied

        :param str arg_flag: the single letter flag for the argument (-d, -t or -k)
        :return: bool
        """
        # initialise dict of arg flags and set required to false
        kwd_arg_names = {'-d': "--disease", '-t': "--tissue", '-k': "--keyword"}
        req = False

        # if no input file was provided, then at least one of the flags is required to be present
        if '-i' not in sys.argv and '--infile' not in sys.argv:
            req = True
            # a gene list is required if there is no infile supplied
            if arg_flag == '-g':
                return req
            else:
                # for each flag in disease, tissue, keyword
                for k, v in kwd_arg_names.items():
                    # do not check self
                    if not k == arg_flag:
                        # if any of the other flags was supplied, then arg_flag is not required (return false)
                        if k in sys.argv or v in sys.argv:
                            req = False
                            return req
        return req

    @staticmethod
    def determine_default(value):
        """
        provide a default value only if there is no input file provided

        :param string value: the default value, if required
        :return: None or value
        """
        if '-i' in sys.argv or '--infile' in sys.argv:
            return None
        else:
            return value

    def make_base_parser(self):
        """
        makes the base parser, containing arguments except the EPMC search parameters

        :return: parser for genes, disease, tissue/organ and logging level arguments
        :rtype: argparse.ArgumentParser
        """
        # initialise parser and add usage
        parser = argparse.ArgumentParser(
            usage="\n[-l] logging level (optional. accepted values are 'critical', 'error', 'warning' (default), "
                  "'info', 'debug')"
                  "\n[-i] full path to input file (optional)"
                  "\n[-o] filepath for an excel results file to be created or overwritten (optional)"
                  "\n[-g] space-separated list of genes (required if no input file given)"
                  "\n[-u] type of gene/protein ID (optional. accepted values are 'gene_exact' (default) or 'accession')"
                  "\n[-s] the Uniprot taxonomic ID (see https://www.uniprot.org/taxonomy/) for the organism of "
                  "interest (optional. default '9606' for human)"
                  "\n[-k] space-separated keywords/terms. Use quotation marks for terms that contain spaces. "
                  "(required if no input file, disease or tissue arguments given)"
                  "\n[-d] a single disease name for synonym expansion. Use quotation marks if the disease name contains"
                  " spaces. (required if no input file, keywords or tissue arguments given)"
                  "\n[-t] a single tissue, organ or cell-type name for synonym expansion. Use quotation marks if the "
                  "name contains spaces. (required if no input file, keywords or disease arguments given)"
                  "\n[-f] turn on logging to file (optional)"
                  "\n[-q] turn off automatic results display (optional)"
                  "\n[-m] number of threads to use (optional. default is to use the maximum available)"
                  "\n[-e] turn on synonym expansion for keywords"
                  "\n[-c] control creation of charts about results (pub years, common keywords, wordclouds"
                  "\n[-w] control creation of top ten lists of words in abstracts (e.g. to substitute for wordclouds)"
                  "\nand any valid Europe PMC search parameters, preceded by -- (e.g. --PUB_YEAR 2000). For a list of "
                  "available parameters and their search syntax, visit https://europepmc.org/Help#SSR"
        )
        parser.add_argument(
            '-i', '--infile', dest='infile',
            help="Full filepath to an excel file with the columns gene ID (uniprot geneID or accession number), IDType "
                 "(gene or accession), TaxID (uniprot number code for a species), Keyword (comma-separated keywords to "
                 "search alongside the gene). Other columns tolerated, but ignored"
        )
        parser.add_argument(
            '-o', '--outfile', dest='outfile',
            help="Full filepath for an excel results file to be created or overwritten"
        )
        parser.add_argument(
            '-g', '--genes', nargs='*', dest='genes', required=self.determine_required('-g'),
            help="List of either HGNC symbols or Uniprot IDs"
        )
        parser.add_argument(
            "-u", "--uniprot_id_type", dest='type', choices=["gene_exact", "accession"],
            default=self.determine_default("gene_exact"), type=str.lower,
            help="Type of supplied gene IDs ('accession' for Uniprot IDs, or 'gene_exact' for HGNC symbols)"
        )
        parser.add_argument(
            "-s", "--taxid", dest='taxid', type=int, default=self.determine_default("9606"),
            help="The Uniprot taxonomic ID for the organism of interest. The default value is 9606 (human). For other "
                 "taxonomic IDs and more information, see https://www.uniprot.org/taxonomy/"
        )
        parser.add_argument(
            '-d', '--disease', nargs=1, dest='disease', required=self.determine_required("-d"),
            help="One disease of interest (synonym expansion wil be performed). It is recommended to use the full name "
                 "rather than abbreviations. Please use quote marks if the term contains spaces"
        )
        parser.add_argument(
            '-t', '--tissue', nargs=1, dest='tissue', required=self.determine_required('-t'),
            help="One tissue of interest (synonym expansion will be performed). Please use quote marks around the term "
                 "if it contains spaces"
        )
        parser.add_argument(
            "-k", "--keyword", nargs='+', dest="kwd", required=self.determine_required('-k'),
            help="Additional terms of interest for the search, e.g. a process or molecule. Terms are split on spaces, "
                 "so use quotation marks for multi-word phrases e.g. \"alzheimers disease\" brain"
        )
        parser.add_argument(
            '-l', '--log-level', dest='log',
            help="console logging level (critical, error, warning, info, debug). 'warning' level used by default if no "
                 "valid logging level supplied"
        )
        parser.add_argument(
            '-f', '--log-file', dest='log_file', default=False, action='store_true',
            help="Specify the -f flag to create a log file (with logging at the info level) in addition to the default "
                 "console logging"
        )
        parser.add_argument(
            '-q', '--quiet-results', dest='quiet_results', default=False, action='store_true',
            help="Specify the -q flag to turn off automatic display of results in the default browser "
                 "(note: automatic browser display is not available in the Windows Subsystem for Linux)"
        )
        parser.add_argument(
            '-e', '--expand-keywords', dest='expand', default=False, action='store_true',
            help="Specify the -e flag to attempt synonym expansion of supplied keywords (please note, this may "
                 "introduce noise)"
        )
        parser.add_argument(
            '-c', '--make-charts', dest='charts', choices=["y", "yes", "no", "n", "auto", "a"], default="a",
            help="Whether to create charts for the output files (yes/y, no/n, auto/a). May increase querying time. "
                 "Not recommended for long gene lists due to potential for creation of high number of chart images. "
                 "Automatic level (charts only created for 30 or fewer genes) used by default"
        )
        parser.add_argument(
            '-w', '--top-ten', dest='top_ten', choices=["y", "yes", "no", "n", "auto", "a"], default="a",
            help="Whether to determine the ten most common words in the abstracts of returned publications for display "
                 "in results pages. May increase querying time."
        )
        parser.add_argument(
            '-m', '--multithread', dest='n_threads', type=int, default=0,
            help="Number of cores to use (if left blank, the number of available cores will be determined automatically"
                 " and all will be used)"
        )
        return parser

    @staticmethod
    def get_soup_from_existing_fields_xml():
        """
        collect a soup object from the epmc_fields xml file (and create file if it doesn't exist)

        :return: beautifulsoup object of the parsed file
        :rtype: BeautifulSoup
        """
        # create and initialise the epmc_fields.xml file if it doesn't exist
        if not os.path.exists('epmc_fields.xml'):
            with open('epmc_fields.xml', 'w') as new_fields_file:
                new_fields_file.write("<date>None</date>")

        # open the file and parse the text in to a beautifulsoup object
        with open('epmc_fields.xml', 'r') as epmc_fields:
            text = epmc_fields.read()
        fields_soup = BeautifulSoup(text, 'html.parser')

        return fields_soup

    def get_soup_from_epmc_fields_page(self, file_date, url):
        """
        collect Europe PMC fields xml and save to a file
        :param str file_date: the date in the previous version of the file
        :param str url: url for the Europe PMC search fields page
        :return: BeautifulSoup object for the EPMC fields, or None
        """
        # initialise to None in case of connection issues
        fields_soup = None

        self.log_msgs['info'].append(f"The local copy of Europe PMC fields xml file was created on {file_date}, so will"
                                     f" be updated")
        self.log_msgs['info'].append(f"Fetching Europe PMC fields from {url}")

        # get the EPMC search fields from the url
        try:
            fields_xml = requests.get(url)
            fields_xml.raise_for_status()
            fields_soup = BeautifulSoup(fields_xml.text, 'html.parser')
            # close the connection
            fields_xml.close()
            self.log_msgs['info'].append(f"Done: successfully collected Europe PMC fields")

        # if there's a connection issue, use an older version of the file (unless the file was just created)
        except (requests.exceptions.HTTPError, requests.exceptions.ConnectionError) as err:
            self.log_msgs['warning'].append(f"Unable to retrieve Europe PMC fields due to {err}")

        return fields_soup

    def write_new_fields_to_epmc_fields_file(self, soup, today):
        """
        rewrite the file of EPMC fields to contain today's date and the fields collected today

        :param BeautifulSoup soup: soup of the fields from the EPMC website
        :param datetime.date today: today's date
        :return: None
        """
        self.log_msgs['info'].append(f"Updating epmc_fields.xml, the local copy of the Europe PMC fields")

        # update the file with today's date and the new soup of fields
        with open('epmc_fields.xml', 'w') as epmc_fields:
            epmc_fields.write(f"<date>{today}</date>\n{soup.prettify()}")

        self.log_msgs['info'].append("Done: successfully updated epmc_fields.xml")

    def add_fields_to_parser(self, soup, parser, url):
        """
        add EPMC fields to the argument parser
        :param BeautifulSoup soup: soup of EPMC fields
        :param argparse.ArgumentParser parser: base arg parser
        :param str url: url for the Europe PMC search fields page
        :return: edited parser
        :rtype: argparse.ArgumentParser
        """
        # collect 'term' tags (which contain the parameter names) from the up-to-date soup
        fields = soup.find_all('term')

        # add each parameter name (field) to the parser, with an appropriate help message
        for field in fields:
            field_name = field.contents[0].strip()
            parser.add_argument(f'--{field_name}',
                                help=f'{field_name}, as given in the Europe PMC fields list at {url}')

        self.log_msgs['info'].append(f"Done: successfully added Europe PMC fields to parser")
        return parser

    def add_europe_pmc_field_arguments(self, parser):
        """
        If the local file of EPMC fields has not been updated today, then retrieve the list of Europe PMC search
        parameters and use them to update the file (prevents accessing the website more than once per day). Add the
        up-to-date parameters to the base parser with relevant help messages

        :param argparse.ArgumentParser parser: base parser
        :return: edited parser containing EMPC search parameters
        :rtype: argparse.ArgumentParser
        """
        # initialise and log
        url = "https://www.ebi.ac.uk/europepmc/webservices/rest/fields"
        today = datetime.date.today()
        self.log_msgs['info'].append("Adding Europe PMC search fields to the argument parser")
        fields_soup = self.get_soup_from_existing_fields_xml()

        # check the date on the file, and if the local copy of the file is not up-to-date then update it
        file_date = fields_soup.find('date').text
        if file_date != str(today):
            # get a soup of the fields on the website, and use them to re-write the fields file
            new_fields_soup = self.get_soup_from_epmc_fields_page(file_date, url)

            if new_fields_soup:
                # rewrite the file to contain the new values, and re-assign the fields soup to be the new fields soup
                self.write_new_fields_to_epmc_fields_file(new_fields_soup, today)
                fields_soup = new_fields_soup
            else:
                # use the fields from the previous version of the file, unless it was just created and is empty
                if file_date == "None":
                    self.log_msgs['warning'].append("No Europe PMC search fields were added to the argument parser; if"
                                                    "you need to use Europe PMC search fields in your query, please "
                                                    "exit and try again")
                    fields_soup = None
                else:
                    self.log_msgs['warning'].append(f"Using fields obtained on {file_date} instead")

        # if there is a soup, then use it to add the EPMC fields in it to the parser, else return the parser unedited
        if fields_soup:
            parser = self.add_fields_to_parser(fields_soup, parser, url)
        return parser

    def make_parser(self):
        """
        create the base parser, and add Europe PMC fields to it

        :return: complete parser
        :rtype: argparse.ArgumentParser
        """
        base_parser = self.make_base_parser()
        parser = self.add_europe_pmc_field_arguments(base_parser)
        return parser

    @staticmethod
    def collect_other_args(arg_dict):
        """
        if any other other arguments (e.g. EPMC fields) have been specified, collect them in to a dict

        :param dict arg_dict: full dict of input args
        :return: dict of other args (e.g. EPMC settings)
        """
        other_args = {}

        # for each supplied arg, if it is not in the list of base args, add it to the other args dict
        for k, v in arg_dict.items():
            if v is not None:
                if k not in ["infile", "outfile", "genes", "type", "taxid", "charts", "top_ten", "n_threads",
                             "disease", "tissue", "kwd", "log", 'log_file', 'quiet_results', 'expand']:
                    other_args[k] = v
        return other_args

    def initialise_args(self, args):
        """
        assign parsed args to relevant class attributes

        :param argparse.Namespace args: parsed args
        :return: assign parsed args to class attributes
        :rtype: None
        """
        # create a dictionary of parameters and args
        arg_dict = vars(args)
        self.log_msgs['info'].append(f"The supplied arguments were {arg_dict}")

        # assign input args to relevant class attributes
        self.logging = arg_dict['log']
        self.outfile = arg_dict['outfile']
        self.infile = arg_dict['infile']

        if not self.infile:
            self.tax_id = arg_dict['taxid']
            self.id_type = arg_dict['type']
            self.kwds = arg_dict['kwd']
            if arg_dict['genes']:
                self.gene_ids = [gene.strip(',') for gene in arg_dict['genes']]
            else:
                self.gene_ids = None

        self.disease = arg_dict['disease']
        self.tissue = arg_dict['tissue']
        self.other_args = self.collect_other_args(arg_dict)
        self.log_file = arg_dict['log_file']
        self.expand = arg_dict['expand']
        self.n_threads = arg_dict['n_threads']

        if 'Linux' in uname().system and 'Microsoft' in uname().release:
            self.quiet_results = True
        else:
            self.quiet_results = arg_dict['quiet_results']

        if arg_dict['top_ten'] in ["y", "yes"]:
            self.top_ten = True

        if arg_dict['charts'] in ["a", "auto"]:
            if self.gene_ids:
                if len(self.gene_ids) <= 30:
                    self.charts = True
            else:
                self.charts = "auto"
        elif arg_dict['charts'] in ["y", "yes"]:
            self.charts = True
            if self.gene_ids and len(self.gene_ids) >= 30:
                self.log_msgs['warning'].append(
                    f"You have specified to create charts for {len(self.gene_ids)} genes, "
                    f"which may result in creation of up to {(3 * len(self.gene_ids))} "
                    f"images. If this is undesirable, please cancel and re-run without the "
                    f"-c option")


    def validate_or_add_extension(self):
        """
        ensure that the optional input file path and required output file path both have the correct file extension
        (xlsx), and add the extension to the path if it is not present

        :return: output and input filepath attributes with .xslx extensions
        :rtype: None
        """
        # seemingly redundant code is not redundant, because class attributes are overwritten

        # if an input file was provided, check that it ends with .xlsx and add the extension if not
        if self.infile:
            if self.infile.endswith('.xlsx'):
                self.log_msgs['info'].append(f"File extension successfully checked for input file '{self.infile}'")
            else:
                self.log_msgs['info'].append(f"The supplied input file path '{self.infile}' did not have the required "
                                             f"extension '.xlsx'")
                self.log_msgs['warning'].append(f"The supplied input file path '{self.infile}' did not have the "
                                                f"required extension '.xlsx'. The appropriate extension has been added "
                                                f"automatically")
                self.infile = self.infile + ".xlsx"

        # if an output file was provided, check that it ends with .xlsx and add the extension if not
        if self.outfile:
            if self.outfile.endswith('.xlsx'):
                self.log_msgs['info'].append(f"File extension successfully checked for output file '{self.outfile}'")
            else:
                self.log_msgs['info'].append(f"The supplied output file path '{self.outfile}' did not have the "
                                             f"required extension '.xlsx'")
                self.log_msgs['warning'].append(f"The supplied output file path '{self.outfile}' did not have the "
                                                f"required extension '.xlsx'. The appropriate extension has been added "
                                                f"automatically")
                self.outfile = self.outfile + ".xlsx"

    def validate_outfile_path(self):
        """
        validate that the output file exists or can be created, and warn of overwriting if exists

        :return: relevant logging messages
        :rtype: None
        """
        # check if the file exists and warn for overwriting
        if os.path.exists(self.outfile):
            self.log_msgs['info'].append(f"Done: successfully validated output file path")
            self.log_msgs['warning'].append(f"The output file '{self.outfile}' exists, and will be overwritten")

        # if file doesn't exist, try to create a test file using the provided path
        else:
            self.log_msgs['info'].append(f"The supplied output file path '{self.outfile}' does not exist, so will be "
                                         f"created")
            df = pd.DataFrame(['test'])
            try:
                df.to_excel(f"{self.outfile}")
                self.log_msgs['info'].append("Done: successfully created and validated output file")

            # log an error if unable to create the file
            except PermissionError:
                self.log_msgs['error'].append(f"Permission Denied: Output file path '{self.outfile}' is not accessible."
                                              f" Please ensure that you have appropriate permissions to write to this"
                                              f" location and try again")

    def validate_input_path(self):
        """
        validate the path of the input file

        :return: relevant log messages
        :rtype: None
        """
        if os.path.exists(self.infile):
            self.log_msgs['info'].append(f"Done: successfully validated input file path")
        else:
            self.log_msgs['error'].append(f"The input file path '{self.infile}' could not be found; please check the"
                                          f"path and try again")

    def validate_files(self):
        """
        validate the input and output file paths, and log relevant messages

        :return: None
        """
        # check for 'xlsx' file extensions, and add them if necessary
        self.validate_or_add_extension()

        # if an outfile was specified, validate that it exists for overwriting or that it is possible to create it
        if self.outfile:
            self.log_msgs['info'].append(f"Validating the supplied output filepath '{self.outfile}'")
            self.validate_outfile_path()

        # in an infile was specified, verify that it exists or log an appropriate error
        if self.infile:
            self.log_msgs['info'].append(f"Validating the supplied input filepath '{self.infile}'")
            self.validate_input_path()

    def log_gene_input_redundancy(self):
        """
        if an input file was supplied and values it contains are supplied via command line, then warn that only the
        input file will be used

        :return: relevant logging messages
        :rtype: None
        """
        for k, v in {'gene list': self.gene_ids, 'keyword list': self.kwds,
                     'taxonomy/species id': self.tax_id, 'uniprot ID type': self.id_type}.items():
            if self.infile and v:
                self.log_msgs['warning'].append(f"An input file and {k} were both supplied. Using only the input file.")
                self.log_msgs['info'].append(f"The {k} supplied at the command line will not be used, in favour of the "
                                             f"supplied input file. The {k} was '{v}'")

    def get_args(self):
        """
        creates the parser, parses arguments, performs some basic argument validation, removes empty arguments

        :return: sets object attributes
        :rtype: None
        """
        parser = self.make_parser()
        args = parser.parse_args()
        self.initialise_args(args)
        self.validate_files()
        self.log_gene_input_redundancy()
