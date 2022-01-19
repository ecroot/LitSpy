# LitSpy
## DESCRIPTION
This package allows the user to search through all titles and abstracts available in Europe PMC for articles containing 
co-occurence of the following:
* supplied genes and their synonyms
* an optionally supplied disease and its synonyms
* an optionally supplied tissue or organ and its synonyms
* optionally other supplied keywords
* and any other Europe PMC advanced search parameters desired

and produces summaries of the search results in HTML and optionally Excel format.

## INSTALLATION
The tool has been packaged in to a wheel to make installation easy; simply ensure that you have Python3.6 or higher
installed on your machine, and run the command

`pip install path/to/litspy-version-py3-none-any.whl`

(or, if your default Python version is Python2, replace `pip` with `pip3`)

## USAGE
### Getting started
**Note:** replace `python` with `python3` in all commands if your default python version is not python3

#### Obtaining usage information
To read about the different arguments you can supply to LitSpy (such as input files, output files, logging 
settings and more), run

`python -m litspy -h`
 
#### Basic commands (examples)
Unless specified otherwise, LitSpy assumes a UniProt gene ID for a human gene has been entered. 
Therefore, a minimal search can simply include a human gene ID and a disease, tissue and/or keyword. For example: 

`python -m litspy -g Cftr -d arthritis`

The command can include multiple genes and multiple keywords (separated by spaces), but only one disease and one 
tissue/organ. To search for phrases that include spaces, you should use quotes around the phrases:

`python -m litspy -g APOE APOC1 APOC2 -d "cystic fibrosis" -t liver -k "normal diet" biomarker`

To supply a Uniprot accession number instead of a gene ID, specify the UniProt ID type as accession using `-u` or 
`--uniprot_id_type`, for example:

`python -m litspy -g P13569 -u accession -d arthritis`

To search for non-human genes, specify the appropriate UniProt taxonomy ID using `s` or `--taxid`, for example:

`python -m litspy -g CFTR -d arthritis -s 10090`

#### Using input files
To supply an input file, use `-i` (or `--infile`) followed by the path to the input file, for example:

`python -m litspy -i /home/myusername/documents/litspy_inputs/input1.xlsx`

LitSpy will preferentially take an input file over values supplied at the command line. If the following
command was entered, only the input file would be used and an appropriate warning would be logged

`python -m litspy -i path/to/input.xlsx -g CFTR APOE`
 
Input files should be spreadsheets in xlsx format, containing the following columns in the following order:

UniProtID | IDType    | TaxID | Keyword
--------- | --------- | ----- | -------
CFTR      | gene      | 9606  | normal diet, biomarker
CFTR      | gene      | 9606  | children
APOE      | gene      | 9606  | normal diet, biomarker
P13569    | accession | 9606  | biomarker

Additional columns to the right are tolerated, but ignored. (They may be used for notes, for example)

LitSpy runs one query per row in the spreadsheet (ignoring the headings).

To query for a gene co-occurring with multiple keywords, enter the keywords as a comma-separated list in the 4th column 
of the input file (note: quote marks are not needed in input files). To query for a gene that co-occurs with **any** of 
a list of keywords, you should enter separate rows for each keyword. So, for example, the input file above would return 
documents containing CFTR AND normal diet AND biomarker, and also documents containing CFTR AND children.

#### Create an Excel output file
To get an output file of the results summary and details in Excel format, use `-o` (or `--outfile`) followed by
the full path to an output file:

`python -m litspy -CFTR -d arthritis -o /home/myusername/documents/litspy_outputs/output1.xlsx`
 
If you provide a path to an existing file, the file will be overwritten. If you provide a path to a new file, the new 
file will be created.

The Excel output file is created **in addition** to the HTML results files, which are always created.

#### Create output charts and top 10 list
Output charts are generated automatically if fewer than 30 genes are supplied. If more than 30 genes are supplied, 
charting is turned off to prevent potential creation of 90+ images. To turn charting on or off despite input size, use 
the `-c` or `--make-charts` flag, with 'y' or 'n' to turn charting on or off respectively.

A list of the ten most common terms (excluding search terms) can be generated for each result set by specifying 
the `-w` or `--top-ten` flag with 'y'.

For example:

`python -m litspy -i 'path/to/input.xlsx -w y -c y`

will create result pages with charts and top ten lists no matter the number of entries in the supplied file

#### Logging
By default, LitSpy logs to the console only at the warning level. To change the console logging level, supply a log 
level argument by entering `-l` (or `--log-level`) followed by a valid logging level, for example:

`python -m litspy -g CFTR -d arthritis -l info`

To turn on logging to a log file, specify the `-f` (or `--log-file`) flag:

`python -m litspy -g CFTR -d arthritis -l warning -f`

The log file is always logged at the info level.

#### Multithreading
Some parts of the tool are able to run in parallel on many threads. You can specify the number of threads to use with 
the `-m` or `--multithread` flag, for example:

`python -m litspy -g CFTR -d arthritis -m 4` 

If you do not specify the number of threads to use, then the tool will automatically determine the number of available 
cores and use this as the maximum number of threads to use. It is possible that this will affect performance of 
other tasks that the machine is running at the same time.

#### Keyword synonym expansion
To attempt collection of synonyms for supplied keywords, specify the `-e` (or `--expand-keywords`) flag in the command:

`python -m litspy -g CFTR -d arthritis -k biomarker -e`

This will return the synonyms from the first node matched (the best match) when searching for the term in the EBI OLS. 

**Warning:** It is expected that this generic searching will be far noisier on average than the more specific searches for
disease or tissue synonyms, and therefore **this option should be used cautiously**.

## BACKGROUND
The aim of this tool is to form part of my PhD project, where it will help to inform target novelty within particular 
contexts, and be combined with other tools. However, it can also be useful as a standalone search tool; for example, to
identify articles of interest more comprehensively than by keyword searches that do not use synonym expansion.

## AUTHOR
Emma Croot, 2020, ec339@le.ac.uk, github.com/ec339

## LICENCE
MIT (see LICENCE.txt)

## CHANGE LOG
This is the first version of the tool, so there are currently no changes

## TESTS
[test suite is under development]

## ACKNOWLEDGEMENTS
**OmixLitMiner**

Although the design and creation of this tool began before [Pascal et al's 2020 publication](https://europepmc.org/article/MED/32092871)
(DOI: 10.3390/ijms21041374, PMID: 32092871), their tool OmixLitMiner provided inspiration for parts of this project;
particularly the use of UniProt to obtain alternative official gene names, and the output charts and styles.

**Supervision team**

My PhD supervisors, Dr. Thanos Didangelos, Dr. Richard Badge, and Prof. Louise Wain

**Others**

Thomas Rowlands, co-developer of rtgo (used in LitSpy for multithreading)