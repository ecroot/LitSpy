from bs4 import BeautifulSoup


class HtmlResults:
    def __init__(self, summary, date):
        """
        initialise
        """
        self.summary = summary
        self.date = date
        self.summary_path = f'html_results/{self.date}/results_summary.html'

    @staticmethod
    def edit_results_table_html(raw_html):
        """
        taking the raw html for a table from the dataframe, edit the content to include links and add a title and
        wrapper tags
        :param raw_html: the html generated using the .to_html method on a dataframe containing results
        :return: html table with links
        """
        table_soup = BeautifulSoup(raw_html, 'html.parser')
        for row in table_soup.find_all('tr'):
            vals = row.find_all('td')
            if vals:
                try:
                    # remove the content from the first row and assign it to the id variable
                    id = vals[0].string.extract()
                    url = vals[-1].contents[0]
                    # create a link tag with the extracted id and URL from above
                    link_tag = table_soup.new_tag("a", href=url)
                    link_tag.append(id)
                    # add the link tag to the first column
                    vals[0].append(link_tag)
                    # convert the keyword lists in to strings
                    kwds = vals[-2].contents[0]
                    kwd_string = kwds.strip('[]')
                    vals[-2].string.replace_with(kwd_string)

                    # wrap the cells with long contents with div tags to allow for scrolling in the table
                    for i in [1, 3, 5, 6]:
                        vals[i].string.wrap(table_soup.new_tag("div"))

                except AttributeError:
                    row.extract()

            # delete the URL cell/header: decompose causes unresolvable navigable string error, so use extract
            all_vals = row.find_all(['td', 'th'])
            if all_vals:
                all_vals[-1].extract()

        return table_soup

    def add_results_title_to_soup_body(self, soup, gene_key):
        """
        add a title to the html soup

        :param BeautifulSoup soup: beautifulsoup object containing the edited table and html & body tags
        :param str gene_key: gene-based key for the summary df
        :return: soup with added title tag
        :rtype: BeautifulSoup
        """
        # create header tag
        title_tag = soup.new_tag("h1")

        # collect the gene game and search terms
        gene_name = self.summary[gene_key]['gene name']
        search_terms = self.summary[gene_key]["search terms"]

        # add the gene name and search terms as a string to the header tag
        title_tag.string = f"Results for '{gene_name}, {search_terms}'"

        # add the header tag to the top of the soup body
        soup.body.insert(0, title_tag)
        return soup

    @staticmethod
    def add_top_ten(soup, top_ten):
        """
        add the top ten list to the soup
        :param BeautifulSoup soup: beautifulsoup object containing table and title
        :param top_ten: list of top ten terms
        :return: soup with top ten words added
        :rtype: BeautifulSoup
        """
        top_ten_strings = []
        for each in top_ten:
            top_ten_strings.append(f"{each[0]} ({each[1]})")

        top_ten_title_tag = soup.new_tag("h2")
        top_ten_title_tag.string = "Top ten most frequently-used terms"

        top_ten_tag = soup.new_tag("p")
        top_ten_tag.string = f"{', '.join(top_ten_strings)}"

        # add below the title
        soup.body.insert(1, top_ten_title_tag)
        soup.body.insert(2, top_ten_tag)
        return soup

    @staticmethod
    def add_charts_to_soup(soup, gene_key):
        """
        add the chart images to the HTML soup

        :param BeautifulSoup soup: beautifulsoup object containing table and title
        :param str gene_key: str - gene-based key for the summary df
        :return: soup with added image tags for each chart
        :rtype: BeautifulSoup
        """
        # for each type of chart, create and add a relevant tag with the path to the chart image as the source
        for chart_name in ["wordcloud", "top20_kwds", "by_year"]:
            chart_tag = soup.new_tag("img")
            chart_tag["src"] = f"../charts/{gene_key.replace('*', '')}_{chart_name}.png"
            chart_tag["alt"] = f"{chart_name} chart"
            chart_tag["width"] = "32.5%"
            chart_tag["height"] = "auto"
            # add the chart above the table and below the title
            soup.body.insert(1, chart_tag)
        return soup

    @staticmethod
    def add_style_to_soup(soup, style):
        """
        add head-wrapped style tag and relevant styling information to the soup

        :param BeautifulSoup soup: beautifulsoup object containing table, title and chart images
        :param str style: str containing relevant styling information
        :return: soup with added style information
        :rtype: BeautifulSoup
        """
        # create style tag and populate with the style string
        style_tag = soup.new_tag("style")
        style_tag.append(style)

        # create meta tag with a character set specified
        charset_tag = soup.new_tag("meta", charset='UTF-8')

        # add the new tags to the soup
        soup.html.insert(0, style_tag)
        soup.html.insert(0, charset_tag)

        # add a new "head" tag around the style tag
        soup.style.wrap(soup.new_tag("head"))
        return soup

    @staticmethod
    def create_results_style_string():
        """
        create a string of style information for the results pages

        :return: style info (str)
        """
        style = (
            # set font for page
            "body {font-family: Arial, Helvetica, sans-serif;}\n"
            # centre headings
            "h1 {text-align: center;}\n"
            # format table
            "table {width: 100%; table-layout: fixed; border-collapse: collapse; background-color: white;}\n"
            # specify width of index and year columns
            "table th:nth-child(1) {width: 2%;}\ntable th:nth-child(4) {width: 3%;}\n"
            # specify width of ID and review columns
            "table th:nth-child(2), table th:nth-child(6) {width: 10%;}\n"
            # specify widths of title and abstract columns
            "table th:nth-child(3) {width: 20%;}\ntable th:nth-child(7) {width: 30%;}"
            # specify widths of author and keyword columns
            "table th:nth-child(5), table th:nth-child(8), {width: 15%;}\n"
            "th {font-weight: bold; text-align: left; background-color: #f2f2f2; white-space: nowrap;}\n"
            "th, td {padding: 4px;}\n"
            # cell height is 3.6 the height of the line height, allowing for 3 lines + padding per cell
            "td div {height: 3.6em; overflow: auto;}\n"
            # alternate grey background on rows
            "tr:nth-child(even) {background-color: #f2f2f2;}\n"
            # highlight row that mouse is hovering over
            "tr:hover {background-color: #E0ECF8;}\n"
            "img {border-style: solid; border-width: thin;}"
        )
        return style

    @staticmethod
    def wrap_soup_in_body_and_html_tags(soup):
        """
        add html and body tags to the html soup

        :param BeautifulSoup soup: soup of the edited table
        :return: soup with added html and body
        :rtype: BeautifulSoup
        """
        # wrap the highest-level tag in a body tag
        soup.contents[0].wrap(soup.new_tag("body"))

        # wrap body in html tag
        soup.body.wrap(soup.new_tag("html"))
        return soup

    def make_html_page_for_results(self, gene_key, results_df, charts, top_ten):
        """
        convert the results df in to HTML, edit the table, add title, charts and style to the html and save this to
        an output HTML file

        :param str gene_key: key based on gene name
        :param pandas.DataFrame results_df: dataframe of results
        :param charts: whether or not to include charts in the output
        :param top_ten: list of the top 10 words found in abstracts
        :return: save a HTML file of results
        :rtype: None
        """
        # create "no results" soup if there are no results for the query
        if results_df.empty:
            # create soup for empty results page
            basic_soup = BeautifulSoup('<h2>No results</h2>', 'html.parser')

            # wrap in html and body tags
            wrapped_basic_soup = self.wrap_soup_in_body_and_html_tags(basic_soup)

            # add title
            full_soup = self.add_results_title_to_soup_body(wrapped_basic_soup, gene_key)

        else:
            # create a html table of the detailed results df
            raw_table_html = results_df.to_html(justify='center')

            # edit the table and return a soup
            table_soup = self.edit_results_table_html(raw_table_html)

            # wrap the table in html and body tags
            wrapped_table_soup = self.wrap_soup_in_body_and_html_tags(table_soup)

            # add a title to the soup
            title_and_table_soup = self.add_results_title_to_soup_body(wrapped_table_soup, gene_key)

            if top_ten:
                pre_chart_soup = self.add_top_ten(title_and_table_soup, top_ten)
            else:
                pre_chart_soup = title_and_table_soup

            if charts:
                # add the chart images to the soup
                full_soup = self.add_charts_to_soup(pre_chart_soup, gene_key)
            else:
                full_soup = pre_chart_soup

        # add style to the soup
        style = self.create_results_style_string()
        final_soup = self.add_style_to_soup(full_soup, style)

        # write a html file from the soup
        with open(f"html_results/{self.date}/results/{gene_key.replace('*', '')}.html", 'w+', encoding='utf-8') as f:
            f.write(final_soup.prettify(formatter='html'))

    @staticmethod
    def edit_summary_table_html(table_html):
        """
        make the raw html table in to a soup object and edit to include links and remove unwanted information

        :param str table_html: the html created using the df.to_html method on the summary df
        :return: html soup of edited table
        :rtype: BeautifulSoup
        """
        # make a soup of the html table
        table_soup = BeautifulSoup(table_html, 'html.parser')

        # for each row in the table, get all normal (non-header) cells and edit appropriately
        for row in table_soup.find_all('tr'):
            vals = row.find_all('td')
            if vals:
                # create a link to the results page
                page_name = vals[-2].string.extract()
                link_tag = table_soup.new_tag("a", href=f"results/{page_name.replace('*', '')}.html")
                link_tag.append(page_name)
                vals[-2].append(link_tag)
                # wrap the query string and search terms divs to allow for scrolling to be applied in style
                vals[1].string.wrap(table_soup.new_tag("div"))
                vals[-1].string.wrap(table_soup.new_tag("div"))
                # split the queries and wrap each in paragraph tags, to separate in output
                queries = vals[-1].div.string.extract()
                query_strings = queries.split(" split_here ")
                for query_string in query_strings:
                    qs_para = table_soup.new_tag('p')
                    qs_para.append(query_string)
                    vals[-1].div.append(qs_para)

        return table_soup

    @staticmethod
    def create_summary_style_string():
        """
        create a string of style information for the summary page html
        :return: summary style info (str)
        """
        style = (
            "body {font-family: Arial, Helvetica, sans-serif;}\n"
            "h1 {text-align: center;}\n"
            "p {margin-top: 0em}\n"
            "table {width: 100%; table-layout: fixed; border-collapse: collapse; background-color: white;}\n"
            "th, td {padding: 4px; vertical-align: top;}\n"
            # specify width of index and gene name columns
            "table th:nth-child(1) {width: 2%;}\ntable th:nth-child(2) {width: 10%;}\n"
            # specify width of search terms and results columns
            "table th:nth-child(3) {width: 30%;}\ntable th:nth-child(4) {width: 8%;}\n" 
            # specify widths of title and abstract columns
            "table th:nth-child(5) {width: 10%;}\ntable th:nth-child(6) {width: 40%;}\n"
            # format table headings
            "th {font-weight: bold; text-align: left; background-color: #f2f2f2; white-space: nowrap;}\n"
            # add automatic scrolling to potentially long fields
            "td div {height: 3em; overflow: auto; padding: 4px;}\n"
            # alternate grey background on rows
            "tr:nth-child(even) {background-color: #f2f2f2;}\n"
            # highlight row that mouse is hovering over
            "tr:hover {background-color: #E0ECF8;}\n"
        )
        return style

    def make_html_page_for_summary(self, summary_df):
        """
        use the summary dataframe to create a html summary page
        :param summary_df: dataframe of summary information
        :return: html for summary page
        """
        # create a html table of the summary df
        raw_summary_html = summary_df.to_html(justify='center')

        # edit the html table
        table_soup = self.edit_summary_table_html(raw_summary_html)

        # wrap in html and body tags
        wrapped_table_soup = self.wrap_soup_in_body_and_html_tags(table_soup)

        # add a title to the soup, above the table
        title_tag = wrapped_table_soup.new_tag("h1")
        title_tag.string = "Results summary"
        wrapped_table_soup.body.insert(0, title_tag)

        # add style to the soup
        style = self.create_summary_style_string()
        final_soup = self.add_style_to_soup(wrapped_table_soup, style)

        # write a html file from the soup
        with open(self.summary_path, 'w+', encoding='utf-8') as f:
            f.write(final_soup.prettify(formatter='html'))

    def get_summary_path(self):
        """
        return the summary html file path as a string

        :return: str file path
        """
        return self.summary_path
