import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="litspy",
    version="0.0.6.2",
    author="E Croot",
    author_email="ec339@le.ac.uk",
    description="Searches through all titles and abstracts available in Europe PMC for co-occurrence of supplied terms "
                "(and their synonyms where available), and produces summaries of the search results.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/ec339/litspy/tree/master/litspy",
    packages=setuptools.find_packages(),
    license='LICENCE.txt',
    classifiers=[
        "Programming Language :: Python",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    install_requires=[
        'beautifulsoup4',
        'requests',
        'lxml',
        'wordcloud',
        'pandas',
        'matplotlib',
        'xlrd',
        'openpyxl',
        'rtgo'
    ],
    scripts=['litspy/__main__.py']
)