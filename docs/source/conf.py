# Configuration file for the Sphinx documentation builder.

import os
import sys

# Add the project root to the Python path
sys.path.insert(0, os.path.abspath('../..'))

# Import version from package
try:
    from hbat import __version__
except ImportError:
    __version__ = '0.0.0+unknown'

# -- Project information -----------------------------------------------------

project = 'HBAT'
copyright = '2025, Abhishek Tiwari'
author = 'Abhishek Tiwari'

# The short X.Y version
version = __version__
# The full version, including alpha/beta/rc tags
release = __version__

# -- General configuration ---------------------------------------------------

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.viewcode',
    'sphinx.ext.napoleon',
    'sphinx.ext.intersphinx',
    'sphinx.ext.coverage',
    'sphinx.ext.githubpages',
    'sphinx_toolbox.shields',
    'sphinx_sitemap',
    'sphinx_copybutton',
    'sphinxcontrib.bibtex',
]

bibtex_bibfiles = ['hbat.bib']
bibtex_encoding = 'utf-8'
bibtex_default_style = 'unsrt'

# Mock imports for modules that might not be available in CI
autodoc_mock_imports = ['tkinter', 'matplotlib']

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
exclude_patterns = []

# The suffix(es) of source filenames.
source_suffix = '.rst'

# The master toctree document.
master_doc = 'index'

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.
html_theme = 'sphinx_book_theme'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

# Custom CSS files
html_css_files = [
    'custom.css',
]

# -- Extension configuration -------------------------------------------------

# autodoc configuration
autodoc_default_options = {
    'members': True,
    'member-order': 'bysource',
    'special-members': '__init__',
    'undoc-members': True,
    'exclude-members': '__weakref__'
}

# autosummary configuration
autosummary_generate = True

# Napoleon settings for Google/NumPy style docstrings
napoleon_google_docstring = True
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = False
napoleon_include_private_with_doc = False
napoleon_include_special_with_doc = True
napoleon_use_admonition_for_examples = False
napoleon_use_admonition_for_notes = False
napoleon_use_admonition_for_references = False
napoleon_use_ivar = False
napoleon_use_param = True
napoleon_use_rtype = True
napoleon_type_aliases = None

# Intersphinx mapping
intersphinx_mapping = {
    'python': ('https://docs.python.org/3/', None),
    'numpy': ('https://numpy.org/doc/stable/', None),
    'matplotlib': ('https://matplotlib.org/stable/', None),
}

# HTML theme options for sphinx_book_theme
html_theme_options = {
    "repository_url": "https://github.com/abhishektiwari/hbat",
    "repository_provider": "github",
    "repository_branch": "main",
    "path_to_docs": "docs/source",
    "use_issues_button": True,
    "use_repository_button": True,
    "use_edit_page_button": True, 
    "use_download_button": False,
    "use_fullscreen_button": True,
    "use_sidenotes": True,
    "icon_links_label": "Quick Links",
    "icon_links": [
        {
            "name": "GitHub",
            "url": "https://github.com/abhishektiwari/hbat",
            "icon": "fa-brands fa-square-github",
            "type": "fontawesome",
        },
        {
            "name": "Abhishek Tiwari",
            "url": "https://www.abhishek-tiwari.com",
            "icon": "https://www.abhishek-tiwari.com/images/logo.svg",
            "type": "local",
        },
        {
            "name": "PyPI",
            "url": "https://pypi.org/project/hbat/",
            "icon": "fa-brands fa-python",
            "type": "fontawesome",
        },
   ]
}

# The version info for the project you're documenting, acts as replacement for
# |version| and |release|, also used in various other places throughout the
# built documents.
html_title = f"{project} Latest Version"
html_short_title = f"{project} Latest Version"

# Add project logo
html_logo = '../../hbat.svg'
html_favicon = '../../hbat.ico'

# Output file base name for HTML help builder.
htmlhelp_basename = 'HBATdoc'
html_baseurl = 'https://hbat.abhishek-tiwari.com/'
sitemap_url_scheme = '{link}'

# -- Edit on GitHub configuration --------------------------------------------

# Variables for Edit on GitHub links
html_context = {
    'display_github': True,  # Enable GitHub integration
    'github_user': 'abhishektiwari',  # GitHub username
    'github_repo': 'hbat',  # GitHub repository name
    'github_version': 'main/',  # Git branch
    'conf_py_path': 'docs/source/',  # Path to documentation source
}

# -- Shield badges configuration ---------------------------------------------

# GitHub and PyPI information for badges
github_username = 'abhishektiwari'
github_repository = 'hbat'
pypi_name = 'hbat'