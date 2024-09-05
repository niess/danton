try:
    import danton
except ModuleNotFoundError:
    # Import local build.
    from pathlib import Path
    import sys
    path = Path(__file__).parent.parent.parent / "src/python"
    sys.path.append(str(path))

project = "Danton"
copyright = "Université Clermont Auvergne, CNRS/IN2P3, LPCA"
author = "Valentin Niess"
release = "1.4.0"

highlight_language = "python3"

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosectionlabel",
    "sphinx.ext.intersphinx",
    'sphinx.ext.mathjax'
]

numfig = True

autodoc_member_order = "groupwise"
autosectionlabel_prefix_document = True
intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable/", None)
}

templates_path = ["_templates"]
exclude_patterns = []

rst_prolog = """
.. |nbsp| unicode:: 0xA0
   :trim:

.. role:: python(code)
    :language: python
    :class: highlight
"""

html_theme = "sphinx_book_theme"
html_static_path = ["_static"]
