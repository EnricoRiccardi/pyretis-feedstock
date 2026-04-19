---
description: Custom Pydocstyle Settings
---
1. PyRETIS project enforces a strict 79-character line limit for all python code, including docstrings.
2. Docstring parameters use full words like `integer` or `object like pandas.DataFrame` rather than abbreviations like `int`.
3. Wrap long lines in docstrings to ensure no line exceeds 79 characters.
4. When writing `list or dict or object like pandas.DataFrame or object like pandas.Series` wrap it onto multiple lines correctly for Sphinx/NumPy style.
