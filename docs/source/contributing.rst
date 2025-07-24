Contributing
============

We welcome contributions to ``pyfraglib``! This document outlines how to contribute to the project.

Getting Started
---------------

See :doc:`installation` for installation instructions. After the initial setup, please create a feature branch and follow our coding standards outlined below. Please make sure that your commit messages roughly follow our style. After every commit, use

   .. code-block:: bash

      ./tools/dev_install.sh

to run our unit test suite, linting, and type checking. If all tests pass, you can create a pull request and we will merge your code after a quick review.


Coding Standards
----------------

Code Style
~~~~~~~~~~

* Follow PEP 8 style guidelines
* Use double quotes for strings (``"`` not ``'``)
* Maximum line length: < 80 characters (this is also enforced by the linter)
* Use descriptive variable and function names (``snake_case`` for functions, methods, and variables; ``CamelCase`` for class names)

Type Annotations
~~~~~~~~~~~~~~~~

All functions must have complete type annotations! Please use ``typing.Final`` for constants and ``typing.NoReturn`` for functions that don't return to the caller. We are trying to follow very strict mypy settings as defined in ``pyproject.toml``, but because some Python libraries are not typed and others rely on dynamic typing a lot, we sometimes need to disable ``mypy``. Nonetheless, please use ``# type: ignore`` very, very sparingly

Documentation and Testing
~~~~~~~~~~~~~~~~~~~~~~~~~

* All public functions and classes must have docstrings
* Include parameter descriptions and return value information
* Add examples for complex functions
* Update relevant ``.rst`` files when implementing new features
* Write as many unit tests for new functionality as possible
* Use the existing test fixtures in ``tests/test_fixtures.py`` and Follow the naming convention ``test_<functionality>``

Thank you for contributing to ``pyfraglib``!
