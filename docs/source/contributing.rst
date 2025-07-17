Contributing
============

We welcome contributions to pyfraglib! This document outlines how to contribute to the project.

Getting Started
---------------

See :doc:`installation` for installation instructions. After initial setup, please create a feature branch and follow our coding standards outlined below. Please make sure that your commit messages roughly follow our style. After every commit, use

   .. code-block:: bash

      ./tools/dev_install.sh

to run our unit test suite, linting, and type checking. If all tests pass, you can create a pull request and we will merge your code after a quick review.


Coding Standards
----------------

Code Style
~~~~~~~~~~

* Follow PEP 8 style guidelines
* Use double quotes for strings (``"`` not ``'``)
* Maximum line length: < 80 characters (should be enforced by the linter)
* Use descriptive variable and function names (snake_case for functions, methods, and variables; CamelCase for class names)

Type Annotations
~~~~~~~~~~~~~~~~

* All functions must have complete type annotations
* Use ``typing.Final`` for constants
* Use ``typing.NoReturn`` for functions that don't return
* Follow strict mypy settings as defined in ``pyproject.toml``
* use ``# type: ignore`` very, very sparingly

Documentation and Testing
~~~~~~~~~~~~~~~~~~~~~~~~~

* All public functions and classes must have docstrings
* Use Google-style docstrings for consistency
* Include parameter descriptions and return value information
* Add examples for complex functions
* Update relevant ``.rst`` files for new features
* Update API documentation if adding new modules
* Write unit tests for all new functionality
* Use the existing test fixtures in ``tests/test_fixtures.py``
* Follow the naming convention ``test_<functionality>``
* Include edge cases and error conditions

Thank you for contributing to pyfraglib!
