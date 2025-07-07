#!/usr/bin/env python3
# This file is part of `pyfraglib`, a software suite to calculate fragmentomics
# features from cfDNA and perform downstream analyses.
#
# Copyright (C) 2024 Daniel Schütte, daniel.schuette@iccb-cologne.org
#
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version. This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
# more details. You should have received a copy of the GNU General Public
# License along with this program. If not, see <https://www.gnu.org/licenses/>.
import sys
import unittest
import os
import argparse


def discover_tests(
    test_dir: str = "tests", pattern: str = "test_*.py"
) -> unittest.TestSuite:
    """
    Discover and load all test modules in the test directory.

    Args:
        test_dir: Directory containing test modules
        pattern: Pattern to match test files

    Returns:
        TestSuite containing all discovered tests
    """
    loader = unittest.TestLoader()
    start_dir = os.path.join(os.path.dirname(__file__), "..", test_dir)
    suite = loader.discover(start_dir, pattern=pattern)
    return suite


def run_specific_tests(test_modules: list[str]) -> unittest.TestResult:
    """
    Run specific test modules.

    Args:
        test_modules: List of test module names to run

    Returns:
        TestResult object with test execution results
    """
    loader = unittest.TestLoader()
    suite = unittest.TestSuite()

    for module_name in test_modules:
        try:
            module_suite = loader.loadTestsFromName(f"tests.{module_name}")
            suite.addTest(module_suite)
        except ModuleNotFoundError:
            print(f"Warning: Test module '{module_name}' not found")

    runner = unittest.TextTestRunner(verbosity=2, buffer=True)
    return runner.run(suite)


def run_all_tests(verbosity: int = 2) -> unittest.TestResult:
    """
    Run all discovered tests.

    Args:
        verbosity: Level of test output detail

    Returns:
        TestResult object with test execution results
    """
    suite = discover_tests()
    runner = unittest.TextTestRunner(verbosity=verbosity, buffer=True)
    return runner.run(suite)


def create_argparser() -> argparse.ArgumentParser:
    """Create argument parser for test runner."""
    parser = argparse.ArgumentParser(
        description="PyFragLib test runner",
        epilog="Examples:\n"
               "  python test_runner.py               # Run all tests\n"
               "  python test_runner.py -m test_core  # Run specific module\n"
               "  python test_runner.py -v 1          # Less verbose output\n",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "-m", "--modules", nargs="+", dest="test_modules",
        help="Specific test modules to run (e.g., test_core test_fragment)"
    )
    parser.add_argument(
        "-v", "--verbosity", type=int, choices=[0, 1, 2], default=2,
        help="Test output verbosity level (0=quiet, 1=normal, 2=verbose)"
    )
    parser.add_argument(
        "--failfast", action="store_true",
        help="Stop testing after first failure"
    )
    parser.add_argument(
        "--list", action="store_true",
        help="List all available test modules"
    )
    return parser


def list_test_modules() -> None:
    """List all available test modules."""
    test_dir: str = os.path.join(os.path.dirname(__file__))
    test_files: list[str] = [
        f for f in os.listdir(test_dir)
        if f.startswith("test_") and f.endswith(".py")
    ]

    print("Available test modules:")
    for test_file in sorted(test_files):
        module_name = test_file[:-3]
        print(f"  {module_name}")


def print_test_summary(result: unittest.TestResult) -> None:
    """Print a summary of test results."""
    print("Test Summary")
    print("-"*30)

    total_tests = result.testsRun
    failures = len(result.failures)
    errors = len(result.errors)
    skipped = len(result.skipped) if hasattr(result, "skipped") else 0
    success_rate = ((total_tests - failures - errors) / total_tests * 100) \
        if total_tests > 0 else 0

    print(f"Tests run: {total_tests}")
    print(f"Failures: {failures}")
    print(f"Errors: {errors}")
    print(f"Skipped: {skipped}")
    print(f"Success rate: {success_rate:.1f}%")

    if result.wasSuccessful():
        print("\nAll tests passed.")
    else:
        print("\nSome tests failed.")
        if result.failures:
            print(f"\nFailures ({len(result.failures)}):")
            for test, traceback in result.failures:
                print(f"  - {test}")

        if result.errors:
            print(f"\nErrors ({len(result.errors)}):")
            for test, traceback in result.errors:
                print(f"  - {test}")


def main() -> int:
    """Main entry point for test runner."""
    project_root = os.path.join(os.path.dirname(__file__), "..")
    sys.path.insert(0, os.path.abspath(project_root))

    parser = create_argparser()
    args = parser.parse_args()
    list_modules: bool = args.list
    if list_modules:
        list_test_modules()
        return 0

    failfast: bool = args.failfast
    if failfast:
        unittest.TextTestRunner.failfast = True

    try:
        test_modules: list[str] = args.test_modules
        result: unittest.TestResult
        if test_modules:
            modules: str = ", ".join(test_modules)
            print(f"Running test modules: {modules}")
            result = run_specific_tests(test_modules)
        else:
            print("Running all tests...")
            verbosity: int = args.verbosity
            result = run_all_tests(verbosity)

        print_test_summary(result)
        return 0 if result.wasSuccessful() else 1

    except KeyboardInterrupt:
        print("\n\nTest execution interrupted by user")
        return 130
    except Exception as e:
        print(f"\nError during test execution: {e}")
        return 1


if __name__ == "__main__":
    sys.exit(main())
