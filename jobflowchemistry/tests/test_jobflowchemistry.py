"""
Unit and regression test for the jobflowchemistry package.
"""

# Import package, test suite, and other packages as needed
import sys

import pytest

import jobflowchemistry


def test_jobflowchemistry_imported():
    """Sample test, will always pass so long as import statement worked."""
    assert "jobflowchemistry" in sys.modules
