import os
import subprocess
import sys
from sys import platform
from typing import Sequence

from copier import run_copy
import pytest

IS_WINDOWS = platform.startswith('win')

@pytest.fixture
def defaults():
    return {
        "project_name": "my_project",
        "package_description": "my_project",
        "author_name": "author",
        "author_email": "author@example.com",
        "python_version": "3.10",
        "license": "MIT",
        "username": "pyfar",
        "publication_slug": "year_where_shorttitle",
        "version": "0.1.0",
    }

def test_copier_defaults(tmpdir, defaults):
    run_copy('.', tmpdir, defaults)
    assert (tmpdir / 'pyproject.toml').exists()


