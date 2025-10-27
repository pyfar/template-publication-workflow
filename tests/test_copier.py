import os
import subprocess
import sys
from sys import platform
from typing import Sequence

from copier import run_copy
import pytest

IS_WINDOWS = platform.startswith('win')


@pytest.fixture(scope='session')
def copier_project_defaults():
    return {
        "project_name": "my_project",
        "package_description": "my_project",
        "author_name": "author",
        "author_email": "author@example.com",
        "python_version": "3.10",
        "license": "MIT",
        "username": "pyfar",
        "version": "0.1.0",
        "year": "2024",
        }

def test_project_folder(copie, copier_project_defaults):
    project_defaults = copier_project_defaults
    project = copie.copy(extra_answers=project_defaults)

    assert project.exit_code == 0
    assert project.exception is None
    assert project.project_dir.is_dir()


def test_readme(copie, copier_project_defaults):
    # create project
    project_defaults = copier_project_defaults
    project = copie.copy(extra_answers=project_defaults)

    # test existence
    assert os.path.exists(os.path.join(project.project_dir, "README.md"))

    # test content
    content = project.project_dir.joinpath("README.md").read_text()
    assert '# Welcome to my_project' in content


def test_license_default(copie, copier_project_defaults):
    # create project
    project_defaults = copier_project_defaults
    project = copie.copy(extra_answers=project_defaults)

    # test LICENSE file
    assert os.path.exists(os.path.join(project.project_dir, "LICENSE"))
    content = project.project_dir.joinpath("LICENSE").read_text()
    assert 'MIT License' in content
    assert '2024 author' in content

