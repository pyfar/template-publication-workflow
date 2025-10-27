import os
import pytest

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


@pytest.mark.parametrize("file_name", [
    "README.md",
    "LICENSE",
    "CHANGELOG.md",
    ".gitignore",
    "workflow/envs/environment.yaml",
    "Snakefile",
    "workflow/Snakefile",
])
def test_generated_file_exists(copie, copier_project_defaults, file_name):
    # create project
    project_defaults = copier_project_defaults
    project = copie.copy(extra_answers=project_defaults)

    # test generated file
    assert os.path.exists(os.path.join(project.project_dir, file_name))


def test_readme(copie, copier_project_defaults):
    # create project
    project_defaults = copier_project_defaults
    project = copie.copy(extra_answers=project_defaults)

    # test README file content
    content = project.project_dir.joinpath("README.md").read_text()
    assert '# Welcome to my_project' in content


def test_license_default(copie, copier_project_defaults):
    # create project
    project_defaults = copier_project_defaults
    project = copie.copy(extra_answers=project_defaults)

    # test LICENSE file content
    content = project.project_dir.joinpath("LICENSE").read_text()
    assert 'MIT License' in content
    assert '2024 author' in content


def test_environment_file(copie, copier_project_defaults):
    # create project
    project_defaults = copier_project_defaults
    project = copie.copy(extra_answers=project_defaults)

    # test environment file content
    content = project.project_dir.joinpath(
        "workflow", "envs", "environment.yaml").read_text()
    assert 'name: my_project' in content
    assert 'python=3.10' in content
