import os
import pytest
import subprocess
import shutil
import json

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
])
def test_generated_file_exists(copie, copier_project_defaults, file_name):
    # create project
    project_defaults = copier_project_defaults
    project = copie.copy(extra_answers=project_defaults)

    # test generated file
    assert project.project_dir.joinpath(file_name).exists()


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


@pytest.mark.parametrize(('license_short', 'expected_content'), [
    ("MIT", "MIT License"),
    ("CC-BY-4.0", "Attribution 4.0 International"),
    ("CC-BY-NC-4.0", "Attribution-NonCommercial 4.0 International"),
    (
        "CC-BY-NC-ND-4.0",
        "Attribution-NonCommercial-NoDerivatives 4.0 International"),
    (
        "CC-BY-NC-SA-4.0",
        "Attribution-NonCommercial-ShareAlike 4.0 International"),
    ("CC-BY-ND-4.0", "Attribution-NoDerivatives 4.0 International"),
    ("CC-BY-SA-4.0", "Attribution-ShareAlike 4.0 International"),
])
def test_license_all_license_files(
        copie, copier_project_defaults, license_short, expected_content):
    # create project
    project_defaults = copier_project_defaults
    project_defaults["license"] = license_short
    project = copie.copy(extra_answers=project_defaults)

    # test LICENSE file content
    content = project.project_dir.joinpath("LICENSE").read_text()
    assert expected_content in content

def test_environment_file(copie, copier_project_defaults):
    # create project
    project_defaults = copier_project_defaults
    project = copie.copy(extra_answers=project_defaults)

    # test environment file content
    content = project.project_dir.joinpath(
        "workflow", "envs", "environment.yaml").read_text()
    assert 'name: my_project' in content
    assert 'python=3.10' in content


def test_conda_env_create_dry_run(copie, copier_project_defaults):
    # skip if conda is not available in the environment running the tests
    if shutil.which("conda") is None:
        pytest.skip(
            "conda executable not found on PATH; "
            "skipping conda validation test")

    # create copier project
    project = copie.copy(extra_answers=copier_project_defaults)
    env_path = project.project_dir.joinpath(
        "workflow", "envs", "environment.yaml")

    # use --dry-run so the test validates conda can parse/solve the
    # environment without creating it
    # if conda returns a non-zero exit code the test will fail and
    # include stderr for debugging
    try:
        subprocess.run(
            ["conda", "env", "create", "-f", str(env_path), "--dry-run"],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
    except subprocess.CalledProcessError as exc:
        pytest.fail(
            "conda failed to validate environment.yaml\n"
            f"exit code: {exc.returncode}\n"
            f"stdout: {exc.stdout}\n"
            f"stderr: {exc.stderr}",
        )


@pytest.mark.parametrize("add_zenodo", [True, False])
def test_add_zenodo_file(copie, copier_project_defaults, add_zenodo):
    # create project
    project_defaults = copier_project_defaults
    project_defaults["add_zenodo"] = add_zenodo
    project = copie.copy(extra_answers=project_defaults)

    # test .zenodo.json file existence
    assert project.project_dir.joinpath(".zenodo.json").exists() == add_zenodo
    assert project.project_dir.joinpath(
        '.github', 'workflows','next_steps_zenodo_issue.md',
        ).exists() == add_zenodo


def test_zenodo_file_content(copie, copier_project_defaults):
    # create project
    project_defaults = copier_project_defaults
    project_defaults["add_zenodo"] = True
    project = copie.copy(extra_answers=project_defaults)

    # test .zenodo.json file content
    content = json.load(project.project_dir.joinpath(".zenodo.json").open())
    assert content["title"] == "author"
    assert content["version"] == "0.1.0"
    assert content["license"] == "MIT"
    assert content["creators"][0]["affiliation"] == "author"
    assert content["creators"][0]["name"] == "author"