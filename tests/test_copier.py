import pytest

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
    "Snakefile",
    "config/config.yaml",
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
    assert '# my_project' in content
    assert 'my_project_description' in content


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


@pytest.mark.parametrize("workshop_tutorial", [True, False])
def test_tutorial_parameter_file_tutorial(
    copie, copier_project_defaults, workshop_tutorial):
    # create project
    copier_project_defaults["workshop_tutorial"] = workshop_tutorial
    project = copie.copy(extra_answers=copier_project_defaults)

    # test if tutorial.md file exists or not
    file_name = "tutorial.md"
    assert workshop_tutorial == project.project_dir.joinpath(
        file_name).exists()


@pytest.mark.parametrize("workshop_tutorial", [True, False])
def test_tutorial_parameter_file_CHANGELOG(
    copie, copier_project_defaults, workshop_tutorial):
    # create project
    copier_project_defaults["workshop_tutorial"] = workshop_tutorial
    project = copie.copy(extra_answers=copier_project_defaults)

    # test if Changelog file exists or not
    file_name = "CHANGELOG.md"
    assert workshop_tutorial != project.project_dir.joinpath(
        file_name).exists()


@pytest.mark.parametrize("workshop_tutorial", [True, False])
def test_tutorial_parameter_file_environment_content(
    copie, copier_project_defaults, workshop_tutorial):
    # create project
    copier_project_defaults["workshop_tutorial"] = workshop_tutorial
    project = copie.copy(extra_answers=copier_project_defaults)
    # test environment file content
    content = project.project_dir.joinpath(
        "workflow", "envs", "environment.yaml").read_text()
    if workshop_tutorial:
        assert '  - pip' in content
        assert '  - numpy' not in content
        assert '    - pyrato' not in content
    else:
        assert '  - numpy' in content
        assert '    - pyrato' in content
