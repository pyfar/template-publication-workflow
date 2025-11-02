import pytest

@pytest.fixture(scope='session')
def copier_project_defaults():
    return {
        "project_name": "my_project",
        "project_description": "my_project_description",
        "author_name": "author",
        "python_version": "3.10",
        "license": "MIT",
        "version": "0.1.0",
        "year": "2024",
        }
