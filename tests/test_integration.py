import pytest
import subprocess
import shutil


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
