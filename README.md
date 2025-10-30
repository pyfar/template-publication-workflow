# Template: Publication workflow

A simple template intended for data processing workflows alongside scientific publications.

This template can be used with [copier](https://copier.readthedocs.io) to initialize a
new project structure.
Workflows are implemented using [snakemake](https://snakemake.readthedocs.io/en/stable/).

# Getting started

To use this template:

1. Install copier using pip within a virtual environment

    ```console
    pip install copier
    ```

2. Run the below copier command to create a target directory at `path_to_target_folder`.
   Note: copier will overwrite files in the given directory, if it already exists.

    ```console
    copier copy https://github.com/pyfar/template-publication-workflow path_to_target_folder
    ```

3. Open the *README.md* in the created project directory to get started with your project.


# Contributing

Check out [pyfar.org](https://pyfar-gallery.readthedocs.io/en/latest/contribute/index.html) if you want to become part of pyfar.