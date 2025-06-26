# Copyright (c) 2025.
# Authors: Tamás K. Stenczel
# This program is distributed under the MIT License, see LICENSE.md.

"""Add or update copyright notice to all python files."""

from datetime import date
import os
import re
import subprocess

# git users -> real name mapping for current contributors
USER_TO_REAL_NAME = {
    "stenczelt": "Tamás K. Stenczel",
    "Tamas K Stenczel": "Tamás K. Stenczel",
    "Tamas Stenczel": "Tamás K. Stenczel",
    "T. K. Stenczel": "Tamás K. Stenczel",
    "Adam Fekete": "Ádám Fekete",
    "ElliottKasoar": "Elliott Kasoar",
    "gelzinyte": "Elena Gelzinyte",
    "gabor1": "Gábor Csányi",
    "gabor": "Gábor Csányi",
}


def get_contributors(file_path):
    """Get list of contributors for a given file."""
    command = f"git log --pretty=format:'%an' {file_path} | sort | uniq"
    contributors = subprocess.check_output(command, shell=True, text=True)
    contributors = contributors.strip().split("\n")

    # Map usernames to real names, if possible
    real_name_contributors = []
    for contributor in contributors:
        if contributor == "":
            continue
        real_name = USER_TO_REAL_NAME[contributor]
        if real_name not in real_name_contributors:
            real_name_contributors.append(real_name)

    return real_name_contributors


def insert_copyright(file_path, contributors):
    """Insert a copyright notice at the top of a file."""
    # Create the copyright notice string
    author_names = ", ".join(contributors)
    copyright_notice = (
        f"# Copyright (c) {date.today().year}.\n"
        f"# Authors: {author_names}\n"
        f"# This program is distributed under the MIT License, see LICENSE.md.\n\n"
    )

    with open(file_path) as file:
        content = file.read()

    # replace existing notice if needed
    pattern = re.compile(
        r"^# Copyright \(c\) [\d-]+\.\n"
        r"# Authors:\s*(.*?)\n"
        r"# This program is distributed under the MIT License, see LICENSE\.md\.\n\n",
        re.MULTILINE,
    )

    if m := pattern.match(content):
        if m.group(1) == author_names:
            print(file_path, "Matched! contributors unchanged")
        else:
            print(file_path, "Updating contributors: ", m.group(1), " ->", author_names)
        content = pattern.sub(copyright_notice, content)
    else:
        print(file_path, "adding contributors", author_names)
        content = copyright_notice + content

    # Insert the copyright notice at the top
    with open(file_path, "w") as file:
        file.write(content)


def process_files(directory):
    """Process all Python files in the given directory."""
    for root, _, files in os.walk(directory):
        if root.startswith("./.venv"):
            continue
        for file in files:
            if file.endswith(".py") and file != "__init__.py":
                file_path = os.path.join(root, file)
                contributors = get_contributors(file_path)
                insert_copyright(file_path, contributors)


if __name__ == "__main__":
    # Replace with your folder path
    process_files(".")
