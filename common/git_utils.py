"""Git utilities for QHAT tools."""
import os
import subprocess


def get_git_hash(reference_file=None):
    """
    Get current git commit hash with dirty flag if uncommitted changes exist.

    Parameters
    ----------
    reference_file : str, optional
        Path to a file in the repository. If provided, uses that file's directory
        as the repository root. If None, uses current working directory.

    Returns
    -------
    str
        Git hash, with "-dirty" suffix if uncommitted changes exist
    """
    if reference_file:
        dirpath = os.path.dirname(os.path.realpath(reference_file))
    else:
        dirpath = os.getcwd()

    commands = ";".join([
        f"pushd {dirpath} > /dev/null",
        "if [[ $(git diff --stat) != '' ]]",
        "then echo $(git rev-parse HEAD)-dirty",
        "else git rev-parse HEAD",
        "fi",
        "popd > /dev/null"
    ])
    output = subprocess.run(commands, shell=True, capture_output=True)
    return output.stdout.decode("utf-8")[:-1]
