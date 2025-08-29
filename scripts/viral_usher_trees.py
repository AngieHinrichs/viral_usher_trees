"""
Utility functions and constants for configuring and executing viral_usher tree builds.
"""

import os
import subprocess
import sys


trees_dir = "trees"
config_name = "config.toml"


def check_top_level_dir():
    """Make sure the expected trees directory exists, error out if it doesn't."""
    if not os.path.isdir(trees_dir):
        print("Can't find {trees_dir} directory.  This script must be run in top level of repo.", file=sys.stderr)
        sys.exit(1)


def generate_config(subdir_path: str, tree_name: str, refseq_acc: str, refseq_assembly: str, taxid: str) -> bool:
    """Make subdir_path/config.toml using viral_usher init with command-line args
    to skip the interactive process."""
    config_path = subdir_path + "/" + config_name
    config_path_tmp = config_path + ".tmp"
    command = ["viral_usher", "init",
               "--refseq", refseq_acc,
               "--taxonomy_id", taxid,
               # Nextclade dataset search is not working well enough, and command fails if there are multiple matches,
               # so just say no to nextclade by default; add back manually to config.toml files where applicable.
               "--nextclade_dataset", "",
               "--workdir", subdir_path,
               "--config", config_path_tmp]
    print(" ".join(command))
    try:
        subprocess.run(command, check=True)
        # Replace the absolute path for workdir
        with open(config_path_tmp, "r") as cfg_in, open(config_path, "w") as cfg_out:
            for line in cfg_in:
                if line.startswith("workdir"):
                    ix = line.index("/" + trees_dir + "/")
                    cfg_out.write("workdir = '." + line[ix:])
                else:
                    cfg_out.write(line)
        os.remove(config_path_tmp)
        return True
    except FileNotFoundError as e:
        print(f"viral_usher needs to be in your path. ({e})", file=sys.stderr)
        sys.exit(1)
    except subprocess.CalledProcessError:
        return False


def job_file_name(idx):
    """For now, just hardcode the output file name pattern and use the same pattern in the job script."""
    return f"job_list.{idx}.txt"
