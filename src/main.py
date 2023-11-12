# main.py

import argparse, yaml
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns
import profile, bootstrap

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", type=str, help="Input directory", required=True)
    parser.add_argument("-o", "--output", type=str, help="Output filename", default="./out_tree.png")
    parser.add_argument("-c", "--config", type=str, help="Config file", default="./config.yaml")

    args = parser.parse_args()

    return args


def read_settings(config_f: str):
    """
    Generate settings dictionary from config YAML.
    """
    with open(config_f, "r") as f:
        print(f"Reading config file {config_f} ...")
        settings = yaml.safe_load(f)

    return settings


def reconstruct_lineage_tree(input_d, output_f, config_f: str) -> None:
    """
    Reconstruct lineage tree from clonal recordings in a directory.
    """
    # Get all FASTQ files in directory
    input_d = Path(input_d)
    fastq_files = input_d.glob("*.fastq")
    # Convert to list of strings
    fastq_files = [str(fastq_file) for fastq_file in fastq_files]

    # Read settings from config file
    settings = read_settings(config_f)
    spacer = settings["spacer"]
    h1, h2 = settings["h1"], settings["h2"]
    s1, s2 = settings["s1"], settings["s2"]
    l = settings["l"]
    D = settings["D"]
    parity = settings["parity"]

    # Retrieve and encode clonal recordings
    clones = []
    names = []
    for fastq_file in fastq_files:
        fastq_name = fastq_file.split("/")[-1].split(".")[0]
        names.append(fastq_name)
        print(f"Profiling {fastq_name} ...")
        clone = profile.get_recordings(fastq_file, spacer, h1, h2, s1, s2, l, D, parity)
        clones.append(clone)
    # Sort clones alpha-numerically by name
    clones = [x for _,x in sorted(zip(names, clones), key = lambda x: x[0])]

    # Bootstrap and visualize lineage tree
    print("Generating lineage matrix ...")
    A = bootstrap.generate_lineage_matrix(clones, parity)

    print("Bootstrapping and visualizing tree ...")
    plt.figure()
    sns.heatmap(A, cmap = "mako")
    plt.savefig(output_f, dpi=1000)

    print("Done.")


if __name__ == "__main__":
    print("@@@@@@@@@@@@@@@@@@@@ tree-builder @@@@@@@@@@@@@@@@@@@@")
    print("------------------------------------------------------")
    print("Reconstruct a lineage tree from clonal DNA recordings.")
    print("------------------------------------------------------")
    args = parse_args()
    reconstruct_lineage_tree(args.input, args.output, args.config)