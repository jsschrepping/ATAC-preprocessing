"""This file contains helper python functions used in the Snakefile.

"""

from pathlib import Path


def get_samples(fastqdir):
    """Given a path `fastqdir` returns a dictionary `samples`, such that
    `samples[sample_name][1]` (or `2`) is a list of files containing
    R1 or R2 respectively.

    """

    if not Path(fastqdir).is_dir():
        raise Exception(f"""{fastqdir} does not exist""")

    files = Path(fastqdir).glob('**/*_R1_*.fastq.gz')
    samples = dict()

    for f in files:
        sample = f.name.split("_")[0]

        if sample not in samples:
            samples[sample] = {'1': [], '2': []}

        R1 = str(f)
        R2 = str(f).replace("_R1_", "_R2_")

        samples[sample]['1'].append(R1)
        samples[sample]['2'].append(R2)

    if not samples:
        raise Exception("Could not find any samples")

    return samples


def get_mode(samples):
    # check if read 2 exists for all files.  If at least one
    # R1 has no associated R2 switch to single end mode
    for sample in samples:
        for R2 in samples[sample]["2"]:
            if not Path(R2).is_file():
                return "single"
    return "paired"


def single_or_paired_input(wildcards, mode):
    """Selects between either a single or paired input based on the
mode"""

    sample = wildcards.sample
    R1, R2 = [f"""trimmed/{mode}/{sample}_R1.fastq""",
              f"""trimmed/{mode}/{sample}_R2.fastq"""]

    if mode == "single":
        return [R1]
    elif mode == "paired":
        return [R1, R2]
