#!/bin/python3

import random
import gzip
import argparse
import os

random.seed(0)

fastq_format = """@{i}
{sequence}
+
{phred}"""


def get_reference(fasta):
    with open(fasta) as f:
        # strip the first line (contains only >chrM)
        return "".join(map(lambda s: s.replace("\n",""), f.readlines()[1:]))


def get_random_substring(string, length=10):
    start = random.randrange(len(string)-length)
    return string[start:start+length]


def get_fastqs(ref, fq_len=4, seq_len=100, read_len=10):
    fq = [[], []]
    for i in range(fq_len):
        seq = get_random_substring(ref, length=seq_len)
        phred = "K"*read_len
        fq[0].append(fastq_format.format(
            i=i, sequence=seq[:read_len], phred=phred))
        fq[1].append(fastq_format.format(
            i=i, sequence=seq[-read_len:], phred=phred))

    return list(map(lambda f: "\n".join(f), fq))


def generate_samples(ref, outputdir, nsamples=10, nsplit=2, **kwargs):
    if not os.path.exists(outputdir):
        os.makedirs(outputdir)

    for s in range(nsamples):
        fq = get_fastqs(ref, **kwargs)
        for r in [1, 2]:
            for split in range(nsplit):
                fname = "{outputdir}/{s:04}_S1_L001_R{r}_{split:03}.fastq.gz".format(
                    outputdir=outputdir, s=s, r=r, split=split)
                with gzip.open(fname, "wb") as f:
                    f.write((fq[r-1]+"\n").encode('ascii'))


def main(args):
    generate_samples(
        ref=get_reference(args.fasta),
        outputdir=args.outputdir,
        nsamples=args.nsamples,
        nsplit=args.nsplit,
        fq_len=args.fq_len,
        seq_len=args.seq_len,
        read_len=args.read_len
    )


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--nsamples', default=10, type=int)
    parser.add_argument('--nsplit', default=2, type=int)
    parser.add_argument('--fq_len', default=4, type=int)
    parser.add_argument('--seq_len', default=100, type=int)
    parser.add_argument('--read_len', default=10, type=int)
    parser.add_argument('fasta', type=str)
    parser.add_argument('outputdir', type=str)
    args = parser.parse_args()
    main(args)
