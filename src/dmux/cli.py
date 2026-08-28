#!/usr/bin/env python

import argparse

from dmux.dmux import demultiplex

__version__ = "0.1.0"


def get_args():
    parser = argparse.ArgumentParser(
        description="""
    A program to demultiplex a set of paired-end, dual-indexed FASTQ files. Reads are classified as matched, index-hopped, or unknown.
    The program takes as input the four FASTQ files (R1, R2, I1, I2) and a text file containing the list of known indexes. The output consists of separate FASTQ files for each matched index, as well as files for unknown and hopped reads.
    """
    )
    parser.add_argument("--version", action="version", version=f"dmux {__version__}")
    parser.add_argument(
        "-i",
        "--indexes",
        help="Text file containing list of known indexes",
        required=True,
    )
    parser.add_argument(
        "-r1", "--read1", help="Zipped read 1 FASTQ file (forward read)", required=True
    )
    parser.add_argument(
        "-i1", "--index1", help="Zipped index 1 FASTQ file", required=True
    )
    parser.add_argument(
        "-i2", "--index2", help="Zipped index 2 FASTQ file", required=True
    )
    parser.add_argument(
        "-r2", "--read2", help="Zipped read 2 FASTQ file (reverse read)", required=True
    )
    parser.add_argument(
        "-o",
        "--outputpath",
        help="Path to output directory",
        required=False,
        default="./",
    )
    return parser.parse_args()


def main():
    args = get_args()

    demultiplex(
        index_file=args.indexes,
        r1_file=args.read1,
        i1_file=args.index1,
        i2_file=args.index2,
        r2_file=args.read2,
        output_path=args.outputpath,
    )


if __name__ == "__main__":
    main()
