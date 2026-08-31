#!/usr/bin/env python

import argparse
import sys
from pathlib import Path

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

def validate_inputs(args):
    files = {
        "R1": args.read1,
        "I1": args.index1,
        "I2": args.index2,
        "R2": args.read2,
        "indexes": args.indexes
    }

    for name, filepath in files.items():
        if not Path(filepath).is_file():
            raise FileNotFoundError(
                f"{name} does not exist at: {filepath}."
            )

def main():
    args = get_args()

    try:
        validate_inputs(args)

        demultiplex(
            index_file=args.indexes,
            r1_file=args.read1,
            i1_file=args.index1,
            i2_file=args.index2,
            r2_file=args.read2,
            output_path=args.outputpath,
        )
    
    except FileNotFoundError as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1
    
    except ValueError as e:
        print(f"Error: {e}", file=sys.stderr)
        return 2
    
    return 0


if __name__ == "__main__":
    sys.exit(main())
