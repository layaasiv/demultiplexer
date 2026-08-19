#!/usr/bin/env python

import gzip

def get_fastq_length(fastq_file: str) -> int:
    """
    Given path to FASTQ file, return the number of lines in the file

    Input:
        fastq_file - str: Path to the FASTQ file.
    Ouput:
        int: The number of lines in the FASTQ file.
    """
    with gzip.open(fastq_file, "r") as fq:
        line_count = sum(1 for line in fq)
        return line_count

def verify_equal_fastq_lengths(r1_len: int, r2_len: int, r3_len: int, r4_len: int) -> bool:
    """
    Check if all FASTQ files contain the same number of lines.

    Input:
        r1_len - int: Number of lines in read 1 FASTQ file.
        r2_len - int: Number of lines in read 2 FASTQ file.
        r3_len - int: Number of lines in read 3 FASTQ file.
        r4_len - int: Number of lines in read 4 FASTQ file.
    Output:
        bool: Whether all input values are equal or not.
    """
    return r1_len == r2_len == r3_len == r4_len

def verify_complete_records(fastq_len: int) -> bool:
    """
    Determine whether the length of the FASTQ is divisible by 4, indicating all the records it contains are complete.

    Input:
        fastq_len - int: Length of the FASTQ file.
    Output:
        bool: True is divisible by 4, False otherwise.
    """
    return isinstance(fastq_len/4, int)

def count_headers(fastq_file: str) -> int:
    """
    Count the number of header lines in the input FASTQ file.

    Input:
        fastq_file - str: Path to FASTQ file.
    Ouput:
        int: The number of headers in the file.
    """
    header_counter = 0
    with gzip.open(fastq_file, "r") as fq:
        while True:
            line = fq.readline()
            if line == "":
                break
            if line.startswith("@"):
                header_counter += 1

def verify_seqlen_equal_qscore(fastq_file: str) -> bool:
    """
    Checks whether sequence and quality score lines are equal in all records of the FASTQ file.

    Input:
        fastq_file - str: Path to the FASTQ file.
    Output:
        bool: False if unequal record encountered, True otherwise.
    """
    with gzip.open(fastq_file, "r") as fq:
        while True:
            line = fq.readline()
            if line == "":
                break
            elif line.startswith("@"):
                seq = fq.readline().strip("\n")
                fq.readline()
                qscores = fq.readline().strip("\n")
                if len(seq) != len(qscores):
                    print("Some records in the FASTQ file have sequence and quality score lines of differing lengths.")
                    return False
    return True