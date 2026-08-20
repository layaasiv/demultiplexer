#!/usr/bin/env python
from pathlib import Path
import gzip

indexes = [line.split("\t")[1] for line in [line.strip() for line in open("data/indexes.txt", "r").readlines()]]

# check file existance
def check_file_exists(file_path: str) -> bool:
    """
    Check if the specified file exists.

    Input:
        file_path (str): The path to the file to check.
    Output:
        bool: True if the file exists, False otherwise.
    """
    return Path(file_path).is_file()

# number of headers in a file
def count_headers(fastq_file: str) -> int:
    """
    Count the number of header lines in the input FASTQ file.

    Input:
        fastq_file - str: Path to FASTQ file.
    Ouput:
        int: The number of headers in the file.
    """
    header_count = 0
    with gzip.open(fastq_file, "rt") as fq:
        while True:
            line = fq.readline()
            if line == "":
                break
            if line.startswith("@"):
                header_count += 1
    return header_count

# number of lines in a file
def get_fastq_length(fastq_file: str) -> int:
    """
    Count the number of lines in the input FASTQ file.

    Input:
        fastq_file - str: Path to FASTQ file.
    Ouput:
        int: The number of lines in the file.
    """
    with gzip.open(fastq_file, "rt") as fq:
        line_count = sum(1 for line in fq)
        return line_count

# verify indexes
def unknown_indexes(fastq_file: str) -> bool:
    """
    Check if the indexes in the input FASTQ file are unknown (i.e., not specified in the indexes file, or contain 'N' bases).

    Input:
        fastq_file - str: Path to FASTQ file.
    Ouput:
        bool: True if unknown indexes are found, False otherwise.
    """
    with gzip.open(fastq_file, "rt") as fq:
        while True:
            line = fq.readline()
            if line == "":
                break
            if line.startswith("@"):
                indexes = line.strip("\n").split(" ")[-1]
                if "N" not in indexes:
                    return False
                elif indexes.split("-")[0] not in indexes or indexes.split("-")[1] not in indexes:
                    return False
    return True

def hopped_indexes(fastq_file: str) -> bool:
    """
    Check if the indexes in the input FASTQ file are hopped (i.e., the two index sequences do not match).

    Input:
        fastq_file - str: Path to FASTQ file.
    Output:
        bool: True if hopped indexes are found, False otherwise.
    """
    with gzip.open(fastq_file, "rt") as fq:
        while True:
            line = fq.readline()
            if line == "":
                break
            if line.startswith("@"):
                indexes = line.strip("\n").split(" ")[-1]
                if indexes.split("-")[0] == indexes.split("-")[1]:
                    return False
                elif indexes.split("-")[0] not in indexes or indexes.split("-")[1] not in indexes or "N" in indexes:
                    return False
    return True

def matched_indexes(fastq_file: str) -> bool:
    """
    Check if the indexes in the input FASTQ file are matched (i.e., the two index sequences match).

    Input:
        fastq_file - str: Path to FASTQ file.
    Output:
        bool: True if matched indexes are found, False otherwise.
    """
    with gzip.open(fastq_file, "rt") as fq:
        while True:
            line = fq.readline()
            if line == "":
                break
            elif line.startswith("@"):
                indexes = line.strip().split(" ")[-1]
                if indexes.split("-")[0] != indexes.split("-")[1]:
                    return False
                elif indexes.split("-")[0] not in indexes or indexes.split("-")[1] not in indexes or "N" in indexes:
                    return False
    return True

def test_dmux(output_path:str) -> bool:
    """
    Test the demultiplexing process by checking the output files for expected properties.

    Output:
        None
    """
    # Check if output files exist
    for index in indexes:
        assert check_file_exists(f"{output_path}/{index}_R1.fastq.gz"), f"Output file for {index} R1 does not exist."
        assert check_file_exists(f"{output_path}/{index}_R2.fastq.gz"), f"Output file for {index} R2 does not exist."

    # Check if unknown and hopped files exist
    assert check_file_exists(f"{output_path}/unknown_R1.fastq.gz"), "Unknown R1 output file does not exist."
    assert check_file_exists(f"{output_path}/unknown_R2.fastq.gz"), "Unknown R2 output file does not exist."
    assert check_file_exists(f"{output_path}/hopped_R1.fastq.gz"), "Hopped R1 output file does not exist."
    assert check_file_exists(f"{output_path}/hopped_R2.fastq.gz"), "Hopped R2 output file does not exist."

    # Check if the number of headers matches the expected count
    for index in indexes:
        r1_header_count = count_headers(f"{output_path}/{index}_R1.fastq.gz")
        r2_header_count = count_headers(f"{output_path}/{index}_R2.fastq.gz")
        assert r1_header_count == r2_header_count, f"Header count mismatch for {index}: R1 has {r1_header_count}, R2 has {r2_header_count}."

    # Check if the number of lines matches the expected count
    for index in indexes:
        r1_line_count = get_fastq_length(f"{output_path}/{index}_R1.fastq.gz")
        r2_line_count = get_fastq_length(f"{output_path}/{index}_R2.fastq.gz")
        assert r1_line_count == r2_line_count, f"Line count mismatch for {index}: R1 has {r1_line_count}, R2 has {r2_line_count}."

    # Check matched, hopped, and unknown indexes are correctly identified
    for index in indexes:
        r1_line_count = get_fastq_length(f"{output_path}/{index}_R1.fastq.gz")
        r2_line_count = get_fastq_length(f"{output_path}/{index}_R2.fastq.gz")
        if r1_line_count > 0 and r2_line_count > 0:
            assert matched_indexes(f"{output_path}/{index}_R1.fastq.gz"), f"Matched indexes not found for {index} R1."
            assert matched_indexes(f"{output_path}/{index}_R2.fastq.gz"), f"Matched indexes not found for {index} R2."

    assert unknown_indexes(f"{output_path}/unknown_R1.fastq.gz"), "Unknown indexes not found in unknown R1 output."
    assert unknown_indexes(f"{output_path}/unknown_R2.fastq.gz"), "Unknown indexes not found in unknown R2 output."
    assert hopped_indexes(f"{output_path}/hopped_R1.fastq.gz"), "Hopped indexes not found in hopped R1 output."
    assert hopped_indexes(f"{output_path}/hopped_R2.fastq.gz"), "Hopped indexes not found in hopped R2 output."