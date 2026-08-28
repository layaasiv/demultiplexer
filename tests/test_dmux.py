from pathlib import Path
import gzip
from dmux.dmux import demultiplex

DATA_DIR = Path(__file__).parent / "data"
INDEX_FILE = DATA_DIR / "indexes.txt"
with open(INDEX_FILE, "r") as fh:
    indexes = [
        line.strip().split("\t")[1]
        for line in fh
    ]

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
    Check whether all records in a FASTQ file have unknown indexes.

    Returns:
        bool: True if all indexes are unknown, False otherwise.
    """
    with gzip.open(fastq_file, "rt") as fq:
        for line in fq:
            if line.startswith("@"):
                index_pair = line.strip().split(" ")[-1]
                index1, index2 = index_pair.split("-")

                if (
                    "N" in index1
                    or "N" in index2
                    or index1 not in indexes
                    or index2 not in indexes
                ):
                    continue

                # Found a known, non-N index
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
        for line in fq:
            if line.startswith("@"):
                index_pair = line.strip().split(" ")[-1]
                index1, index2 = index_pair.split("-")
                if index1 != index2:
                    continue
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
        for line in fq:
            if line.startswith("@"):
                index_pair = line.strip().split(" ")[-1]
                index1, index2 = index_pair.split("-")
                if (
                    index1 == index2
                    and index1 in indexes
                    and index2 in indexes
                    and "N" not in index1
                    and "N" not in index2
                ):
                    continue
                return False
    return True

def test_dmux(tmp_path: str) -> bool:
    """
    Test the demultiplexing process by checking the output files for expected properties.

    Output:
        None
    """

    output_path = Path(tmp_path) / "output"

    demultiplex(
        index_file=INDEX_FILE,
        r1_file=DATA_DIR / "R1.fastq.gz",
        i1_file=DATA_DIR / "R2.fastq.gz",
        i2_file=DATA_DIR / "R3.fastq.gz",
        r2_file=DATA_DIR / "R4.fastq.gz",
        output_path=output_path,
    )

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