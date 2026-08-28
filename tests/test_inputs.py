import gzip
from pathlib import Path

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

def count_headers(fastq_file: str) -> int:
    """
    Count the number of header lines in the input FASTQ file.

    Input:
        fastq_file - str: Path to FASTQ file.
    Ouput:
        int: The number of headers in the file.
    """
    header_counter = 0
    with gzip.open(fastq_file, "rt") as fq:
        while True:
            line = fq.readline()
            if line == "":
                break
            if line.startswith("@"):
                header_counter += 1
    return header_counter

def verify_seqlen_equal_qscore(fastq_file: str) -> bool:
    """
    Checks whether sequence and quality score lines are equal in all records of the FASTQ file.

    Input:
        fastq_file - str: Path to the FASTQ file.
    Output:
        bool: False if unequal record encountered, True otherwise.
    """
    with gzip.open(fastq_file, "rt") as fq:
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

def test_input_files(input_path):
    input_path = Path(input_path)
    r1_len = get_fastq_length(input_path / "R1.fastq.gz")
    r2_len = get_fastq_length(input_path / "R2.fastq.gz")
    r3_len = get_fastq_length(input_path / "R3.fastq.gz")
    r4_len = get_fastq_length(input_path / "R4.fastq.gz")

    assert(r1_len != 0)
    assert(r1_len == r2_len == r3_len == r4_len)
    assert(r1_len % 4 == 0)

    r1_headers = count_headers(input_path / "R1.fastq.gz")
    r2_headers = count_headers(input_path / "R2.fastq.gz")
    r3_headers = count_headers(input_path / "R3.fastq.gz")
    r4_headers = count_headers(input_path / "R4.fastq.gz")

    assert(r1_headers == r2_headers == r3_headers == r4_headers)

    assert(verify_seqlen_equal_qscore(input_path / "R1.fastq.gz"))
    assert(verify_seqlen_equal_qscore(input_path / "R2.fastq.gz"))
    assert(verify_seqlen_equal_qscore(input_path / "R3.fastq.gz"))
    assert(verify_seqlen_equal_qscore(input_path / "R4.fastq.gz"))