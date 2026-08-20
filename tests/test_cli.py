#!/usr/bin/env python

import subprocess
from pathlib import Path
import gzip

def count_fastq_records(fastq_file: str) -> int:
    """
    Count the number of records in a FASTQ file.

    Input:
        fastq_file (str): Path to the FASTQ file.
    Output:
        int: The number of records in the FASTQ file.
    """
    with gzip.open(fastq_file, "rt") as fq:
        return sum(1 for _ in fq) // 4

def get_record_ids(fastq_file: str) -> list:
    """
    Get the record IDs from a FASTQ file.

    Input:
        fastq_file (str): Path to the FASTQ file.
    Output:
        list: A list of record IDs from the FASTQ file.
    """
    ids = []
    with gzip.open(fastq_file, "rt") as fq:
        while True:
            line = fq.readline()
            if line == "":
                break
            if line.startswith("@"):
                ids.append(line.strip().split(":")[5])
    return ids

def cli_end_to_end_pipeline(tmp_path: str):
    """
    Run the CLI end-to-end pipeline and verify the output files.

    Input:
        tmp_path (str): Path to a temporary directory for output files.
    Output:
        None
    """
    output_dir = Path(tmp_path) / "output"

    # run the CLI command
    result = subprocess.run(
        [
            "python",
            "new_dmux.py",
            "-r1 tests/data/R1.fastq.gz",
            "-i1 tests/data/R2.fastq.gz",
            "-i2 tests/data/R3.fastq.gz",
            "-r2 tests/data/R4.fastq.gz",
            "-i tests/data/indexes.txt",
            "-o", str(output_dir),
        ],
        capture_output=True,
        text=True,
    )

    # check that the command ran successfully
    assert result.returncode == 0, f"CLI command failed with error: {result.stderr}"

    # --------------------- check output file lengths -----------------------
    # hopped records
    assert count_fastq_records(output_dir / "hopped_R1.fastq.gz") == 3, "Expected 3 hopped records in hopped_R1.fastq.gz."
    assert count_fastq_records(output_dir / "hopped_R2.fastq.gz") == count_fastq_records(output_dir / "hopped_R1.fastq.gz"), "Expected equal number of hopped records in hopped_R1.fastq.gz and hopped_R2.fastq.gz."

    # unknown records
    assert count_fastq_records(output_dir / "unknown_R1.fastq.gz") == 2, "Expected 2 unknown records in unknown_R1.fastq.gz."
    assert count_fastq_records(output_dir / "unknown_R2.fastq.gz") == count_fastq_records(output_dir / "unknown_R1.fastq.gz"), "Expected equal number of unknown records in unknown_R1.fastq.gz and unknown_R2.fastq.gz."

    # matched records
    assert count_fastq_records(output_dir / "GTAGCGTA_R1.fastq.gz") == 3, "Expected 3 matched records in GTAGCGTA_R1.fastq.gz."
    assert count_fastq_records(output_dir / "GTAGCGTA_R2.fastq.gz") == count_fastq_records(output_dir / "GTAGCGTA_R1.fastq.gz"), "Expected equal number of matched records in GTAGCGTA_R1.fastq.gz and GTAGCGTA_R2.fastq.gz."

    assert count_fastq_records(output_dir / "CGATCGAT_R1.fastq.gz") == 2, "Expected 2 matched records in CGATCGAT_R1.fastq.gz."
    assert count_fastq_records(output_dir / "CGATCGAT_R2.fastq.gz") == count_fastq_records(output_dir / "CGATCGAT_R1.fastq.gz"), "Expected equal number of matched records in CGATCGAT_R1.fastq.gz and CGATCGAT_R2.fastq.gz."

    for index in ["GATCAAGG", "AACAGCGA"]:
        assert count_fastq_records(output_dir / f"{index}_R1.fastq.gz") == 1, f"Expected 1 matched record in {index}_R1.fastq.gz."
        assert count_fastq_records(output_dir / f"{index}_R2.fastq.gz") == count_fastq_records(output_dir / f"{index}_R1.fastq.gz"), f"Expected equal number of matched records in {index}_R1.fastq.gz and {index}_R2.fastq.gz."
    
    for index in ["TAGCCATG", "CGGTAATC", "CTCTGGAT", "TACCGGAT", "CTAGCTCA", "CACTTCAC"]:
        assert count_fastq_records(output_dir / f"{index}_R1.fastq.gz") == 0, f"Expected 0 matched records in {index}_R1.fastq.gz."
        assert count_fastq_records(output_dir / f"{index}_R2.fastq.gz") == count_fastq_records(output_dir / f"{index}_R1.fastq.gz"), f"Expected equal number of matched records in {index}_R1.fastq.gz and {index}_R2.fastq.gz."

    # ---------- direct verification of classification of records -------------
    GTAGCGTA_ids = ['1265', '1682', '1775'] 
    CGATCGAT_ids = ['1286', '1721'] 
    GATCAAGG_ids = ['1347']
    AACAGCGA_ids = ['1367']
    hopped_ids = ['1401', '1450', '1512']
    unknown_ids = ['1574', '1620']
    matched_indexes = ["GTAGCGTA", "CGATCGAT", "GATCAAGG", "AACAGCGA"]

    # matched records
    for index in matched_indexes:
        r1_ids = get_record_ids(output_dir / f"{index}_R1.fastq.gz")
        r2_ids = get_record_ids(output_dir / f"{index}_R2.fastq.gz")
        assert r1_ids == r2_ids, f"Record IDs in {index}_R1.fastq.gz and {index}_R2.fastq.gz do not match."
        assert set(r1_ids) == set(eval(f"{index}_ids")), f"Record IDs in {index}_R1.fastq.gz do not match expected IDs."
    
    # hopped records
    hopped_r1_ids = get_record_ids(output_dir / "hopped_R1.fastq.gz")
    hopped_r2_ids = get_record_ids(output_dir / "hopped_R2.fastq.gz")
    assert hopped_r1_ids == hopped_r2_ids, "Record IDs in hopped_R1.fastq.gz and hopped_R2.fastq.gz do not match."
    assert set(hopped_r1_ids) == set(hopped_ids), "Record IDs in hopped_R1.fastq.gz do not match expected IDs."

    # unknown records
    unknown_r1_ids = get_record_ids(output_dir / "unknown_R1.fastq.gz")
    unknown_r2_ids = get_record_ids(output_dir / "unknown_R2.fastq.gz")
    assert unknown_r1_ids == unknown_r2_ids, "Record IDs in unknown_R1.fastq.gz and unknown_R2.fastq.gz do not match."
    assert set(unknown_r1_ids) == set(unknown_ids), "Record IDs in unknown_R1.fastq.gz do not match expected IDs."
