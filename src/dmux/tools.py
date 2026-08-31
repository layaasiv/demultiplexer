def fastq_parser(fastq_file: str) -> tuple:
    """
    Parses a FASTQ file one record at a time, and returns a tuple containing the sequence identifier, the sequence, and the quality scores.

    Input:
        fastq_file (str): The opened file handle to the FASTQ file to be parsed.
    Ouput:
        tuple: A tuple containing the sequence identifier, the sequence, and the quality scores.
    """
    while True:
        header = fastq_file.readline().rstrip("\n")
        if not header:
            break
        seq = fastq_file.readline().rstrip("\n")
        plus_line = fastq_file.readline().rstrip("\n")
        qscores = fastq_file.readline().rstrip("\n")
        yield header, seq, plus_line, qscores


def create_new_header(
    header: str, index_seq_1: str, reverse_complement_index_seq_2: str
) -> str:
    """
    Creates a new header for the demultiplexed FASTQ file implementing the format: @<original_header> <index_seq_1>-<reverse_complement_index_seq_2>.

    Input:
        header (str): The original header from the FASTQ file.
        index_seq_1 (str): The index sequence from the first read.
        reverse_complement_index_seq_2 (str): The reverse complement of the index sequence from the second read.

    Output:
        str: A new header string.
    """
    return f"{header} {index_seq_1}-{reverse_complement_index_seq_2}"


def write_record_to_file(
    output_file_handle: str,
    new_header: str,
    sequence: str,
    plus_line: str,
    qscores: str,
) -> None:
    """
    Writes the demultiplexed records to a new FASTQ file with the new header format: @<original_header> <index_seq_1>-<reverse_complement_index_seq_2>.

    Input:
        output_file_handle (str): The opened file handle to the output FASTQ file.
        new_header (str): The new header for the demultiplexed record.
        sequence (str): The sequence from the FASTQ file.
        plus_line (str): The '+' line from the FASTQ file.
        qscores (str): The quality scores from the FASTQ file.

    Output:
        None
    """
    output_file_handle.write(f"{new_header}\n{sequence}\n{plus_line}\n{qscores}\n")


def reverse_complement(seq: str) -> str:
    """
    Returns the reverse complement of a given DNA sequence.

    Input:
        seq (str): The DNA sequence to be reverse complemented.
    Output:
        str: The reverse complement of the input DNA sequence.
    """
    complement = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}
    for base in seq:
        if base not in complement:
            raise ValueError(
                "Index contains a value that is not a DNA nucleotide (ACGTN)."
            )
    return "".join(complement[base] for base in reversed(seq))
