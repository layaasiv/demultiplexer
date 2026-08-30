![CI](https://github.com/layaasiv/demultiplexer/actions/workflows/ci.yml/badge.svg)

# dmux

A command-line tool for demultiplexing paired-end, dual-indexed FASTQ data. The tool handles three possible mutually exclusive cases:

  1. Records with unknown indexes: Either index 1 or index 2 contain "N" nucleotides or are not included in the list of possible indexes specified by the user.
  2. Records with hopped indexes: Index 1 and the reverse complement of index 2 are not the same.
  3. Records with matched indexes: Index 1 and the reverse complement of index 2 are the same. This is the ideal situation because there is no ambiguity about which index these records belong to. 

## Features
- Streaming FASTQ processing
- Dual-index matching
- Reverse-complement handling
- Unknown-index filtering
- Index-hopping detection
- Compressed FASTQ input/output
- Command-line interface
- Automated tests
- CI with GitHub Actions

## Installation

pip install .

## Usage

dmux \
  -r1 ... \
  -i1 ... \
  -i2 ... \
  -r2 ... \
  -i indexes.txt \
  -o output/

## Input

Four Illumina-sequencing-style FASTQ files (compressed) in which:
  - read 1: Forward read sequence
  - read 2: Index 1 sequence
  - read 3: Index 2 sequence
  - read 4: Reverse read sequence

A text file containing indexes used in the experiment, representing the index sequences that can be expected in the sequencing data.

Path to directory where output files will be saved.

## Output

Compressed FASTQ files for hopped, unknown and each of the indexes in the given ```indexes.txt``` file. Each file contains the corresponding reads from the multiplexed file. Each category (hopped, unknown, and all the indexes) will have two output files, one each for the forward read and reverse read. 

Additionally, three TSV files:
  - ```match_percents.tsv```: Counts and proportions of reads identified as belonging to each index.
  - ```hop_percents.tsv```: Counts and proportions of indexes belonging to a hopped pair.
  - ```total_percents.tsv```: Count and proportion of reads in each category (hopped, unknown and matched).

Also output to stdout, which shows the number of records categorized as each index, hopped, and unknown.

## Testing

Thorough tests of the input data, output, and CLI are included in `tests/`. Tests are run using a synthetic test dataset stored in `tests/data/`. Run tests with this command:

```
pytest
```

## Development

ruff check .
ruff format --check .
pytest