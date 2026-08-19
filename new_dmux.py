#!/usr/bin/env python

import dmux_tools as dt
import argparse
import gzip

__version__ = "1.0.0"

def get_args():
    parser = argparse.ArgumentParser(description="""
    A program to demultiplex a set of paired-end, dual-indexed FASTQ files. Reads are classified as matched, index-hopped, or unknown.
    The program takes as input the four FASTQ files (R1, R2, I1, I2) and a text file containing the list of known indexes. The output consists of separate FASTQ files for each matched index, as well as files for unknown and hopped reads.
    """)
    parser.add_argument("--version", action="version", version=f"dmux {__version__}")
    parser.add_argument("-i", "--indexes", help="Text file containing list of known indexes", required=True)
    parser.add_argument("-r1", "--read1", help="Zipped read 1 FASTQ file (forward read)", required=True)
    parser.add_argument("-i1", "--index1", help="Zipped index 1 FASTQ file", required=True)
    parser.add_argument("-i2", "--index2", help="Zipped index 2 FASTQ file", required=True)
    parser.add_argument("-r2", "--read2", help="Zipped read 2 FASTQ file (reverse read)", required=True)
    parser.add_argument("-o", "--outputpath", help="Path to output directory", required=False, default="./")
    return parser.parse_args()

def main():
    args = get_args()
    index_file = args.indexes
    r1_file = args.read1
    i1_file = args.index1
    i2_file = args.index2
    r2_file = args.read2
    output_path = args.outputpath
    indexes = []

    with open(index_file, "r") as fh:
        for line in fh:
            line = line.strip('\n').split('\t')
            indexes.append(line[1])

    # initialize counters for matched and hopped reads
    matched_counter = {k:0 for k in indexes}
    hopped_counter = {}
    unknown_counter = 0

    # create the output files
    # output fastq files for correctly matched indexes
    output_files = {}
    for ind in indexes:
        output_files[ind] = [gzip.open(f"{output_path}/{ind}_R1.fastq.gz", "wt"),
                            gzip.open(f"{output_path}/{ind}_R2.fastq.gz", "wt")]
    # output fastq files for unknown indexes (those not specified in the indexes file, or those with 'N' bases)
    unk_r1 = gzip.open(f"{output_path}/unknown_R1.fastq.gz", "wt")
    unk_r2 = gzip.open(f"{output_path}/unknown_R2.fastq.gz", "wt")
    # output fastq files for hopped indexes (those with index sequences that do not match)
    hopped_r1 = gzip.open(f"{output_path}/hopped_R1.fastq.gz", "wt")
    hopped_r2 = gzip.open(f"{output_path}/hopped_R2.fastq.gz", "wt")
    # output tsv files for matched, hopped, and total counts
    mat_per = open(f"{output_path}/mat_percents.tsv", 'w')
    hop_per = open(f"{output_path}/hop_percents.tsv", 'w')
    tot_per = open(f"{output_path}/total_percents.tsv", 'w')


    with gzip.open(r1_file, "rt") as r1, gzip.open(r2_file, "rt") as r2, gzip.open(i1_file, "rt") as i1, gzip.open(i2_file, "rt") as i2:
        for r1_rec in dt.fastq_parser(r1):
            i1_rec = next(dt.fastq_parser(i1))
            i2_rec = next(dt.fastq_parser(i2))
            r2_rec = next(dt.fastq_parser(r2))

            index_seq_1 = i1_rec[1]
            index_seq_2 = i2_rec[1]
            rc_index_seq_2 = dt.reverse_complement(index_seq_2)

            header_r1 = dt.create_new_header(r1_rec[0], index_seq_1, rc_index_seq_2)
            header_r2 = dt.create_new_header(r2_rec[0], index_seq_1, rc_index_seq_2)

            if 'N' in index_seq_1 or 'N' in index_seq_2 or index_seq_1 not in indexes or rc_index_seq_2 not in indexes:
                dt.write_record_to_file(unk_r1, header_r1, r1_rec[1], r1_rec[2], r1_rec[3])
                dt.write_record_to_file(unk_r2, header_r2, r2_rec[1], r2_rec[2], r2_rec[3])
                unknown_counter += 1

            elif index_seq_1 == rc_index_seq_2:
                dt.write_record_to_file(output_files[index_seq_1][0], header_r1, r1_rec[1], r1_rec[2], r1_rec[3])
                dt.write_record_to_file(output_files[index_seq_1][1], header_r2, r2_rec[1], r2_rec[2], r2_rec[3])
                matched_counter[index_seq_1] += 1
            
            elif index_seq_1 != rc_index_seq_2:
                dt.write_record_to_file(hopped_r1, header_r1, r1_rec[1], r1_rec[2], r1_rec[3])
                dt.write_record_to_file(hopped_r2, header_r2, r2_rec[1], r2_rec[2], r2_rec[3])
                counter_key = sorted([index_seq_1, rc_index_seq_2])
                counter_key = "-".join(counter_key)
                if counter_key not in hopped_counter:
                    hopped_counter[counter_key] = 1
                else:
                    hopped_counter[counter_key] += 1

            else:
                raise Exception("Unexpected case encountered during demultiplexing.")

    print("Matched reads per index:")
    for ind in matched_counter: 
        print(ind, '\t', matched_counter[ind])

    print("Hopped reads per index:")
    for ind in hopped_counter:
        print(ind, '\t', hopped_counter[ind])

    print(f'Unknown reads = {unknown_counter}')

    # assigning more variables (for the output data)
    total_hop = sum(hopped_counter.values())
    total_matched = sum(matched_counter.values())
    grand_total = unknown_counter + total_matched + total_hop


    # calculations for matched reads and writing into tsv file
    print(f'Index name \t Count \t % Sample in Matched \t % Sample in Total', file=mat_per)
    for item in matched_counter: 
        per_ind_sam = (matched_counter[item]/total_matched) * 100
        per_ind_tot = (matched_counter[item]/grand_total) * 100 
        print(f'{item} \t {matched_counter[item]} \t {per_ind_sam} \t {per_ind_tot}', file=mat_per)

    # calculations for hopped reads and writing into tsv file 
    print(f'Index pair \t Count \t % Index pair in Total', file=hop_per)
    for item in hopped_counter:
        per_hop_tot = (hopped_counter[item]/grand_total) * 100 
        print(f'{item} \t {hopped_counter[item]} \t {per_hop_tot}', file=hop_per)

    # creating totals tsv file
    print(f'Index pair \t Count \t % Index pair in Total', file=tot_per)
    print(f'Matched \t {total_matched} \t {(total_matched/grand_total)*100}', file=tot_per)
    print(f'Hopped \t {total_hop} \t {(total_hop/grand_total)*100}', file=tot_per)
    print(f'Unknown \t {unknown_counter} \t {(unknown_counter/grand_total)*100}', file=tot_per)

    # close all files that were opened 
    for index in output_files:
        output_files[index][0].close()
        output_files[index][1].close()

    unk_r1.close()
    unk_r2.close()
    hopped_r1.close()
    hopped_r2.close()
    mat_per.close()
    hop_per.close()
    tot_per.close()

if __name__ == "__main__":
    main()