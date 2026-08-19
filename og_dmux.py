#!/usr/bin/env python

# import modules
import argparse
import gzip
import bioinfo

# set up argparse 
def get_args():
    parser = argparse.ArgumentParser(description="A program to k-merize our data set at different sizes for k and plot the result")
    parser.add_argument("-r1", "--r1", help="Name of r1 file", required=True)
    parser.add_argument("-r2", "--r2", help="Name of r2 file", required=True)
    parser.add_argument("-r3", "--r3", help="Name of r3 file", required=True)
    parser.add_argument("-r4", "--r4", help="Name of r4 file", required=True)
    parser.add_argument("-q", "--q_cutoff", help="Quality score cutoff for indexes", required=False)
    parser.add_argument("-i", "--indexestxt", help="Name of text tile with all possible indexes", required=True)
    return parser.parse_args()

# define arg variables
args = get_args()
fr1 = args.r1
fr2 = args.r2
fr3 = args.r3
fr4 = args.r4
q_cutoff = args.q_cutoff
indexes_txt = args.indexestxt

# define functions
def rev_comp(seq: str) -> str:
    '''Takes a DNA sequence and returns its reverse complement as a string'''
    rc_seq = ''
    comp_bases = {'A':'T', 'T':'A', 'G':'C', 'C':'G', 'N':'N'}
    rev = seq[::-1]
    for  base in rev:
        rc_seq += comp_bases[base]
    return rc_seq

# define variables
r1_rec = ['','','','']
r2_rec = ['','','','']
r3_rec = ['','','','']
r4_rec = ['','','','']

indexes = [] 

with open(indexes_txt, 'r') as fh:
    num_lines = 0
    for line in fh:
        num_lines += 1
        if num_lines > 1:
            line = line.strip('\n').split('\t')
            ind = line[4]
            indexes.append(line[4])

matched_counter = {}
hopped_counter = {}
unknown_counter = 0

# open all output files 
output_files = {}
for index in indexes: 
    fhr1 = open(index + '_R1.fq', 'w')
    fhr2 = open(index + '_R2.fq', 'w')
    output_files[index] = [fhr1, fhr2]

unk_r1 = open('unknown_R1.fq', 'w')
unk_r2 = open('unknown_R2.fq', 'w')
hopped_r1 = open('hopped_R1.fq', 'w')
hopped_r2 = open('hopped_R2.fq', 'w')

mat_per = open('mat_percents.tsv', 'w')
hop_per = open('hop_percents.tsv', 'w')
tot_per = open('total_percents.tsv', 'w')

# read all 4 files simultaneously, record by record, and save 1 record at a time from each file into its corresponding list
with gzip.open(fr1, 'rt') as r1, gzip.open(fr2, 'rt') as r2, gzip.open(fr3, 'rt') as r3, gzip.open(fr4, 'rt') as r4:
    while True:
        r1_rec[0] = r1.readline().strip('\n')
        if r1_rec[0] == '':
            break
        r1_rec[1] = r1.readline().strip('\n')
        r1_rec[2] = r1.readline().strip('\n')
        r1_rec[3] = r1.readline().strip('\n')

        r2_rec[0] = r2.readline().strip('\n')
        r2_rec[1] = r2.readline().strip('\n')
        r2_rec[2] = r2.readline().strip('\n')
        r2_rec[3] = r2.readline().strip('\n')

        r3_rec[0] = r3.readline().strip('\n')
        r3_rec[1] = r3.readline().strip('\n')
        r3_rec[2] = r3.readline().strip('\n')
        r3_rec[3] = r3.readline().strip('\n')

        r4_rec[0] = r4.readline().strip('\n')
        r4_rec[1] = r4.readline().strip('\n')
        r4_rec[2] = r4.readline().strip('\n')
        r4_rec[3] = r4.readline().strip('\n')

        # reassign some variables for ease of access moving forward
        header = r1_rec[0] 
        index1 = r2_rec[1]
        i1_qscore = r2_rec[3]
        index2 = r3_rec[1]
        i2_qscore = r3_rec[3]

        rc_index2 = rev_comp(index2)
        ind_pair = (index1, rc_index2)
        # create the new headers with the indexes attached to the end and save to variable for ease of access
        header_r1 = r1_rec[0] + ' ' + index1 + '-' + rc_index2
        header_r4 = r4_rec[0] + ' ' + index1 + '-' + rc_index2

        if 'N' in index1 or 'N' in rc_index2:
            unk_r1.write(header_r1 +'\n' + r1_rec[1] +'\n'+ r1_rec[2] +'\n'+ r1_rec[3] +'\n')
            unk_r2.write(header_r4 + '\n' + r4_rec[1] + '\n' + r4_rec[2] + '\n' + r4_rec[3] + '\n')
            unknown_counter += 1
        
        # include q_cutoff condition 

        elif index1 not in indexes or rc_index2 not in indexes: 
            unk_r1.write(header_r1 +'\n' + r1_rec[1] +'\n'+ r1_rec[2] +'\n'+ r1_rec[3] +'\n')
            unk_r2.write(header_r4 + '\n' + r4_rec[1] + '\n' + r4_rec[2] + '\n' + r4_rec[3] + '\n')
            unknown_counter += 1

        elif index1 == rc_index2:
            output_files[index1][0].write(header_r1 +'\n' + r1_rec[1] +'\n'+ r1_rec[2] +'\n'+ r1_rec[3] +'\n')
            output_files[index1][1].write(header_r4 + '\n' + r4_rec[1] + '\n' + r4_rec[2] + '\n' + r4_rec[3] + '\n')

            if index1 in matched_counter:
                matched_counter[index1] += 1
            else:
                matched_counter[index1] = 1
        
        elif index1 != rc_index2: 
            hopped_r1.write(header_r1 +'\n' + r1_rec[1] +'\n'+ r1_rec[2] +'\n'+ r1_rec[3] +'\n')
            hopped_r2.write(header_r4 + '\n' + r4_rec[1] + '\n' + r4_rec[2] + '\n' + r4_rec[3] + '\n')

            if ind_pair in hopped_counter:
                hopped_counter[ind_pair] += 1
            else:
                hopped_counter[ind_pair] = 1

        else:
            raise Exception('Impossible')
        
        # reset variables 
        r1_rec = ['','','','']
        r2_rec = ['','','','']
        r3_rec = ['','','','']
        r4_rec = ['','','','']

for ind in matched_counter: 
    print(ind, '\t', matched_counter[ind])

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
