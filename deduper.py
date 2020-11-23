#!/usr/bin/env python

import re
import argparse
import subprocess

def get_args():
    '''adds argparse arguments to allow script to run from command line'''
    parser = argparse.ArgumentParser(description="Susan's Deduper: Removes PCR duplicates and returns three output files: deduped.sam, removed_duplicates.sam, unmatched_umi.sam. Sorts SAM file and returns first occurence of a duplicate by position. Requires list of known UMIs, and samtools to be installed.")
    parser.add_argument("-f", "--file", help="input sam file", type=str, required=True)
    parser.add_argument("-u", "--umis", help="file containing list of UMIs", type=str, required=True)
    parser.add_argument("-p", "--paired_end", help="set for paired-end data", action="store_true")
    return parser.parse_args()

args = get_args()

#store argparse arguments
input_file = args.file
umi_file = args.umis
paired_end = args.paired_end

#print error message if paired-end data is selected
if paired_end == True:
    raise Exception("Paired-end data not currently supported. Coming soon in version 2.0!")

#use subprocess to sort samfile
subprocess.call("samtools sort {} -o temp.sorted.sam".format(input_file), shell=True)

#save output to sorted_sam
sorted_sam = "./temp.sorted.sam"

def create_umi_dict(umi_file):
    '''Populates a dictionary with UMIs as keys and 0 as values, using a file containing UMI sequences'''
    umi_dict = {}
    with open(umi_file) as uf:
        for line in uf:
            line = line.strip()
            umi_dict.setdefault(line, 0)
    uf.close()
    return umi_dict

#save umi information in dictionary
umi_dict = create_umi_dict(umi_file)

def reverse_strand(cigar):
    '''Takes cigar string of reverse strand and returns number to add to position for extract_read_info'''
    #find matches, soft clipping on the right, deletions and gaps to add to position
    position_adjust = 0
    match = re.findall(r"\d+S$|\d+M|\d+D|\d+N", cigar)
    if match:
        for i in match:
            position_adjust += int(i[:-1])
    print(position_adjust)
    return position_adjust

def extract_read_info(line):
    '''Takes a SAM file read and extracts UMI, RNAME (chromosome), POS adjusted for soft clipping, and strand. Puts info into a list'''
    soft_clip = 0
    split_line = line.split('\t')
    read_umi = (split_line[0].split(':'))[-1]
    position = int(split_line[3])
    if ((int(split_line[1]) & 16)) == 16:
        strand = 'reverse'
        adj_pos = position + reverse_strand(split_line[5])
        print(adj_pos)
    else:
        strand = 'forward'
        match = re.match(r"^(\d+)S", split_line[5])
        if match:
            soft_clip = int(match.group(1))
        adj_pos = position - soft_clip

    read_info = (read_umi, split_line[2], adj_pos, strand)

    return read_info


#open all files
open_sam = open(sorted_sam)
dedup_out = open('deduped.sam', 'w')
duplicates_out = open('removed_duplicates.sam', 'w')
no_match_out = open('unmatched_umi.sam', 'w')


def deduper(sam_file):
    '''Takes sorted sam file and outputs: a file with PCR duplicates removed (output_file), a file of removed duplicates (dup_file), and a file of unmatched UMIs (unmatched_file)'''
    #could add a global duplicate counter, or counter per chromosome?
    duplicate_check = {}
    current_chromosome = 1
    read_counter = 0
    duplicate_counter = 0
    unmatched_counter = 0
    output_counter = 0
    for line in sam_file:
        line = line.strip()
        #check if line is a header line
        if line.startswith('@'):
            dedup_out.write(line + '\n')
        #check if line is empty to avoid errors at end of file
        elif line != '':
            read_counter += 1
            read_info = extract_read_info(line)
            #check if chromosome matches current chromosome
            if read_info[1] != current_chromosome:
                current_chromosome = read_info[1]
                duplicate_check = {}
            #check if UMI is known
            if read_info[0] not in umi_dict:
                no_match_out.write(line + '\n')
                unmatched_counter += 1
            else:
                #check if read is a duplicate
                if read_info not in duplicate_check:
                    dedup_out.write(line + '\n')
                    duplicate_check.setdefault(read_info, 0)
                    output_counter += 1
                else:
                    duplicates_out.write(line + '\n')
                    duplicate_counter += 1
    return read_counter, duplicate_counter, unmatched_counter, output_counter

read_counter, duplicate_counter, unmatched_counter, output_counter = deduper(open_sam)

#close all opened files
open_sam.close()
dedup_out.close()
duplicates_out.close()
no_match_out.close()

#remove sorted sam file here
subprocess.call("rm temp.sorted.sam", shell=True)

print("Number of input reads: " + str(read_counter))
print("Percent PCR duplicates: " + str(duplicate_counter/read_counter*100) + "%")
print("Unmatched UMIs: " + str(unmatched_counter))
print("Output deduplicated reads: " + str(output_counter))

