#! /usr/bin/env python

#Will have list of known UMI's for comparisons [umi_known]
#regex, numpy, math, argparse are okay to import
#Major dataframe: PCRdup_dict
    #-key = tuple of [UMI, pos, chrom_num, strand]
    #-value = [count]

import argparse
import re
import bioinfo

#NEED TO CREATE A MORE USEFUL MESSAGE FOR HELP ARGPARSE
def get_args():
    parser = argparse.ArgumentParser(description="Deduper code will parse SAM file input and remove duplicates. Adjusting for soft clipping, will output SAM file containing only unique single-paired reads.", add_help=False) #PUT HELP MESSAGE HERE
    parser.add_argument("-f", "--file", help="Absolute path for path to input sorted SAM", required=True, type=str)
    parser.add_argument("-o", "--outfile", help="Absolute Path for output sorted SAM", required=True, type=str)
    parser.add_argument("-u", "--umi", help="Txt file with list of UMIs", required=True, type=str)
    parser.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS, help='Deduper code assumes that SAM file inputs are sorted and only contain mapped single-paired reads.')
    return parser.parse_args()
args = get_args()

#open the sorted SAM file:
with open(args.file) as input_:
    #open output file:
    with open(args.outfile, 'w') as output:

    #global variables
        curr_entry = []
        PCRdup_dict = {}
        umi_known = [] #UMI in UMI_known -> continue
        strand = ""
        header = 0
        unknown_umi = 0
        dup_rm = 0
        chr_dict={}
        uniq_counter = 0

        #read in the UMI file into global variable for future comparisons
        with open(args.umi) as umi_file:
            for i in umi_file:
                i=i.strip('\n')
                umi_known.append(i)

        #parse through each line of the SAM file (for loop)
            #-if the line does not start with "@": (or could check if starts with "NS")
            #    *split the line by tabs to create the columns that we're used to looking through to check flags
            #    *save entry in a string for later (string will be reset at the end of every loop) [entry_str]
        for line_ in input_:
            if(line_.startswith("@")):
                output.write(line_)
                header+=1
            else:                
                #Check col 0 (QNAME) UMI: #ex:NS500451:154:HWKTMBGXX:1:11101:15364:1139:GAACAGGT
            #       *split the QNAME by ":"
            #       *string following the last colon is [UMI] (len(split) - 1 = 7)
            #       *if UMI in UMI_known -> continue
            #        *else break out of this loop and go back to original SAM file                
                line=line_.split('\t') #tab-separated file
                qname = line[0]
                qname = qname.split(":")
                umi = qname[7]

                if(umi not in umi_known):
                    unknown_umi += 1
                elif(umi in umi_known):
                #check if have same alignment position:
                # 1. check POS (col 3)
                #     *save position as a global variable [pos]
                #    *using cigar_clip() function, check CIGAR string (col6) for soft clipping and other letters
                    #       -if CIGAR_str contains S, M, N, D --> pos will need to be adjusted
                    pos = line[3]
                    flag = int(line[1]) #bitwise flag
                    cigar_str = line[5]
                    
                    #2. -check col 2 (RNAME):
                    #    *save value as global variable [chrom_num]
                    chrom_num = line[2]

                    #3. check if on plus or minus strand
                    #parsing the bitwise flag
                    if((flag & 16 == 16)):
                        #checking strand info
                        strand = "-"
                    else:
                        strand = "+"

                    #1.1CHECK CIGAR STRING/soft clipping --> soft_clip() or other outside function
                        #M, D, N, right S
                            # + strand with right s = subtract value from pos
                            # - strand with left s, M, N, D = add value to pos
                        #ignore I,X, and anything else
                        #minus strands should have -1 for correct inclusivity

                        #FUNCTION DEFINITION FOUND IN BIOINFO.PY
                    pos = bioinfo.cigar_clip(cigar_str, strand, pos)

        #-save variables to tuple and compare to dictionary
            #   curr_entry = [UMI, pos, chrom_num, entry_str]
            #      -if(curr_entry in PCRdup_dict):
            #         *increment the counter for PCR duplicates
            #        *do not write anything to output file.
                #   -else if(curr_entry not in PCRdup_dict):
                #     *add to dictionary, with value = [1]
                #    *write entry_str to output file
                # -continue back to the original SAM file
                    curr_entry = (umi, pos, chrom_num, strand)
                    if(curr_entry in PCRdup_dict):
                        PCRdup_dict[curr_entry] += 1
                        dup_rm += 1
                    elif(curr_entry not in PCRdup_dict):
                        PCRdup_dict[curr_entry] = 1
                        uniq_counter+=1
                        output.write(line_)

                        #keep second dictionary to identify the unique reads in each chromsome
                        if(chrom_num in chr_dict):
                            chr_dict[chrom_num] +=1
                        elif(chrom_num not in chr_dict):
                            chr_dict[chrom_num] = 1
                                     

        #remember to CLEAR THE DICTIONARY AFTERWARDS (to save space)
        PCRdup_dict={}
    
        print("Header: ", header)
        print("Unknown Umi's: ", unknown_umi)
        print("Duplicates Removed: ", dup_rm)
        print("Total overall unique reads: ", uniq_counter)
        # for key,value in chr_dict.items():
        #     print(key, value, sep = ',')