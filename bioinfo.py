#!/usr/bin/env python

import argparse
import re

'''This module is a collection of useful bioinformatics functions
written during the Bioinformatics and Genomics Program coursework.
You should update this docstring to reflect what you would like it to say'''

__version__ = "0.3"         # Read way more about versioning here:
                            # https://en.wikipedia.org/wiki/Software_versioning

DNA_bases = set('ATGCNatcgn')
RNA_bases = set('AUGCNaucgn')

def convert_phred(letter: str) -> int:
    """Converts a single character into a phred score"""
    return ord(letter) - 33

def qual_score(phred_score: str) -> float:
    """Take in quality score string. Return average quality score of the entire input string."""
    sum = 0
    for i in phred_score:
        sum+=convert_phred(i)
    return sum/len(phred_score)

def validate_base_seq(seq,RNAflag=False):
    '''This function takes a string. Returns True if string is composed
    of only As, Ts (or Us if RNAflag), Gs, Cs. False otherwise. Case insensitive.'''
    return set(seq)<=(RNA_bases if RNAflag else DNA_bases)

def gc_content(DNA):
    '''Returns GC content of a DNA sequence as a decimal between 0 and 1.'''
    #pass
    assert validate_base_seq(DNA), "String contains invalid characters"  #this is the error message that will be shown
    seq = DNA.upper()                        #ensure string is all uppercase
    GC = seq.count('G') + seq.count('C')     #count number of Gs + Cs
    return GC/len(seq)                       #calculate GC content of entire sequence

def oneline_fasta(input_file:str, output_file:str) -> str:
    '''Take in string parameters of the input file name where sequences have '\n' in the middle, and the output file name. 
        Function removes '\n' from the middle of sequences in the input file. Output file will have headers on one line, 
        entire sequence on one line. Will return the string 'Successful' when the function is complete'''
    with open(input_file) as input:
        with open(output_file, 'w') as output:
            counter = 0 
            for line in input:
                if ">" in line and counter!=0:
                    output.write('\n')
                if ">" not in line:
                    line=line.strip('\n')

                output.write(line)
                counter+=1
    return('Successful')         


def cigar_clip(cigar_: str, strand_: str, pos_: int) -> int:
    #regex for the cigar string
    #regular expression to split up the entire string based on separating numbers from any letters

#regex possible patterns: 
    #r'(\d+)(\w)', '40M25N5M'
    #[0-9]+[MIDNSHPX](,[0-9]+[MIDNSHPX])+
    #r'([0-9]+)([MIDNSHPX=])'
    #r'([0-9]+)([MIDNSX=])' vs r'([0-9]+[MIDNSX=])' --> what would be the difference in output matches?
    #re.findall(r'([0-9]+)([MIDNSX=])', cigar)
    #ignore I, X, = ?

    #1.1CHECK CIGAR STRING/soft clipping --> soft_clip() or other outside function
                        #M, D, N, right S
                            # + strand with right s = subtract value from pos
                            # - strand with left s, M, N, D = add value to pos
                        #minus strands should have -1 for correct inclusivity
    
    pos_ = int(pos_)
    #ignore I, X, = ?
    cigar_values = re.findall(r'([0-9]+)([MDNS])', cigar_)
    #print(cigar_values)
    #will produce list of pairs of tuples (num, letter)

#is this using the right variables???
    counter = 0 #used to keep track of the right vs left S (soft clipping)
    if(strand_ == "-"):
        pos_ = int(pos_) - 1
    for i in cigar_values:
        number = int(i[0])
        letter = i[1]

        if(strand_ == "+"):
            #check for an S on the left side
            #everything other letter won't affect the position
            if(letter == 'S' and counter!=0): #soft clipping was found
                pos_ = int(pos_) - number


        elif(strand_ == "-"):
            #left_S = re.findall(r'([0-9]+)([S]')
            if(letter == 'S' and counter ==0):
                pos_ = int(pos_) + number
                #don't adjust the value?
                #pos_ = int(pos_)
            elif(letter == 'M'):
                pos_ = int(pos_) + number
            elif(letter == 'N'):
                pos_ = int(pos_) + number
            elif(letter == 'D'):
                pos_ = int(pos_) + number

        #how the math works will depend on the strand it's on
        counter+=1
    return pos_

   #     def soft_clip(str: CIGAR_str, int: pos_5) -> str
    #        return pos_5 
     #   soft_clip("2S12M", 113) -> 111
       #     input: CIGAR_str = "2S12M", pos_5 = 113 --> return: pos_5 = 111

if __name__ == "__main__":
    write tests for functions above
    assert convert_phred("I") == 40, "wrong phred score for 'I'"
    assert convert_phred("C") == 34, "wrong phred score for 'C'"
    assert convert_phred("2") == 17, "wrong phred score for '2'"
    assert convert_phred("@") == 31, "wrong phred score for '@'"
    assert convert_phred("$") == 3, "wrong phred score for '$'"
    print("Convert_phred function is working!")
    phred_score: str = "FFHHHHHJJJJIJIJJJIJJJJJJIIIJJJEHJJJJJJJIJIDGEHIJJFIGGGHFGHGFFF@EEDE@C??DDDDDDD@CDDDDBBDDDBDBDD@"
    assert qual_score(phred_score) == 37.62105263157895, "wrong average phred score"
    print("You calculated the correct average phred score")
    assert validate_base_seq("AATAGAT") == True, "DNA string not recognized"
    assert validate_base_seq("AGCTACTGNNNNCTACTG") == True, "DNA string not recognized"
    print("Correctly identified DNA")
    assert validate_base_seq("Coding is fun") == False, "Non-DNA identified as DNA"
    assert validate_base_seq("This week was exhausting") == False, "Non-DNA identified as DNA"
    print("Correctly determined non-DNA")
    assert gc_content("GCGCGC") == 1
    assert gc_content("AATTATA") == 0
    assert gc_content("GCATGCAT") == 0.5
    print("correctly calculated GC content")
    #I will include this bioinfo_test.txt file with my submission
    #testing = oneline_fasta('bioinfo_test.txt', 'bioinfo_result.txt')
    #assert(testing == 'Successful')
    #print('Converted FASTA correctly')
    soft = cigar_clip('3S22M2S', '-', 33)
    #assert(soft == 31)
    #print(soft)
        #my bioinfo_test.txt file:
        # >seq1
        # A
        # C
        # T
        # G
        # >seq2
        # AAAAAAA
        # >seq3
        # ATCTCA
        # GTACAG
        # TGACST
        # >seq4
        # shfgosk;l