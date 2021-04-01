# -*- coding: utf-8 -*-
"""
Created on Tue Mar 23 12:22:26 2021

@author: Dafydd
"""
import csv, os


def read_sequence_file(sequence_file):
    '''
    Takes a sequence.txt file location and returns a list of dictionaries
    detailing the names and locations of the fasta files
    Parameters
    ----------
    sequence_file : file location
        the sequence file used to identify the various fasta files

    Returns
    -------
    sequences : list of dict.
        information from sequences.txt file in the form of a list

    '''
    sequences = []
    try:
        with open(sequence_file, newline='') as seq_file:
            reader = csv.DictReader(seq_file, delimiter='	')
            for row in reader:
                sequences.append(row)
    except:
        print("Error: no such file as" + sequence_file)
        return 0
    
    return sequences

def parsnp(sequences):
    '''
    Takes the sequences as outputted by read_sequence_file() and
    uses parSNP to generate a vcf file for use with vg commands
    Parameters
    ----------
    sequences : list of dict.
        DESCRIPTION.

    Returns
    -------
    reference_file : string
        location of the reference file.
    vcf_file : string
        location of the vcf file.

    '''
    reference = sequences[0]['seq_path']
    directory = os.path.abspath(os.path.join(sequences[1]['seq_path'],os.pardir))
    
    return reference, vcf


if __name__ == "__main__":
    sequences = input("Please enter the location of the sequences.txt file:")