#!/usr/bin/python3

"""
Author: Saul Pierotti
Mail: saulpierotti.bioinfo@gmail.com
Last updated: 21/04/2020

Takes in input a fasta file and returns a random permutation of the contained
sequences
"""

import sys
import random


def parse_fasta(fastapath):
    fastalist = []
    with open(fastapath) as fastafile:
        for line in fastafile:
            line = line.rstrip()
            if line[0] == ">":
                fastalist.append([line])
                fastalist[-1].append("")
            else:
                fastalist[-1][1] += line
    return fastalist


def print_fasta_scrambled(fastalist):
    random.shuffle(fastalist)
    for seq in fastalist:
        print(seq[0] + "\n" + seq[1])


if __name__ == "__main__":
    fastalist = parse_fasta(sys.argv[1])
    print_fasta_scrambled(fastalist)
