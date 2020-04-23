#!/usr/bin/python3

"""
Author: Saul Pierotti
Mail: saulpierotti.bioinfo@gmail.com
Last updated: 21/04/2020

Takes in input a fasta file and a list of IDs
It returns a fasta file containing only the sequences with header
not containing a string from in ID
"""

import sys


def parse_fasta(fastapath):
    headers = []
    seqs = []
    with open(fastapath) as fastafile:
        for line in fastafile:
            line = line.rstrip()
            if line[0] == ">":
                headers.append(line)
                seqs.append("")
            else:
                seqs[-1] += line
    return headers, seqs


def parse_ids(idpath):
    idlist = []
    with open(idpath) as idfile:
        for line in idfile:
            idlist.append(line.rstrip())
    return idlist


def strip_seqs(fastalist, idlist):
    indexes_to_remove = []
    for idname in idlist:
        for i in range(len(fastalist[0])):
            if idname in fastalist[0][i]:
                indexes_to_remove.append(i)
    for i in range(len(fastalist[0])):
        if i not in indexes_to_remove:
            print(fastalist[0][i] + "\n" + fastalist[1][i])


if __name__ == "__main__":
    fastalist = parse_fasta(sys.argv[1])
    idlist = parse_ids(sys.argv[2])
    strip_seqs(fastalist, idlist)
