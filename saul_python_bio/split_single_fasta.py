#!python

"""
This script takes in input a fasta file and produces a set of fasta files, one for
each sequence that it contained.
The first argument is the original fasta file, while the second argument is the
folder where to place the resulting files.
"""

import sys
import os


def main(fa_file, fa_dir):
    """
    Splits fa_file in single-sequence fasta files and puts them in fa_dir
    """
    os.system("mkdir " + fa_dir)
    seq_num = 0
    out_content = ""
    fa_basename = fa_file.split("/")[-1].split(".")[0]
    extension = "." + ".".join(fa_file.split("/")[-1].split(".")[1:])
    first_seq = True
    with open(fa_file) as in_handle:
        for line in in_handle:
            if line != "":
                if line[0] == ">":
                    if not first_seq:
                        with open(out_file, "w") as out_handle:
                            out_handle.write(out_content)
                            out_content = ""  # reset the string to write
                    else:
                        first_seq = False
                    seq_num += 1
                    out_file = fa_dir + fa_basename + "_seq" + str(seq_num) + extension
                    out_content += line
                else:
                    out_content += line
    if out_content != "":  # last sequence
        with open(out_file, "w") as out_handle:
            out_handle.write(out_content)
            out_content = ""  # reset the string to write


if __name__ == "__main__":
    assert len(sys.argv) == 3
    FA_FILE = sys.argv[1]
    FA_FOLDER = sys.argv[2]
    main(FA_FILE, FA_FOLDER)
