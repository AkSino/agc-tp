#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""OTU clustering"""

import argparse
import sys
import os
import gzip
import statistics
from collections import Counter
# https://github.com/briney/nwalign3
# ftp://ftp.ncbi.nih.gov/blast/matrices/
# import nwalign3 as nw

__author__ = "Aurélien"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Aurélien"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Aurélien"
__email__ = "burieaurel@eisti.eu"
__status__ = "Developpement"


def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', '-amplicon_file', dest='amplicon_file', type=isfile, required=True,
                        help="Amplicon is a compressed fasta file (.fasta.gz)")
    parser.add_argument('-s', '-minseqlen', dest='minseqlen', type=int, default = 400,
                        help="Minimum sequence length for dereplication")
    parser.add_argument('-m', '-mincount', dest='mincount', type=int, default = 10,
                        help="Minimum count for dereplication")
    parser.add_argument('-c', '-chunk_size', dest='chunk_size', type=int, default = 100,
                        help="Chunk size for dereplication")
    parser.add_argument('-k', '-kmer_size', dest='kmer_size', type=int, default = 8,
                        help="kmer size for dereplication")
    parser.add_argument('-o', '-output_file', dest='output_file', type=str,
                        default="OTU.fasta", help="Output file")
    return parser.parse_args()


def read_fasta(amplicon_file, minseqlen):
    with gzip.open(amplicon_file) as file:
        sequences = file.readlines()
        seqs = ""
        for sequence in sequences:
            #print("sequence")
            seq = sequence.replace(b"\n", b"")
            seq = seq.decode('utf8')
            #print(seq)
            for character in seq:
                if character not in "TGAC":
                    if len(seqs)>=minseqlen:
                        yield seqs
                        #print(seqs)
                    seq = ""
                    seqs = ""
                    break
            seqs += seq
        #print(seqs)
        yield seqs

def dereplication_fulllength(amplicon_file, minseqlen, mincount):
    dict_derep = {}
    seq = read_fasta(amplicon_file, minseqlen)
    for sequence in seq:
        counter = 0
        seq2 = read_fasta(amplicon_file, minseqlen)
        for sequence2 in seq2:
            if sequence2 == sequence:
                counter += 1
        if counter < mincount:
            continue
        dict_derep[sequence] = counter
    dict_derep = sorted(dict_derep.items(), key=lambda t: t[1], reverse=True)
    for derep in dict_derep:
        yield [derep[0], derep[1]]





def get_chunks(sequence, chunk_size):
    seq_length = len(sequence)
    seq_list = []
    treshold = int(seq_length) // int(chunk_size)
    if treshold <4:
        raise ValueError("Change chunk size")
    for i in range(treshold):
        seq = sequence[i*chunk_size:(i+1)*chunk_size]
        seq_list.append(seq)
    return seq_list

def get_unique(ids):
    return {}.fromkeys(ids).keys()


def common(lst1, lst2):
    return list(set(lst1) & set(lst2))

def cut_kmer(sequence, kmer_size):
    pass

def get_identity(alignment_list):
    pass

def chimera_removal(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
    pass

def abundance_greedy_clustering(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
    pass

def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))

def write_OTU(OTU_list, output_file):
    pass
#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()
    amplicon_file = args.amplicon_file
    minseqlen = args.minseqlen
    mincount = args.mincount
    derep = dereplication_fulllength(amplicon_file, minseqlen, mincount)
    print(next(derep))
    seq = "TGGGGAATATTGCACAATGGGCGCAAGCCTGATGCAGCCATGCCGCGTGTATGAAGAAGGCCTTCGGGTTGTAAAGTACTTTCAGCGGGGAGGAAGGTGTTGTGGTTAATAACCGCAGCAATTGACGTTACCCGCAGAAGAAGCACCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGGAAAGCGCA"
    chunks = get_chunks(seq,args.chunk_size)


if __name__ == '__main__':
    main()
