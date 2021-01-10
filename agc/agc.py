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
__copyright__ = "EISTI"
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
    """Creates a sequence generator.
        :Parameters:
            amplicon_file: fasta.gz file
            minseqlen: minimum length of sequences
        Returns: A generator of sequences
    """
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
    """Creates a sequence, count generator.
        :Parameters:
            amplicon_file: fasta.gz file
            minseqlen: minimum length of sequences
            mincount: minimum count to keep sequence
        Returns: A generator of unique sequences having an occurrence greater
                 than mincount as well as their occurrence. The sequences
                 will be returned in descending order of occurrence
    """
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
    """Transforms a sequence into chunks.
        :Parameters:
            sequence: sequence
            chunk_size: size of each chunk
        Returns: List of chunks
    """
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
    """Get unique ids
    """
    return {}.fromkeys(ids).keys()


def common(lst1, lst2):
    """Get commun objects
    """
    return list(set(lst1) & set(lst2))

def cut_kmer(sequence, kmer_size):
    """Yields the kmer of sequence
        :Parameters:
            sequence : sequence
            kmer_size : size of kmer
        Returns: yields the kmers
    """
    for i in range(len(sequence)-kmer_size+1):
        yield sequence[i:i+kmer_size]

def get_identity(alignment_list):
    """Get identity between 2 alignments
        :Parameters:
            alignment_list: list of 2 alignments
        Returns: % of identity
    """
    min_length = min(len(alignment_list[0]),len(alignment_list[1]))
    count = 0
    for i in range(min_length):
        if alignment_list[0][i] == alignment_list[1][i]:
            count += 1
    percent = count/min_length * 100
    return percent

def chimera_removal(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
    """
        :Parameters:
            amplicon_file: fasta.gz file
            minseqlen: minimum length of sequences
            mincount: minimum count to keep sequence
            chunk_size: size of each chunk
            kmer_size : size of kmer
        Returns: A non-chimeric sequence generator in the format: yield [sequence, count]
    """
    for s, c in dereplication_fulllength(amplicon_file, minseqlen, mincount):
        chunk_dict = {}
        try:
            chunks = get_chunks(s, chunk_size)[:4]
            for chunk in chunks:
                dict_sequence = {}
                sequence_without_chunk = s.replace(chunk,"")
                kmers = cut_kmer(sequence_without_chunk,chunk_size)
                for kmer in kmers:
                    identity = get_identity([kmer,chunk])
                    if len(dict_sequence) < 8 and identity > 0:
                        dict_sequence[kmer] = identity
                    seq_with_min_identity = min(dict_sequence, default=10000000, key=lambda k: dict_sequence[k])
                    if identity > dict_sequence[seq_with_min_identity]:
                        del dict_sequence[seq_with_min_identity]
                        dict_sequence[kmer] = identity
                sorted_dict_sequence = sorted(dict_sequence, key=lambda k: dict_sequence[k])
                chunk_dict[chunk] = {kmer:dict_sequence[kmer] for kmer in sorted_dict_sequence}
            list_sequences = list(chunk_dict.keys())
            seq_parente1, seq_parente2 = "", ""
            for chunk in chunk_dict[list_sequences[0]].keys():
                count = 1
                for i in range(1, len(list_sequences)):
                    for chunkk in chunk_dict[list_sequences[i]].keys():
                        if chunk == chunkk:
                            count += 1
                        if count == len(list_sequences):
                            if seq_parente1 != "" and chunk == chunkk:
                                seq_parente1 = chunk
                            elif seq_parente2 != "" and chunk == chunkk:
                                seq_parente2 = chunkk
                        break
            if seq_parente1 == "":
                max, maxx, chunkk, chunkkk = 0, 0, "", ""
                for i in range(len(list_sequences)):
                    for chunk in chunk_dict[list_sequences[i]].keys():
                        if chunk_dict[list_sequences[i]][chunk] > max:
                            max = chunk_dict[list_sequences[i]][chunk]
                            chunkk = chunk
                        elif chunk_dict[list_sequences[i]][chunk] > maxx:
                            maxx = chunk_dict[list_sequences[i]][chunk]
                            chunkkk = chunk
                seq_parente1 = chunkk
                seq_parente2 = chunkkk
            if seq_parente2 == "":
                max, chunkk = 0, ""
                for i in range(len(list_sequences)):
                    for chunk in chunk_dict[list_sequences[i]].keys():
                        if chunk_dict[list_sequences[i]][chunk] > max:
                            max = chunk_dict[list_sequences[i]][chunk]
                            chunkk = chunk
                seq_parente2 = chunkk
            for i in range(len(list_sequences)):
                list_identity = []
                count = 0
                for chunk in chunk_dict[list_sequences[i]].keys():
                    identity_p1 = get_identity([seq_parente1, chunk])
                    identity_p2 = get_identity([seq_parente2, chunk])
                    list_identity.append(identity_p1)
                    list_identity.append(identity_p2)
                    if identity_p1 == identity_p2:
                        count += 1
                if statistics.stdev(list_identity) < 5 or count < 2:
                    yield [s, c]
                    break
        except KeyError:
            pass

def abundance_greedy_clustering(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
    """List of OTU.
        :Parameters:
            amplicon_file: fasta.gz file
            minseqlen: minimum length of sequences
            mincount: minimum count to keep sequence
            chunk_size: size of each chunk
            kmer_size : size of kmer
        Returns: A list of OTU sequences in the format: [[sequence, count]]
    """
    pass

def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))

def write_OTU(OTU_list, output_file):
    """OTU list to a file with a specific format.
        :Parameters:
            otu_list: list of OTU
            output_file: name of output file
    """
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
    kmer = cut_kmer(seq,args.kmer_size)
    print(next(kmer))


if __name__ == '__main__':
    main()
