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
#from collections import Counter
# https://github.com/briney/nwalign3
# ftp://ftp.ncbi.nih.gov/blast/matrices/
#import nwalign3 as nw

__author__ = "Colin Davidson"
__copyright__ = "EISTI"
__credits__ = ["Colin DAVIDSON"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Colin DAVIDSON"
__email__ = "your@email.fr"
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
            amplicon_file: fasta file to read
            minseqlen: min length of a sequence
        Returns: A generator of sequences
    """
    if isfile(amplicon_file):
        with gzip.open(amplicon_file, 'r') as filename:
            content = []
            for amplicon in filename.read().decode('utf-8').replace('\n','').split('>'):
                if len(amplicon) >= minseqlen:
                    for i in amplicon[::-1]:
                        if i not in 'ATGC':
                            content.append(amplicon.split(i)[-1])
                            break
            for amplicon in content:
                yield amplicon


def dereplication_fulllength(amplicon_file, minseqlen, mincount):
    """Creates a sequence, count generator.
        :Parameters:
            amplicon_file: fasta file to read
            minseqlen: min length of a sequence
            mincount: minimal value to retain this amplicon
        Returns: A generator of [sequence, count]
    """
    dict_amplicons = {}
    generator_amplicons = read_fasta(amplicon_file, minseqlen)

    for amplicon in generator_amplicons:
        if amplicon in dict_amplicons:
            dict_amplicons[amplicon] += 1
        else:
            dict_amplicons[amplicon] = 1

    dict_amplicons = {amplicon:dict_amplicons[amplicon] for amplicon in
    sorted(list(dict_amplicons.keys())) if dict_amplicons[amplicon] >= mincount}

    for amplicon, occurence in dict_amplicons.items():
        yield [amplicon, occurence]


def get_chunks(sequence, chunk_size):
    """Transforms a sequence into chunks.
        :Parameters:
            sequence: str of sequence
            chunk_size: int for each chunk's size
        Returns: An array of chunks
    """
    if chunk_size > len(sequence)/4:
        print('ValueError: chunk_size ({}) > len(sequence)({})/4'.format(chunk_size,len(sequence)))
        raise ValueError
    return [sequence[index:index+chunk_size] for index
    in range(0,len(sequence)-len(sequence)%chunk_size,chunk_size)]

def get_unique(ids):
    """ Gets unique ids
    """
    return {}.fromkeys(ids).keys()


def common(lst1, lst2):
    """ Checks common items
    """
    return list(set(lst1) & set(lst2))

def cut_kmer(sequence, kmer_size):
    """ Reads a sequence and yields the k-mer
        :Parameters:
            sequence : Sequence (str)
            kmer_size : Size of kmer (int)
        Returns: yields the kmers
    """
    for i in range(len(sequence) - kmer_size + 1):
        yield sequence[i:i+kmer_size]

def get_identity(alignment_list):
    """Computes identity between two alignements.
        :Parameters:
            alignment_list: list of two alignments to compare
        Returns: A percentage of identity
    """
    return (sum([[0,1][alignment_list[0][index] == alignment_list[1][index]] for index in
    range(min(len(alignment_list[0]),
    len(alignment_list[1])))])/min(len(alignment_list[0]),
    len(alignment_list[1])))*100


def chimera_removal(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
    """Creates a chimera free sequence generator.
        :Parameters:
            amplicon_file: fasta file to read
            minseqlen: min length of a sequence
            mincount: minimal value to retain this amplicon
            chunk_size: int for each chunk's size
            kmer_size : Size of kmer (int)
        Returns: A generator of [sequence, count]
    """
    for sequence, count in dereplication_fulllength(amplicon_file, minseqlen, mincount):
        sequence_dict = {}
        try:
            chunks = get_chunks(sequence, chunk_size)[:4]
            for chunk in chunks:
                dict_id = {}
                altered_sequence = sequence.replace(chunk,'')
                segments = cut_kmer(altered_sequence, len(chunk))
                for segment in segments:
                    identity = get_identity([segment, chunk])
                    if len(dict_id) < 8 and identity > 0:
                        dict_id[segment] = identity
                    min_segment_identity = min(dict_id, default=666666666666666,
                    key=lambda x: dict_id[x])
                    if identity > dict_id[min_segment_identity]:
                        del dict_id[min_segment_identity]
                        dict_id[segment] = identity
                sorted_segments = sorted(dict_id, key=lambda x: dict_id[x])
                sequence_dict[chunk] = {segment:dict_id[segment] for segment in sorted_segments}
            list_sequence_dict = list(sequence_dict.keys())
            parent1, parent2 = '', ''
            for segment1 in sequence_dict[list_sequence_dict[0]].keys():
                counter = 1
                for index in range(1, len(list_sequence_dict)):
                    for segment2 in sequence_dict[list_sequence_dict[index]].keys():
                        if segment1 == segment2:
                            counter += 1
                        if counter == len(list_sequence_dict):
                            if parent1 != '' and segment1 == segment2:
                                parent1 = segment1
                            elif parent2 != '' and segment1 == segment2:
                                parent2 = segment2
                        break
            if parent1 == '':
                max1, max2, segment1, segment2 = 0, 0, '', ''
                for index in range(len(list_sequence_dict)):
                    for segment in sequence_dict[list_sequence_dict[index]].keys():
                        if sequence_dict[list_sequence_dict[index]][segment] > max1:
                            max1 = sequence_dict[list_sequence_dict[index]][segment]
                            segment1 = segment
                        elif sequence_dict[list_sequence_dict[index]][segment] > max2:
                            max2 = sequence_dict[list_sequence_dict[index]][segment]
                            segment2 = segment
                parent1 = segment1
                parent2 = segment2
            if parent2 == '':
                max2, segment2 = 0, ''
                for index in range(len(list_sequence_dict)):
                    for segment in sequence_dict[list_sequence_dict[index]].keys():
                        if sequence_dict[list_sequence_dict[index]][segment] > max2:
                            max2 = sequence_dict[list_sequence_dict[index]][segment]
                            segment2 = segment
                parent2 = segment2

            for index in range(len(list_sequence_dict)):
                identities = []
                counter = 0
                for segment in sequence_dict[list_sequence_dict[index]].keys():
                    identities += [get_identity([parent1,segment]),get_identity([parent2,segment])]
                    if identities[-2] == identities[-1]:
                        counter += 1
                if statistics.stdev(identities) < 5 or counter < 2:
                    yield [sequence, count]
                    break
        except ValueError:
            pass


def abundance_greedy_clustering(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
    """Regroups all OTU sequences into a list.
        :Parameters:
            amplicon_file: fasta file to read
            minseqlen: min length of a sequence
            mincount: minimal value to retain this amplicon
            chunk_size: int for each chunk's size
            kmer_size : Size of kmer (int)
        Returns: A list of OTU sequences [[sequence, count]]
    """
    otu_list = []
    list_chimera_removal = list(chimera_removal(amplicon_file,
    minseqlen, mincount, chunk_size, kmer_size))
    for index in range(len(list_chimera_removal)):
        add = True
        sequence1, count1 = list_chimera_removal[index]
        for jndex in range(index + 1, len(list_chimera_removal)):
            sequence2, count2 = list_chimera_removal[jndex]
            if sequence1 != sequence2 and get_identity([sequence1, sequence2])>97 and count2>count1:
                add = False
                break
        if add:
            otu_list.append([sequence1, count1])
    return otu_list


def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))

def write_OTU(otu_list, output_file):
    """Writes an OTU list under a specific format.
        :Parameters:
            otu_list: list of OTU sequences [[sequence, count]]
            output_file: name of file to save
    """
    with open(output_file, 'w') as file:
        i = 1
        for otu in otu_list:
            text = ">OTU_{} occurence:{}\n{}\n".format(i, otu[1], fill(otu[0]))
            if i == len(otu_list):
                text = text[:-1]
            file.write(text)
            i += 1
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
    chunk_size = args.chunk_size
    kmer_size = args.kmer_size

    # Get OTU List
    otu_list = abundance_greedy_clustering(amplicon_file,minseqlen, mincount, chunk_size, kmer_size)

    # Write OTU List
    write_OTU(otu_list, "results.fna")


if __name__ == '__main__':
    main()
