#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
(c) 2017 Brant Faircloth || http://faircloth-lab.org/
All rights reserved.

This code is distributed under a 3-clause BSD license. Please see
LICENSE.txt for more information.

Created on 19 January 2017 13:11 CST (-0600)
"""

#import pdb
import textwrap
import string
from Bio.Seq import Seq


def pretty_print(seq):
    '''function to textwrap, to default width, all output'''
    temp_seq = textwrap.wrap(seq)
    return '\n'.join(temp_seq)


def task1(dna):
    '''print the length of the sequence to the screen along w/ explanation'''
    # print DNA out to screen and make it pretty
    pretty_dna = pretty_print(dna)
    print("We are starting with a DNA string of: \n{}\n").format(pretty_dna)
    # now compute length and output that to the screen
    print("1. The length of the DNA string is {}\n".format(len(dna)))


def task2(dna):
    '''given a string representing dna, produce the RNA equivalent'''
    # convert DNA to RNA after first ensuring all chars are lowercase
    rna = dna.replace('t', 'u')
    # make that pretty (wrap the text)
    pretty_rna = pretty_print(rna)
    print("2. Converting the DNA strand to RNA gives us:\n{}\n".format(
        pretty_rna
    ))


def task3(dna):
    '''given a string representing DNA, produce the reverse complement of the
    string'''
    # setup a translation table to convert to complement
    dna_chars = 'acgt'
    complement_chars = 'tgca'
    translation_table = string.maketrans(dna_chars, complement_chars)
    # perform translation -> giving us complement
    comp_dna = string.translate(dna, translation_table)
    # reverse the complement
    revcomp = comp_dna[::-1]
    pretty_revcomp = pretty_print(revcomp)
    # print that out
    print("3. The reverse complement is:\n{}\n".format(pretty_revcomp))


def task4(dna):
    '''convert DNA string to codons, extract 13th and 14th codons, return
    codons'''
    # divide DNA into coding position w/ map, zip, and iter
    codons = map(''.join, zip(*[iter(dna)]*3))
    print("4. The 13th codon is [{}] and the 14th codon is [{}]\n".format(
        codons[12],
        codons[13]
    ))
    return codons


def task5(codons):
    '''convert codons into their AA sequence'''
    # read in the translation table for vert mtDNA
    with open('VertMitTransTable.txt') as infile:
        for line in infile:
            line = line.strip()
            if line.startswith('AAs'):
                aa = line.split('=')[1].strip()
            elif line.startswith('Starts'):
                starts = line.split('=')[1].strip()
            elif line.startswith('Base1'):
                b1 = line.split('=')[1].strip()
            elif line.startswith('Base2'):
                b2 = line.split('=')[1].strip()
            elif line.startswith('Base3'):
                b3 = line.split('=')[1].strip()
    # create a dict to hold codon to AA refs
    amino_acids = {}
    for k, v in enumerate(aa):
        # we can add an asterisk if also start codon
        #if starts[k] == 'M':
        #    v = v + '*'
        # create codon for each amino acid
        aa_name = '{}{}{}'.format(b1[k].lower(), b2[k].lower(), b3[k].lower())
        # because codons unique, use those as key with AA as value
        amino_acids[aa_name] = v
    # create an empty list to hold the translated sequence we're about to make
    translated = []
    # for each codon, add the AA to the list above
    for codon in codons:
        translated.append(amino_acids[codon])
    # join those
    translated_seq = ''.join(translated)
    return translated_seq


def validate_task5(dna, translated):
    '''Use BioPython to validate translation'''
    # drop off the 2 extra BP
    my_seq = Seq(dna[:-2])
    # use if then to validate our results == the one from BioPython
    if translated == str(my_seq.translate(table="Vertebrate Mitochondrial")):
        pass
    else:
        # if it doesn't raise an error
        raise IOError("There is a problem with the translation")


def pretty_print_task5(translated_seq):
    '''Pretty print (line wrap) translated sequence'''
    # send to pretty printing fxn
    pretty_translated_seq = pretty_print(translated_seq)
    print("5. The translated sequence is:\n{}\n".format(pretty_translated_seq))


def main():
    '''the main loop'''
    # read in the contents of the sequence file and strip off the crud
    dna = open('CodingSeq.txt').read().strip()
    # make sure all DNA chars are lowercase
    dna = dna.lower()
    # run task 1
    task1(dna)
    task2(dna)
    task3(dna)
    codons = task4(dna)
    translated = task5(codons)
    validate_task5(dna, translated)
    pretty_print_task5(translated)


if __name__ == '__main__':
    main()
