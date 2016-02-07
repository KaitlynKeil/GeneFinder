# -*- coding: utf-8 -*-
"""
YOUR HEADER COMMENT HERE

@author: YOUR NAME HERE

"""

import random
from amino_acids import aa, codons, aa_table   # you may find these useful
from load import load_seq


def shuffle_string(s):
    """Shuffles the characters in the input string
        NOTE: this is a helper function, you do not
        have to modify this in any way """
    return ''.join(random.sample(s, len(s)))

# YOU WILL START YOUR IMPLEMENTATION FROM HERE DOWN ###


def get_complement(nucleotide):
    """ Returns the complementary nucleotide

        nucleotide: a nucleotide (A, C, G, or T) represented as a string
        returns: the complementary nucleotide
        # make sure all complements work and that non-nucleotides fail
    >>> get_complement('A')
    'T'
    >>> get_complement('C')
    'G'
    >>> get_complement('G')
    'C'
    >>> get_complement('T')
    'A'
    >>> get_complement('J')
    'Not a DNA nucleotide'
    """
    if nucleotide == 'A':
        return 'T'
    elif nucleotide == 'T':
        return 'A'
    elif nucleotide == 'G':
        return 'C'
    elif nucleotide == 'C':
        return 'G'
    else:
        return 'Not a DNA nucleotide'

def get_reverse_complement(dna):
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence

        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string
    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'
    >>> get_reverse_complement("CCGOGTTCA")
    'Strange mutation in the DNA string; check for data errors'
    """
    complimentString = ''
    for x in dna:
        comp = get_complement(x)
        complimentString += comp
    if "Not a DNA nucleotide" in complimentString:
        return 'Strange mutation in the DNA string; check for data errors'
    return complimentString[::-1]

def rest_of_ORF(dna):
    """ Takes a DNA sequence that is assumed to begin with a start
        codon and returns the sequence up to but not including the
        first in frame stop codon.  If there is no in frame stop codon,
        returns the whole string.

        dna: a DNA sequence
        returns: the open reading frame represented as a string
    >>> rest_of_ORF("ATGTGAA")
    'ATG'
    >>> rest_of_ORF("ATGAGATAGG")
    'ATGAGA'
    >>> rest_of_ORF("ATGAGATA")
    'ATGAGATA'
    """

    end_condons = ['TGA', 'TAA', 'TAG']
    ORF = dna[0:3]
    x = 1
    while x <= len(dna)/3:
        current_index = x*3
        checked_codon = dna[current_index:current_index+3]
        if checked_codon in end_condons:
            return ORF
        else:
            ORF += checked_codon
        x += 1

    return ORF

def find_all_ORFs_oneframe(dna):
    """ Finds all non-nested open reading frames in the given DNA
        sequence and returns them as a list.  This function should
        only find ORFs that are in the default frame of the sequence
        (i.e. they start on indices that are multiples of 3).
        By non-nested we mean that if an ORF occurs entirely within
        another ORF, it should not be included in the returned list of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC']
    >>> find_all_ORFs_oneframe("GCAATGCATGAATGTAGATAGATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC']
    """
    checker = 0 # keeps track of which characters have been looked at
    number_orf = 0 # gives the index for the ORFs
    oneframe = []
    start_codon = 'ATG'
    while checker <= len(dna):
        if dna[checker:checker + 3] == start_codon:
            ORFs = rest_of_ORF(dna[checker:])
            checker += len(ORFs)
            oneframe.append(ORFs)
        checker += 3
    return oneframe
    
def find_all_ORFs(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence in
        all 3 possible frames and returns them as a list.  By non-nested we
        mean that if an ORF occurs entirely within another ORF and they are
        both in the same frame, it should not be included in the returned list
        of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

    >>> find_all_ORFs("ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']
    """
    ORF_list = []
    for x in range(3):
        current_dna = dna[x:]
        current_orf_list = find_all_ORFs_oneframe(current_dna)
        ORF_list.extend(current_orf_list)
    return ORF_list
    
def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    """
    all_ORFs = []
    reverse_complement = get_reverse_complement(dna) # stores the reverse as a string
    dna_forward_ORFs = find_all_ORFs(dna)
    all_ORFs.extend(dna_forward_ORFs)
    rc_ORFs = find_all_ORFs(reverse_complement)
    all_ORFs.extend(rc_ORFs)
    return all_ORFs # a list of all the ORFs

def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    """
    possibilities = find_all_ORFs_both_strands(dna)
    if len(possibilities) > 0:
        longest_strand = possibilities[0]
        for x in range(len(possibilities) - 1):
            if len(possibilities[x + 1]) > len(longest_strand):
                longest_strand = possibilities[x + 1]
        return longest_strand
    else:
        return '' # if it does not find any possibilities, it will return an empty string

def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF """
    i = 0
    possible_longest_ORFs = []
    while i < num_trials:
        shuffled_dna = shuffle_string(dna)
        possible_longest_ORFs.append(len(longest_ORF(shuffled_dna)))
        i += 1
    possible_longest_ORFs.sort()
    return possible_longest_ORFs[-1]

def coding_strand_to_AA(dna):
    """ Computes the Protein encoded by a sequence of DNA.  This function
        does not check for start and stop codons (it assumes that the input
        DNA sequence represents an protein coding region).

        dna: a DNA sequence represented as a string
        returns: a string containing the sequence of amino acids encoded by the
                 the input DNA fragment

        >>> coding_strand_to_AA("ATGCGA")
        'MR'
        >>> coding_strand_to_AA("ATGCCCGCTTT")
        'MPA'
        >>> coding_strand_to_AA("AT")
        ''
    """
    import amino_acids
    i = 0
    aa_chain = ''
    while i < len(dna):
        current_codon = dna[i:i+3]
        if len(current_codon) % 3 != 0:
            return aa_chain
        amino_acid = amino_acids.aa_table[current_codon]
        aa_chain += amino_acid
        i +=3
    return aa_chain

def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """
    threshold = longest_ORF_noncoding(dna, 1500)
    ORFs = find_all_ORFs_both_strands(dna)
    sequences = []
    for i in range(len(ORFs)):
        if len(ORFs[i]) >= threshold:
            sequences.append(coding_strand_to_AA(ORFs[i]))
    return sequences
    
if __name__ == "__main__":
    import doctest
    #doctest.testmod()
    doctest.run_docstring_examples(coding_strand_to_AA, globals())
