"""
Barcode class

Author: Ryan Wick
email: rrwick@gmail.com
"""

import math
from .misc import reverse_complement


class Barcode(object):

    def __init__(self, name, sequence, slope, intercept):
        self.name = name
        self.sequence = sequence
        self.rev_comp_sequence = reverse_complement(sequence)
        self.slope = slope
        self.intercept = intercept

    def get_sequence(self, direction):
        if direction == 'forward':
            return self.sequence
        else:  # direction == 'reverse'
            return self.rev_comp_sequence

    def get_perfect_score(self, match_score):
        return len(self.sequence) * match_score

    def get_expected_random_seq_score(self, random_seq_length):
        return self.slope * math.log10(math.log10(random_seq_length)) + self.intercept


# Predefined Oxford Nanopore native barcoding sequences.
# The slope and intercept values were determined by putting these sequences through the random
# alignment mode of nanopore_barcode_binner.py.
BARCODES = {'NB01': Barcode('NB01', 'GGTGCTG' 'AAGAAAGTTGTCGGTGTCTTTGTG' 'TTAACCT', 86.5, -32.4),
            'NB02': Barcode('NB02', 'GGTGCTG' 'TCGATTCCGTTTGTAGTCGTCTGT' 'TTAACCT', 85.3, -29.8),
            'NB03': Barcode('NB03', 'GGTGCTG' 'GAGTCTTGTGTCCCAGTTACCAGG' 'TTAACCT', 79.9, -26.5),
            'NB04': Barcode('NB04', 'GGTGCTG' 'TTCGGATTCTATCGTGTTTCCCTA' 'TTAACCT', 85.8, -30.9),
            'NB05': Barcode('NB05', 'GGTGCTG' 'CTTGTCCAGGGTTTGTGTAACCTT' 'TTAACCT', 86.2, -32.7),
            'NB06': Barcode('NB06', 'GGTGCTG' 'TTCTCGCAAAGGCAGAAAGTAGTC' 'TTAACCT', 80.5, -27.6),
            'NB07': Barcode('NB07', 'GGTGCTG' 'GTGTTACCGTGGGAATGAATCCTT' 'TTAACCT', 83.6, -30.7),
            'NB08': Barcode('NB08', 'GGTGCTG' 'TTCAGGGAACAAACCAAGTTACGT' 'TTAACCT', 80.9, -28.6),
            'NB09': Barcode('NB09', 'GGTGCTG' 'AACTAGGCACAGCGAGTCTTGGTT' 'TTAACCT', 78.8, -26.0),
            'NB10': Barcode('NB10', 'GGTGCTG' 'AAGCGTTGAAACCTTTGTCCTCTC' 'TTAACCT', 83.1, -29.7),
            'NB11': Barcode('NB11', 'GGTGCTG' 'GTTTCATCTATCGGAGGGAATGGA' 'TTAACCT', 82.4, -29.1),
            'NB12': Barcode('NB12', 'GGTGCTG' 'CAGGTAGAAAGAAGCAGAATCGGA' 'TTAACCT', 84.1, -29.4)}
