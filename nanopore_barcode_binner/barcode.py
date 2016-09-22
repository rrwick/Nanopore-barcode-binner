"""
Barcode class

Author: Ryan Wick
email: rrwick@gmail.com
"""

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


# Predefined Oxford Nanopore native barcoding sequences.
BARCODES = {'NB01': Barcode('NB01', 'AAGAAAGTTGTCGGTGTCTTTGTG',
                            71.4332584271663, -22.5656089159852),
            'NB02': Barcode('NB02', 'TCGATTCCGTTTGTAGTCGTCTGT',
                            67.1855037419996, -18.0468313716098),
            'NB03': Barcode('NB03', 'GAGTCTTGTGTCCCAGTTACCAGG',
                            63.5927250003896, -15.9980232349283),
            'NB04': Barcode('NB04', 'TTCGGATTCTATCGTGTTTCCCTA',
                            66.8473788581195, -18.5724631096649),
            'NB05': Barcode('NB05', 'CTTGTCCAGGGTTTGTGTAACCTT',
                            66.5864343955349, -18.8843342040538),
            'NB06': Barcode('NB06', 'TTCTCGCAAAGGCAGAAAGTAGTC',
                            65.4698361318286, -17.6317088824286),
            'NB07': Barcode('NB07', 'GTGTTACCGTGGGAATGAATCCTT',
                            64.9387340902843, -17.4358282814267),
            'NB08': Barcode('NB08', 'TTCAGGGAACAAACCAAGTTACGT',
                            65.5252104348060, -18.0882967050121),
            'NB09': Barcode('NB09', 'AACTAGGCACAGCGAGTCTTGGTT',
                            63.9800483064580, -15.9533793977375),
            'NB10': Barcode('NB10', 'AAGCGTTGAAACCTTTGTCCTCTC',
                            66.5716250549626, -18.8997212170206),
            'NB11': Barcode('NB11', 'GTTTCATCTATCGGAGGGAATGGA',
                            66.8540339246670, -18.6488190718658),
            'NB12': Barcode('NB12', 'CAGGTAGAAAGAAGCAGAATCGGA',
                            70.1581047865076, -20.4691636410478)}
