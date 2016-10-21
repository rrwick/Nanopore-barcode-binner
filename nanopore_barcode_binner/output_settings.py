"""
Controls which columns appear in the output table

Author: Ryan Wick
email: rrwick@gmail.com
"""


class OutputSettings(object):

    def __init__(self, verbosity):
        self.read_name = verbosity > 0
        self.read_lengths = verbosity > 2
        self.all_scores = verbosity > 2
        self.start_and_end = verbosity > 2
        self.seq = verbosity > 2
        self.best_scores = verbosity > 1
        self.call_confidence = verbosity > 2
        self.final_call = verbosity > 0
