"""
Controls which columns appear in the output table

Author: Ryan Wick
email: rrwick@gmail.com
"""


class OutputSettings(object):

    def __init__(self,
                 include_read_lengths=False, include_raw_score=False,
                 include_adjusted_score=False, include_start_and_end=False,
                 include_seq=False, include_before_and_after_seq=False,
                 include_best_scores=False, include_call_confidence=False):
        self.include_read_lengths = include_read_lengths
        self.include_raw_score = include_raw_score
        self.include_adjusted_score = include_adjusted_score
        self.include_start_and_end = include_start_and_end
        self.include_seq = include_seq
        self.include_before_and_after_seq = include_before_and_after_seq
        self.include_best_scores = include_best_scores
        self.include_call_confidence = include_call_confidence
