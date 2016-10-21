"""
Alignment results between one read and one barcode

Author: Ryan Wick
email: rrwick@gmail.com
"""


class OverallBarcodeAlignmentResult(object):

    def __init__(self, single_results):
        assert len(single_results) == 6
        self.single_results = single_results

        forward_2d_score = get_adjusted_score(single_results[0])
        reverse_2d_score = get_adjusted_score(single_results[1])
        forward_template_score = get_adjusted_score(single_results[2])
        reverse_template_score = get_adjusted_score(single_results[3])
        forward_complement_score = get_adjusted_score(single_results[4])
        reverse_complement_score = get_adjusted_score(single_results[5])

        self.best_start_score = max_with_nones(reverse_2d_score, reverse_template_score,
                                               forward_complement_score)
        self.best_end_score = max_with_nones(forward_2d_score, forward_template_score,
                                             reverse_complement_score)
        self.mean_best_score = (self.best_start_score + self.best_end_score) / 2.0

    def get_output_line(self, output_settings, reads_are_2d):
        empty_output_line = []
        if output_settings.all_scores:
            empty_output_line += ['']
        if output_settings.start_and_end:
            empty_output_line += ['', '']
        if output_settings.seq:
            empty_output_line += ['']

        output_line = []
        if reads_are_2d:
            single_results = self.single_results
        else:
            single_results = self.single_results[2:4]
        for single_result in single_results:
            if single_result is None:
                output_line += empty_output_line
            else:
                output_line += single_result.get_output_line(output_settings)
        return output_line


class SingleBarcodeAlignmentResult(object):

    def __init__(self, cpp_result, barcode, expected_end, trimmed_read_seq, trimmed_length,
                 match_score):
        parts = cpp_result.split(',')

        self.read_start_pos = int(parts[0])
        self.read_end_pos = int(parts[1])
        self.before_hit_sequence = trimmed_read_seq[self.read_start_pos-20:self.read_start_pos]
        self.hit_sequence = trimmed_read_seq[self.read_start_pos:self.read_end_pos+1]
        self.after_hit_sequence = trimmed_read_seq[self.read_end_pos+1:self.read_end_pos+21]

        if expected_end == 'END':
            self.read_start_pos += trimmed_length
            self.read_end_pos += trimmed_length

        self.raw_score = int(parts[4])
        self.expected_score_for_random_seq = barcode.get_expected_random_seq_score(len(trimmed_read_seq))
        self.adjusted_score = 100.0 * (self.raw_score - self.expected_score_for_random_seq) / \
            (barcode.get_perfect_score(match_score) - self.expected_score_for_random_seq)

    def get_output_line(self, output_settings):
        output_line = []
        if output_settings.all_scores:
            output_line.append(self.adjusted_score)
        if output_settings.start_and_end:
            output_line.append(self.read_start_pos)
            output_line.append(self.read_end_pos)
        if output_settings.seq:
            output_line.append(self.hit_sequence)
        return output_line


def get_adjusted_score(single_result):
    if single_result is None:
        return None
    else:
        return single_result.adjusted_score


def max_with_nones(*args):
    return max(x for x in args if x is not None)
