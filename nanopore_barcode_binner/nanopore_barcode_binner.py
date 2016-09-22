#!/usr/bin/env python3
"""
Nanopore barcode binner

Author: Ryan Wick
email: rrwick@gmail.com
"""

import argparse
import os
import sys
from .cpp_function_wrappers import barcode_alignment
from .misc import load_fasta_or_fastq
from .barcode import BARCODES
from .barcode_alignment_results import SingleBarcodeAlignmentResult, OverallBarcodeAlignmentResult
from .output_settings import OutputSettings
from .version import __version__


MATCH_SCORE = 3
MISMATCH_SCORE = -6
GAP_OPEN_SCORE = -5
GAP_EXTEND_SCORE = -2


def get_arguments():
    """
    Parse the command line arguments.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', '--reads', required=True,
                        help='FASTA or FASTQ of reads')
    parser.add_argument('-b', '--barcodes', default='NB01,NB02',
                        help='Names of barcodes, separated by commas')
    parser.add_argument('-o', '--out_dir', default='.',
                        help='Names of barcodes, separated by commas')
    parser.add_argument('-d', '--end_size', type=int, default=100,
                        help='Number of bases on the ends of reads within barcodes will be '
                             'searched for')
    parser.add_argument('--min_score', type=float, default=20.0,
                        help='Reads must have at least this score for their best barcode to be '
                             'classified')
    parser.add_argument('--gap', type=float, default=15.0,
                        help="A read's best score must be better than its second best score by "
                             "this much to be classified")
    return parser.parse_args()


def main():
    """
    Script execution starts here.
    """
    args = get_arguments()
    barcodes = [BARCODES[x] for x in args.barcodes.split(',')]

    output_settings = OutputSettings()

    # Load the reads and group them up by their type (2D/template/complement).
    reads = load_fasta_or_fastq(args.reads)
    grouped_reads = {}
    read_names = []
    for read_name, read_seq in reads:
        short_name = read_name.split('_Basecall_2D')[0]
        if short_name not in grouped_reads:
            read_names.append(short_name)
            grouped_reads[short_name] = NanoporeRead(short_name)
        grouped_reads[short_name].add_read_seq(read_name, read_seq)
    reads = [grouped_reads[x] for x in read_names]

    header_string = get_header(args, output_settings)
    print(header_string)

    for read in reads:
        output_line = [read.name]

        if output_settings.include_read_lengths:
            output_line += read.get_read_lengths()

        for barcode in barcodes:
            if output_settings.include_raw_score:
                output_line.append(barcode.get_perfect_score(MATCH_SCORE))

            read.align_barcode(barcode, args.end_size)
            output_line += read.get_output_line(barcode, output_settings)

        output_line += read.determine_best_barcode_match(barcodes, args.min_score, args.gap,
                                                         output_settings)
        print('\t'.join([str(x) for x in output_line]), flush=True)


def get_header(args, output_settings):
    header = ['Read name']
    if output_settings.include_read_lengths:
        header += ['2D read length', 'template read length', 'complement read length',
                   'longest read length']

    for barcode_name in args.barcodes.split(','):
        if output_settings.include_raw_score:
            header.append(barcode_name + ' perfect score')

        directions = ['2D, forward', '2D, reverse',
                      'template, forward', 'template, reverse',
                      'complement, forward', 'complement, reverse']
        for i, direction in enumerate(directions):

            if output_settings.include_raw_score:
                header.append(barcode_name + ' raw score (' + direction + ' barcode)')
                header.append(barcode_name + ' random score (' + direction + ' barcode)')
            if output_settings.include_adjusted_score:
                header.append(barcode_name + ' adjusted score (' + direction + ' barcode)')
            if output_settings.include_start_and_end:
                header.append(barcode_name + ' read start (' + direction + ' barcode)')
                header.append(barcode_name + ' read end (' + direction + ' barcode)')
            if output_settings.include_before_and_after_seq:
                header.append(barcode_name + ' before seq (' + direction + ' barcode)')
            if output_settings.include_seq:
                header.append(barcode_name + ' read seq (' + direction + ' barcode)')
            if output_settings.include_before_and_after_seq:
                header.append(barcode_name + ' after seq (' + direction + ' barcode)')

        if output_settings.include_best_scores:
            header.append(barcode_name + ' best start read barcode score')
            header.append(barcode_name + ' best end read barcode score')
            header.append(barcode_name + ' mean best barcode score')

    if output_settings.include_best_scores:
        header.append('Best barcode match')
    if output_settings.include_call_confidence:
        header.append('Call confidence')

    header += ['Barcode call']  # Always included
    header_string = '\t'.join(header)

    return header_string


def get_expected_end(read_type, direction):
    if direction == 'forward':
        if read_type == 'complement':
            return 'START'
        else:  # 2d and template reads
            return 'END'
    else:  # reverse barcodes
        if read_type == 'complement':
            return 'END'
        else:  # 2d and template reads
            return 'START'


class NanoporeRead(object):

    def __init__(self, name):
        self.name = name
        self.sequence_2d = ''
        self.sequence_template = ''
        self.sequence_complement = ''
        self.barcode_alignment_results = {}

    def add_read_seq(self, full_name, sequence):
        if full_name.endswith('_Basecall_2D_2d'):
            self.sequence_2d = sequence
        elif full_name.endswith('_Basecall_2D_template'):
            self.sequence_template = sequence
        elif full_name.endswith('_Basecall_2D_complement'):
            self.sequence_complement = sequence
        else:  # Assume 1D
            self.sequence_template = sequence

    def has_2d_sequence(self):
        return self.sequence_2d != ""

    def has_template_sequence(self):
        return self.sequence_template != ""

    def has_complement_sequence(self):
        return self.sequence_complement != ""

    def get_read_lengths(self):
        lengths = [len(self.sequence_2d), len(self.sequence_template),
                   len(self.sequence_complement)]
        return [x if x > 0 else '' for x in lengths] + [max(lengths)]

    def get_sequence(self, read_type):
        if read_type == '2d':
            return self.sequence_2d
        elif read_type == 'template':
            return self.sequence_template
        else:  # read_type == 'complement'
            return self.sequence_complement

    def has_type(self, read_type):
        if read_type == '2d' and self.has_2d_sequence():
            return True
        elif read_type == 'template' and self.has_template_sequence():
            return True
        elif read_type == 'complement' and self.has_complement_sequence():
            return True
        return False

    def align_barcode(self, barcode, end_size):
        all_results = []

        for read_type in ['2d', 'template', 'complement']:

            if not self.has_type(read_type):
                all_results += [None, None]
                continue

            read_seq = self.get_sequence(read_type)

            for direction in ['forward', 'reverse']:
                expected_end = get_expected_end(read_type, direction)

                # Trim the read sequence to the region worth examining.
                region_size = min(len(read_seq), end_size)
                if expected_end == 'START':
                    trimmed_read_seq = read_seq[:region_size]
                else:  # expected_end == 'END'
                    trimmed_read_seq = read_seq[-region_size:]

                cpp_result = barcode_alignment(trimmed_read_seq, barcode.get_sequence(direction),
                                               MATCH_SCORE, MISMATCH_SCORE, GAP_OPEN_SCORE,
                                               GAP_EXTEND_SCORE)
                all_results.append(SingleBarcodeAlignmentResult(cpp_result, barcode, expected_end,
                                                                trimmed_read_seq, len(read_seq),
                                                                end_size, MATCH_SCORE))
        self.barcode_alignment_results[barcode.name] = OverallBarcodeAlignmentResult(all_results)

    def get_output_line(self, barcode, output_settings):
        output_line = self.barcode_alignment_results[barcode.name].get_output_line(output_settings)

        return output_line

    def determine_best_barcode_match(self, barcodes, min_score, min_gap, output_settings):
        output_line = []

        scores = [(b.name, self.barcode_alignment_results[b.name].mean_best_score)
                  for b in barcodes]
        scores = sorted(scores, key=lambda x: x[1], reverse=True)
        best_name = scores[0][0]
        best_score = scores[0][1]

        if output_settings.include_best_scores:
            output_line.append(best_name)

        # If the best score is too low, we don't make a barcode call.
        if best_score < min_score:
            final_call = 'None'
            confidence = 0.0

        else:  # Best score is decent
            confidence = 100 * (best_score - min_score) / (100.0 - min_score)

            # If there's only one possible barcode (an odd situation), then it gets the call.
            if len(scores) == 1:
                final_call = best_name

            # If there are multiple possible barcodes, then the best must be sufficiently better
            # than the second best to be called.
            else:
                second_best_score = scores[1][1]
                gap = best_score - second_best_score
                if gap >= min_gap:
                    final_call = best_name
                    confidence *= min(1.0, gap / 100.0)
                else:
                    final_call = 'None'
                    confidence = 0.0

        if output_settings.include_call_confidence:
            output_line.append(confidence)
        output_line.append(final_call)

        return output_line
