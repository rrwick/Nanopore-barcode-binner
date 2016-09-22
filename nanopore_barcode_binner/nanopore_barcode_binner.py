#!/usr/bin/env python3
"""
Nanopore barcode binner

Author: Ryan Wick
email: rrwick@gmail.com
"""

import argparse
import os
import collections
from .cpp_function_wrappers import barcode_alignment
from .misc import load_fasta_or_fastq, strip_read_extensions
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
                        help='Output directory for binned reads')
    parser.add_argument('--best', action='store_true',
                        help='Only output the best type for a read')
    parser.add_argument('-d', '--end_size', type=int, default=100,
                        help='Number of bases on the ends of reads within barcodes will be '
                             'searched for')
    parser.add_argument('--min_score', type=float, default=20.0,
                        help='Reads must have at least this score for their best barcode to be '
                             'classified')
    parser.add_argument('--gap', type=float, default=15.0,
                        help="A read's best score must be better than its second best score by "
                             "this much to be classified")
    parser.add_argument('--version', action='version', version=__version__)
    return parser.parse_args()


def main():
    """
    Script execution starts here.
    """
    args = get_arguments()
    barcodes = [BARCODES[x] for x in args.barcodes.split(',')]

    output_settings = OutputSettings()

    # Load the reads and group them up by their type (2D/template/complement).
    reads, read_type = load_fasta_or_fastq(args.reads)
    grouped_reads = {}
    read_names = []
    for read in reads:
        read_name = read[0]
        short_name = read_name.split('_Basecall_2D')[0]
        if short_name not in grouped_reads:
            read_names.append(short_name)
            grouped_reads[short_name] = NanoporeRead(short_name)
        grouped_reads[short_name].add_read(read)
    reads = [grouped_reads[x] for x in read_names]

    header_string = get_header(args, output_settings)
    print(header_string)

    reads_by_barcode = collections.defaultdict(list)
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
        reads_by_barcode[read.final_call].append(read)

    out_dir = os.path.abspath(args.out_dir)
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)

    barcode_names = [b.name for b in barcodes] + ['None']
    for barcode_name in barcode_names:
        if not reads_by_barcode[barcode_name]:
            continue
        binned_reads_filename = strip_read_extensions(os.path.basename(args.reads)) + '_' + \
            barcode_name
        if read_type == 'FASTA':
            binned_reads_filename += '.fasta'
        else:  # read_type == 'FASTQ'
            binned_reads_filename += '.fastq'
        binned_reads_path = os.path.join(out_dir, binned_reads_filename)
        with open(binned_reads_path, 'wt') as out_reads:
            if read_type == 'FASTA':
                for read in reads_by_barcode[barcode_name]:
                    out_reads.write(read.get_fasta_lines(args.best))
            else:  # read_type == 'FASTQ'
                for read in reads_by_barcode[barcode_name]:
                    out_reads.write(read.get_fastq_lines(args.best))


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
        self.read_2d = None
        self.read_template = None
        self.read_complement = None
        self.barcode_alignment_results = {}
        self.final_call = 'None'

    def add_read(self, read):
        full_name = read[0]
        if full_name.endswith('_Basecall_2D_2d'):
            self.read_2d = read
        elif full_name.endswith('_Basecall_2D_template'):
            self.read_template = read
        elif full_name.endswith('_Basecall_2D_complement'):
            self.read_complement = read
        else:  # Assume 1D
            self.read_template = read

    def has_2d_sequence(self):
        return self.read_2d is not None

    def has_template_sequence(self):
        return self.read_template is not None

    def has_complement_sequence(self):
        return self.read_complement is not None

    def get_read_lengths(self):
        max_length = 0
        lengths = []
        for read in [self.read_2d, self.read_template, self.read_complement]:
            if read is not None:
                length = len(read[1])
                lengths.append(length)
                max_length = max(max_length, length)
            else:
                lengths.append('')
        return lengths + [max_length]

    def get_sequence(self, read_type):
        if read_type == '2d':
            return self.read_2d[1]
        elif read_type == 'template':
            return self.read_template[1]
        else:  # read_type == 'complement'
            return self.read_complement[1]

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
        return self.barcode_alignment_results[barcode.name].get_output_line(output_settings)

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
            self.final_call = 'None'
            confidence = 0.0

        else:  # Best score is decent
            confidence = 100 * (best_score - min_score) / (100.0 - min_score)

            # If there's only one possible barcode (an odd situation), then it gets the call.
            if len(scores) == 1:
                self.final_call = best_name

            # If there are multiple possible barcodes, then the best must be sufficiently better
            # than the second best to be called.
            else:
                second_best_score = scores[1][1]
                gap = best_score - second_best_score
                if gap >= min_gap:
                    self.final_call = best_name
                    confidence *= min(1.0, gap / 100.0)
                else:
                    self.final_call = 'None'
                    confidence = 0.0

        if output_settings.include_call_confidence:
            output_line.append(confidence)
        output_line.append(self.final_call)

        return output_line

    def get_fasta_lines(self, only_best):
        if only_best:
            return self.get_fasta_lines_one_type(self.get_best_read())
        else:
            fasta_lines = ''
            if self.has_2d_sequence():
                fasta_lines += self.get_fasta_lines_one_type(self.read_2d)
            if self.has_complement_sequence():
                fasta_lines += self.get_fasta_lines_one_type(self.read_complement)
            if self.has_template_sequence():
                fasta_lines += self.get_fasta_lines_one_type(self.read_template)
            return fasta_lines

    @staticmethod
    def get_fasta_lines_one_type(read):
        fasta_lines = '>'
        fasta_lines += read[2]
        fasta_lines += '\n'
        fasta_lines += read[1]
        fasta_lines += '\n'
        return fasta_lines

    def get_fastq_lines(self, only_best):
        if only_best:
            return self.get_fastq_lines_one_type(self.get_best_read())
        else:
            fastq_lines = ''
            if self.has_2d_sequence():
                fastq_lines += self.get_fastq_lines_one_type(self.read_2d)
            if self.has_complement_sequence():
                fastq_lines += self.get_fastq_lines_one_type(self.read_complement)
            if self.has_template_sequence():
                fastq_lines += self.get_fastq_lines_one_type(self.read_template)
            return fastq_lines

    @staticmethod
    def get_fastq_lines_one_type(read):
        fastq_lines = '@'
        fastq_lines += read[4]
        fastq_lines += '\n'
        fastq_lines += read[1]
        fastq_lines += '\n'
        fastq_lines += read[2]
        fastq_lines += '\n'
        fastq_lines += read[3]
        fastq_lines += '\n'
        return fastq_lines

    def get_best_type(self):
        if self.has_2d_sequence():
            return '2d'
        elif not self.has_complement_sequence():
            return 'template'
        elif not self.has_template_sequence():
            return 'complement'

        # If there is both template and complement, but not 2d, then we choose the best using the
        # average Phred score.
        else:
            template_err_rate = self.est_error_rate(self.read_template)
            complement_err_rate = self.est_error_rate(self.read_complement)
            if template_err_rate <= complement_err_rate:
                return 'template'
            else:
                return 'complement'

    def get_best_read(self):
        best_type = self.get_best_type()
        if best_type == '2d':
            return self.read_2d
        elif best_type == 'template':
            return self.read_template
        else:  # best_type == 'complement'
            return self.read_complement

    @staticmethod
    def est_error_rate(read):
        """
        Returns an error rate estimate using the Phred quality scores.
        """
        # If there are no quality scores, we cannot estimate the error rate.
        if len(read) < 4:
            return 1.0
        qualities = read[3]

        error_count = 0.0
        for score in qualities:
            phred = ord(score) - 33
            error_count += 10.0 ** (-phred / 10.0)
        return error_count / len(qualities)
