#!/usr/bin/env python3
"""
Nanopore barcode binner

Author: Ryan Wick
email: rrwick@gmail.com
"""

import argparse
import os
import sys
import collections
import random
import math
from .cpp_function_wrappers import barcode_alignment
from .misc import load_fasta_or_fastq, strip_read_extensions, quit_with_error, get_random_sequence
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
    parser.add_argument('--version', action='version', version=__version__)

    binning_mode_settings = parser.add_argument_group('Binning mode',
                                                      'When these options are used, the program '
                                                      'runs in its normal mode, sorting reads into '
                                                      'bins based on barcode alignments. It will '
                                                      'create FASTQ files containing the reads '
                                                      'and print a table of results to stdout.')
    binning_mode_settings.add_argument('-r', '--reads',
                                       help='FASTA or FASTQ of reads (required)')
    binning_mode_settings.add_argument('-b', '--barcodes', default='NB01,NB02',
                                       help='Names of barcodes, separated by commas (required)')
    binning_mode_settings.add_argument('-o', '--out_dir', default='.',
                                       help='Output directory for binned reads (default = current '
                                            'directory)')
    binning_mode_settings.add_argument('--best', action='store_true',
                                       help='Only save the best sequence for a read (default = '
                                            'save all sequences)')
    binning_mode_settings.add_argument('--no_trim', action='store_true',
                                       help='Leave the barcodes on the reads')
    binning_mode_settings.add_argument('--end_size', type=int, default=200,
                                       help='Number of bases on the ends of reads within barcodes '
                                            'will be searched for (default = 200)')
    binning_mode_settings.add_argument('--min_score', type=float, default=20.0,
                                       help='Reads must have at least this adjusted score for '
                                            'their best barcode to be classified (default = 20)')
    binning_mode_settings.add_argument('--gap', type=float, default=10.0,
                                       help="A read's best score must be better than its second "
                                            "best score by this much to be classified (default = "
                                            "10)")
    binning_mode_settings.add_argument('--verbosity', type=int, required=False, default=1,
                                       help='Level of information in stdout table (0 to 3, '
                                            'default: 1)')

    random_mode_settings = parser.add_argument_group('Random test-alignment mode',
                                                     'When these options are used, the program '
                                                     'runs in a special mode to test barcode '
                                                     'alignments to random sequences.')
    random_mode_settings.add_argument('--random_test_barcode', type=str, default=None,
                                      help='The barcode used in random alignment')
    random_mode_settings.add_argument('--random_test_count', type=int, default=None,
                                      help='When in random-test mode, this is the number of random '
                                           'alignments conducted')
    random_mode_settings.add_argument('--random_output_all', action='store_true',
                                      help='Output all random alignment scores (default = only '
                                           'output the regression results)')

    args = parser.parse_args()

    if args.random_test_barcode or args.random_test_count:
        if not args.random_test_barcode or not args.random_test_count:
            quit_with_error('When running in random-test mode, both --random_test_barcode and '
                            '--random_test_count are required')
    else:
        if not args.reads:
            quit_with_error('--reads are required')
        if not args.reads:
            quit_with_error('--barcodes are required')

    return args


def main():
    """
    Script execution starts here.
    """
    args = get_arguments()
    if args.random_test_barcode:
        random_mode(args.random_test_barcode, args.random_test_count, args.random_output_all)
        sys.exit(0)

    barcodes = [BARCODES[x] for x in args.barcodes.split(',')]

    output_settings = OutputSettings(args.verbosity)

    # Load the reads and group them up by their type (2D/template/complement).
    reads, read_type = load_fasta_or_fastq(args.reads)
    grouped_reads = {}
    read_names = []
    for read in reads:
        read_name = read[0]
        short_name = read_name.split('_Basecall')[0]
        if short_name not in grouped_reads:
            read_names.append(short_name)
            grouped_reads[short_name] = NanoporeRead(short_name)
        grouped_reads[short_name].add_read(read)
    reads = [grouped_reads[x] for x in read_names]

    # If all of the read groups have only a single template read, then we assume we're dealing
    # with a 1D dataset.
    reads_are_2d = any(x.has_2d_sequence() or x.has_complement_sequence() for x in reads)

    if args.verbosity > 0:
        print(get_header(args, output_settings, reads_are_2d))

    reads_by_barcode = collections.defaultdict(list)

    for read in reads:
        for barcode in barcodes:
            read.align_barcode(barcode, args.end_size)
        read.determine_best_barcode_match(barcodes, args.min_score, args.gap)
        reads_by_barcode[read.final_call].append(read)
        if args.verbosity > 0:
            output_line = read.get_output_line(output_settings, barcodes, reads_are_2d)
            print('\t'.join([convert_to_string(x) for x in output_line]), flush=True)

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
            for read in reads_by_barcode[barcode_name]:
                if not args.no_trim:
                    read.trim_read(barcode_name, args.min_score)
                if read_type == 'FASTA':
                    out_reads.write(read.get_fasta_lines(args.best))
                else:  # read_type == 'FASTQ'
                    out_reads.write(read.get_fastq_lines(args.best))


def random_mode(barcode_seq, count, output_all):
    """
    This function runs the program in random test mode, where the given barcode is aligned to
    random sequences of varying length.
    """
    import numpy as np

    if barcode_seq in BARCODES:
        barcode_seq = BARCODES[barcode_seq].sequence

    min_read_length = len(barcode_seq)
    max_read_length = 100000
    regression_interval = 1000

    header = ['Random sequence count']
    if output_all:
        header += ['Sequence length', 'Log log read length', 'Raw alignment score']
    header += ['Slope', 'Intercept']
    print('\t'.join(header), flush=True)

    log_log_read_lengths = []
    raw_scores = []

    for i in range(count):
        # Generate random sequences such that the log log of their lengths is a uniform
        # distribution.
        log_log_rand_seq_length = random.uniform(math.log(math.log(min_read_length, 10.0), 10.0),
                                                 math.log(math.log(max_read_length, 10.0), 10.0))
        rand_seq_length = int(round(10.0 ** (10.0 ** log_log_rand_seq_length)))
        rand_seq = get_random_sequence(rand_seq_length)

        cpp_result = barcode_alignment(rand_seq, barcode_seq, MATCH_SCORE, MISMATCH_SCORE,
                                       GAP_OPEN_SCORE, GAP_EXTEND_SCORE)
        raw_score = int(cpp_result.split(',')[4])

        output_line = [str(i+1)]
        if output_all:
            output_line += [str(rand_seq_length), str(log_log_rand_seq_length), str(raw_score)]

        log_log_read_lengths.append(log_log_rand_seq_length)
        raw_scores.append(raw_score)

        regression_step = (i+1) % regression_interval == 0
        if regression_step:
            regression = np.polyfit(np.array(log_log_read_lengths), np.array(raw_scores), 1)
            slope, intercept = str(regression[0]), str(regression[1])
        else:
            slope, intercept = '', ''
        output_line += [slope, intercept]

        if output_all or regression_step:
            print('\t'.join(output_line))

    regression = np.polyfit(np.array(log_log_read_lengths), np.array(raw_scores), 1)
    slope, intercept = regression[0], regression[1]

    print()
    print('slope =', slope)
    print('intercept =', intercept)
    print()
    full_equation_string = 'expected_score = ' + str(slope) + ' * log10(log10(sequence_length))'
    if intercept > 0.0:
        full_equation_string += ' + ' + str(intercept)
    elif intercept < 0.0:
        full_equation_string += ' - ' + str(-intercept)
    print('Full equation predicting alignment score of the barcode to random sequences:')
    print(full_equation_string)
    print()


def get_header(args, output_settings, reads_are_2d):
    header = []
    if output_settings.read_name:
        header += ['Read name']
    if output_settings.read_lengths:
        if reads_are_2d:
            header += ['2D read length', 'Template read length', 'Complement read length',
                       'Longest read length']
        else:
            header += ['Read length']
    barcode_names = args.barcodes.split(',')
    for barcode_name in barcode_names:
        if reads_are_2d:
            directions = ['2D, forward', '2D, reverse',
                          'template, forward', 'template, reverse',
                          'complement, forward', 'complement, reverse']
        else:
            directions = ['forward', 'reverse']
        for i, direction in enumerate(directions):
            if output_settings.all_scores:
                header.append(barcode_name + ' score (' + direction + ' barcode)')
            if output_settings.start_and_end:
                header.append(barcode_name + ' read start (' + direction + ' barcode)')
                header.append(barcode_name + ' read end (' + direction + ' barcode)')
            if output_settings.seq:
                header.append(barcode_name + ' read seq (' + direction + ' barcode)')

    if output_settings.best_scores:
        for barcode_name in barcode_names:
            header.append(barcode_name + ' read start score')
        header.append('Best start match')
        for barcode_name in barcode_names:
            header.append(barcode_name + ' read end score')
        header.append('Best end match')
        for barcode_name in barcode_names:
            header.append(barcode_name + ' mean score')
        header.append('Best mean match')

    if output_settings.call_confidence:
        header.append('Call confidence')

    if output_settings.final_call:
        header += ['Barcode call']

    return '\t'.join(header)


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


def convert_to_string(value):
    if isinstance(value, float):
        return '%.2f' % value
    else:
        return str(value)


class NanoporeRead(object):

    def __init__(self, name):
        self.name = name
        self.read_2d = None
        self.read_template = None
        self.read_complement = None
        self.barcode_alignment_results = {}
        self.best_start_match_barcode_name = ''
        self.best_end_match_barcode_name = ''
        self.best_mean_match_barcode_name = ''
        self.final_call = 'None'
        self.call_confidence = 0.0

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

    def get_read_lengths(self, reads_are_2d):
        if reads_are_2d:
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
        else:
            return len(self.read_template[1])

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
                trimmed_length = len(read_seq) - len(trimmed_read_seq)

                cpp_result = barcode_alignment(trimmed_read_seq, barcode.get_sequence(direction),
                                               MATCH_SCORE, MISMATCH_SCORE, GAP_OPEN_SCORE,
                                               GAP_EXTEND_SCORE)
                all_results.append(SingleBarcodeAlignmentResult(cpp_result, barcode, expected_end,
                                                                trimmed_read_seq, trimmed_length,
                                                                MATCH_SCORE))
        self.barcode_alignment_results[barcode.name] = OverallBarcodeAlignmentResult(all_results)

    def determine_best_barcode_match(self, barcodes, min_score, min_gap):
        start_scores = [(b.name, self.barcode_alignment_results[b.name].best_start_score)
                        for b in barcodes]
        start_scores = sorted(start_scores, key=lambda x: x[1], reverse=True)
        self.best_start_match_barcode_name = start_scores[0][0]
        best_start_score = start_scores[0][1]

        end_scores = [(b.name, self.barcode_alignment_results[b.name].best_end_score)
                      for b in barcodes]
        end_scores = sorted(end_scores, key=lambda x: x[1], reverse=True)
        self.best_end_match_barcode_name = end_scores[0][0]
        best_end_score = end_scores[0][1]

        scores = [(b.name, self.barcode_alignment_results[b.name].mean_best_score)
                  for b in barcodes]
        scores = sorted(scores, key=lambda x: x[1], reverse=True)
        self.best_mean_match_barcode_name = scores[0][0]
        best_score = scores[0][1]

        # If the best score is too low, we don't make a barcode call.
        if best_score < min_score:
            return

        # If the best matching barcode at the start and end don't agree and both have a decent
        # score, then we don't make a call.
        if self.best_start_match_barcode_name != self.best_end_match_barcode_name and \
                best_start_score >= min_score and best_end_score >= min_score:
            return

        self.call_confidence = 100 * (best_score - min_score) / (100.0 - min_score)

        # If there's only one possible barcode (an odd situation), then it gets the call.
        if len(scores) == 1:
            self.final_call = self.best_mean_match_barcode_name

        # If there are multiple possible barcodes, then the best must be sufficiently better than
        # the second best to be called.
        else:
            second_best_score = scores[1][1]
            gap = best_score - second_best_score
            if gap >= min_gap:
                self.final_call = self.best_mean_match_barcode_name
                self.call_confidence *= min(1.0, gap / 100.0)
            else:
                pass  # No barcode call

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

    def get_output_line(self, output_settings, barcodes, reads_are_2d):
        output_line = []
        if output_settings.read_name:
            output_line.append(self.name)
        if output_settings.read_lengths:
            output_line += self.get_read_lengths(reads_are_2d)
        for barcode in barcodes:
            results = self.barcode_alignment_results[barcode.name]
            output_line += results.get_output_line(output_settings, reads_are_2d)
        if output_settings.best_scores:
            for barcode in barcodes:
                output_line.append(self.barcode_alignment_results[barcode.name].best_start_score)
            output_line.append(self.best_start_match_barcode_name)
            for barcode in barcodes:
                output_line.append(self.barcode_alignment_results[barcode.name].best_end_score)
            output_line.append(self.best_end_match_barcode_name)
            for barcode in barcodes:
                output_line.append(self.barcode_alignment_results[barcode.name].mean_best_score)
            output_line.append(self.best_mean_match_barcode_name)
        if output_settings.call_confidence:
            output_line.append(self.call_confidence)
        if output_settings.final_call:
            output_line.append(self.final_call)
        return output_line

    @staticmethod
    def get_fasta_lines_one_type(read):
        if len(read[1]) == 0:
            return ''
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
        if len(read[1]) == 0:
            return ''
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
        if len(qualities) < 1:
            return 1.0

        error_count = 0.0
        for score in qualities:
            phred = ord(score) - 33
            error_count += 10.0 ** (-phred / 10.0)
        return error_count / len(qualities)

    def trim_read(self, barcode_name, min_score):
        """
        Removes the parts of the read from the barcode alignments onward.
        """
        if self.has_2d_sequence():
            self.trim_read_one_type('2d', barcode_name, min_score)
        if self.has_template_sequence():
            self.trim_read_one_type('template', barcode_name, min_score)
        if self.has_complement_sequence():
            self.trim_read_one_type('complement', barcode_name, min_score)

    def trim_read_one_type(self, read_type, barcode_name, min_score):
        """
        Removes the parts of the read from the barcode alignments onward.
        """
        if barcode_name == 'None':
            return

        results = self.barcode_alignment_results[barcode_name]

        if read_type == '2d':
            read = self.read_2d
            start_results = results.single_results[1]
            end_results = results.single_results[0]
        elif read_type == 'template':
            read = self.read_template
            start_results = results.single_results[3]
            end_results = results.single_results[2]
        else:  # read_type == 'complement'
            read = self.read_complement
            start_results = results.single_results[4]
            end_results = results.single_results[5]

        if start_results.adjusted_score >= min_score:
            start_trim = start_results.read_end_pos + 1
        else:
            start_trim = 0

        if end_results.adjusted_score >= min_score:
            end_trim = end_results.read_start_pos
        else:
            end_trim = len(read[1])

        if len(read) > 3:
            new_read = (read[0], read[1][start_trim:end_trim], read[2],
                        read[3][start_trim:end_trim], read[4])
        else:
            new_read = (read[0], read[1][start_trim:end_trim], read[2])

        if read_type == '2d':
            self.read_2d = new_read
        elif read_type == 'template':
            self.read_template = new_read
        else:  # read_type == 'complement'
            self.read_complement = new_read
