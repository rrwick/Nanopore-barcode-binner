#!/usr/bin/env python3
"""
Nanopore barcoder

Author: Ryan Wick
email: rrwick@gmail.com
"""

import argparse
import os
import sys
import random
import math
from multiprocessing import cpu_count
import statistics
from .cpp_function_wrappers import barcode_alignment
from .misc import load_fasta_or_fastq, reverse_complement
from .version import __version__

BARCODES = {'NB01': {'name': 'NB01', 'sequence': 'AAGAAAGTTGTCGGTGTCTTTGTG', 'slope': 71.4332584271663, 'intercept': -22.5656089159852},
            'NB02': {'name': 'NB02', 'sequence': 'TCGATTCCGTTTGTAGTCGTCTGT', 'slope': 67.1855037419996, 'intercept': -18.0468313716098},
            'NB03': {'name': 'NB03', 'sequence': 'GAGTCTTGTGTCCCAGTTACCAGG', 'slope': 63.5927250003896, 'intercept': -15.9980232349283},
            'NB04': {'name': 'NB04', 'sequence': 'TTCGGATTCTATCGTGTTTCCCTA', 'slope': 66.8473788581195, 'intercept': -18.5724631096649},
            'NB05': {'name': 'NB05', 'sequence': 'CTTGTCCAGGGTTTGTGTAACCTT', 'slope': 66.5864343955349, 'intercept': -18.8843342040538},
            'NB06': {'name': 'NB06', 'sequence': 'TTCTCGCAAAGGCAGAAAGTAGTC', 'slope': 65.4698361318286, 'intercept': -17.6317088824286},
            'NB07': {'name': 'NB07', 'sequence': 'GTGTTACCGTGGGAATGAATCCTT', 'slope': 64.9387340902843, 'intercept': -17.4358282814267},
            'NB08': {'name': 'NB08', 'sequence': 'TTCAGGGAACAAACCAAGTTACGT', 'slope': 65.5252104348060, 'intercept': -18.0882967050121},
            'NB09': {'name': 'NB09', 'sequence': 'AACTAGGCACAGCGAGTCTTGGTT', 'slope': 63.9800483064580, 'intercept': -15.9533793977375},
            'NB10': {'name': 'NB10', 'sequence': 'AAGCGTTGAAACCTTTGTCCTCTC', 'slope': 66.5716250549626, 'intercept': -18.8997212170206},
            'NB11': {'name': 'NB11', 'sequence': 'GTTTCATCTATCGGAGGGAATGGA', 'slope': 66.8540339246670, 'intercept': -18.6488190718658},
            'NB12': {'name': 'NB12', 'sequence': 'CAGGTAGAAAGAAGCAGAATCGGA', 'slope': 70.1581047865076, 'intercept': -20.4691636410478}}

MATCH_SCORE = 3
MISMATCH_SCORE = -6
GAP_OPEN_SCORE = -5
GAP_EXTEND_SCORE = -2

INCLUDE_RAW = False
INCLUDE_START_AND_END = False
INCLUDE_SEQ = False
INCLUDE_BEFORE_AND_AFTER_SEQ = True


def get_arguments():
    """
    Parse the command line arguments.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', '--reads', required=True,
                        help='FASTA or FASTQ of reads')
    parser.add_argument('-b', '--barcodes', default='NB01,NB02',
                        help='Names of barcodes, separated by commas')
    parser.add_argument('-d', '--distance', type=int, default=100,
                        help='Allowed distance between barcode and expected end of read')
    parser.add_argument('--min_score', type=float, default=20.0,
                        help='Reads must have at least this score for their best barcode to be classified')
    parser.add_argument('--gap', type=float, default=10.0,
                        help="A read's best score must be better than its second best score by this much to be classified")
    return parser.parse_args()


def main():
    """
    Script execution starts here.
    """
    for barcode_name in BARCODES:
        barcode_seq = BARCODES[barcode_name]['sequence']
        BARCODES[barcode_name]['rev_comp_sequence'] = reverse_complement(barcode_seq)

    args = get_arguments()
    barcodes = [BARCODES[x] for x in args.barcodes.split(',')]

    # Load the reads and group them up by their type (2D/template/complement).
    reads = load_fasta_or_fastq(args.reads)
    grouped_reads = {}
    read_names = []
    for read_name, read_seq in reads:
        if read_name.endswith('_Basecall_2D_2d'):
            short_name = read_name[:-15]
            read_type = '2d'
        elif read_name.endswith('_Basecall_2D_template'):
            short_name = read_name[:-21]
            read_type = 'template'
        elif read_name.endswith('_Basecall_2D_complement'):
            short_name = read_name[:-23]
            read_type = 'complement'
        else:  # Assume 1D
            short_name = read_name
            read_type = 'template'
        if short_name not in grouped_reads:
            read_names.append(short_name)
            grouped_reads[short_name] = {}
        grouped_reads[short_name][read_type] = read_seq

    header_string, \
       score_columns_forward_2D, score_columns_reverse_2D, \
       score_columns_forward_template, score_columns_reverse_template, \
       score_columns_forward_complement, score_columns_reverse_complement = get_header(args)

    print(header_string, flush=True)

    for read_name in read_names:
        output_line = [read_name]
        read_group = grouped_reads[read_name]

        if '2d' in read_group:
            read_seq_2d = read_group['2d']
            output_line.append(str(len(read_seq_2d)))
        else:
            read_seq_2d = None
            output_line.append('')

        if 'template' in read_group:
            read_seq_template = read_group['template']
            output_line.append(str(len(read_seq_template)))
        else:
            read_seq_template = None
            output_line.append('')

        if 'complement' in read_group:
            read_seq_complement = read_group['complement']
            output_line.append(str(len(read_seq_complement)))
        else:
            read_seq_complement = None
            output_line.append('')

        empty_result_count = 1
        if INCLUDE_RAW:
            empty_result_count += 2
        if INCLUDE_START_AND_END:
            empty_result_count += 2
        if INCLUDE_BEFORE_AND_AFTER_SEQ:
            empty_result_count += 2
        if INCLUDE_SEQ:
            empty_result_count += 1
        empty_results = [''] * empty_result_count

        all_barcode_scores = {}
        for i, barcode in enumerate(barcodes):
            perfect_score = len(barcode['sequence']) * MATCH_SCORE
            output_line.append(str(perfect_score))

            if '2d' in read_group:
                output_line += align_barcode(read_seq_2d, barcode['sequence'], barcode['slope'], barcode['intercept'], perfect_score, 'END', args.distance)
            else:
                output_line += empty_results
            if '2d' in read_group:
                output_line += align_barcode(read_seq_2d, barcode['rev_comp_sequence'], barcode['slope'], barcode['intercept'], perfect_score, 'START', args.distance)
            else:
                output_line += empty_results

            if 'template' in read_group:
                output_line += align_barcode(read_seq_template, barcode['sequence'], barcode['slope'], barcode['intercept'], perfect_score, 'END', args.distance)
            else:
                output_line += empty_results
            if 'template' in read_group:
                output_line += align_barcode(read_seq_template, barcode['rev_comp_sequence'], barcode['slope'], barcode['intercept'], perfect_score, 'START', args.distance)
            else:
                output_line += empty_results

            if 'complement' in read_group:
                output_line += align_barcode(read_seq_complement, barcode['sequence'], barcode['slope'], barcode['intercept'], perfect_score, 'END', args.distance)
            else:
                output_line += empty_results
            if 'complement' in read_group:
                output_line += align_barcode(read_seq_complement, barcode['rev_comp_sequence'], barcode['slope'], barcode['intercept'], perfect_score, 'START', args.distance)
            else:
                output_line += empty_results

            best_start_read_barcode_score = max_with_empty_strings(output_line[score_columns_reverse_2D[i]],
                                                                   output_line[score_columns_reverse_template[i]],
                                                                   output_line[score_columns_forward_complement[i]])

            best_end_read_barcode_score = max_with_empty_strings(output_line[score_columns_forward_2D[i]],
                                                                 output_line[score_columns_forward_template[i]],
                                                                 output_line[score_columns_reverse_complement[i]])
            mean_best_barcode_score = (best_start_read_barcode_score + best_end_read_barcode_score) / 2.0

            output_line += [best_start_read_barcode_score, best_end_read_barcode_score, mean_best_barcode_score]
            all_barcode_scores[barcode['name']] = mean_best_barcode_score

        # Find the best matching barcode
        all_barcode_scores = sorted(all_barcode_scores.items(), key=lambda x: x[1], reverse=True)
        best_name = all_barcode_scores[0][0]
        best_score = all_barcode_scores[0][1]
        output_line.append(best_name)  # Best match

        # If the best score is too low, we don't make a barcode call.
        if best_score < args.min_score:
            final_call = 'None'
            confidence = 0.0

        else:
            confidence = 100 * (best_score - args.min_score) / (100.0 - args.min_score)

            # If there's only one possible barcode (kind of a dumb situation), then it gets the call.
            if len(all_barcode_scores) == 1:
                final_call = best_name

            # If there are multiple possible barcodes, then the best must be sufficiently better than the second best.
            else:  # More than one barcode
                second_best_score = all_barcode_scores[1][1]
                gap = best_score - second_best_score
                if gap >= args.gap:
                    final_call = best_name
                    confidence *= min(1.0, gap / 100.0)
                else:
                    final_call = 'None'
                    confidence = 0.0

        output_line += [final_call, confidence]

        print('\t'.join([str(x) for x in output_line]), flush=True)


def align_barcode(read_seq, barcode_seq, slope, intercept, perfect_score, expected_end, allowed_distance):

    # Trim the read sequence to the region worth examining.
    full_read_length = len(read_seq)
    region_size = min(full_read_length, allowed_distance + len(barcode_seq) + 20)
    if expected_end == 'START':
        read_seq = read_seq[:region_size]
    else:  # expected_end == 'END'
        read_seq = read_seq[-region_size:]

    alignment_result = barcode_alignment(read_seq, barcode_seq, MATCH_SCORE, MISMATCH_SCORE, GAP_OPEN_SCORE, GAP_EXTEND_SCORE)
    parts = alignment_result.split(',')

    read_start_pos = int(parts[0])
    read_end_pos = int(parts[1])
    before_hit_sequence = read_seq[read_start_pos-20:read_start_pos]
    hit_sequence = read_seq[read_start_pos:read_end_pos+1]
    after_hit_sequence = read_seq[read_end_pos+1:read_end_pos+21]
    if expected_end == 'END':
        missing_start = full_read_length - region_size
        read_start_pos += missing_start
        read_end_pos += missing_start

    raw_score = int(parts[4])

    expected_score_for_random_seq = slope * math.log10(math.log10(full_read_length)) + intercept
    adjusted_score = 100.0 * (raw_score - expected_score_for_random_seq) / (perfect_score - expected_score_for_random_seq)

    output_line = []
    if INCLUDE_RAW:
        output_line.append(raw_score)
        output_line.append(expected_score_for_random_seq)
    output_line.append(adjusted_score)
    if INCLUDE_START_AND_END:
        output_line.append(read_start_pos)
        output_line.append(read_end_pos)
    if INCLUDE_BEFORE_AND_AFTER_SEQ:
        output_line.append(before_hit_sequence)
    if INCLUDE_SEQ:
        output_line.append(hit_sequence)
    if INCLUDE_BEFORE_AND_AFTER_SEQ:
        output_line.append(after_hit_sequence)

    return output_line  



def get_header(args):
    score_columns_forward_2D = []
    score_columns_reverse_2D = []
    score_columns_forward_template = []
    score_columns_reverse_template = []
    score_columns_forward_complement = []
    score_columns_reverse_complement = []
    all_scale_columns = [score_columns_forward_2D, score_columns_reverse_2D,
                         score_columns_forward_template, score_columns_reverse_template,
                         score_columns_forward_complement, score_columns_reverse_complement]

    header = ['Read name', '2D read length', 'template read length', 'complement read length']
    for barcode_name in args.barcodes.split(','):
        header.append(barcode_name + ' perfect score')

        directions = ['2D, forward', '2D, reverse', 'template, forward', 'template, reverse', 'complement, forward', 'complement, reverse']

        for i, direction in enumerate(directions):

            if INCLUDE_RAW:
                header.append(barcode_name + ' raw score (' + direction + ' barcode)')
                header.append(barcode_name + ' random score (' + direction + ' barcode)')

            header.append(barcode_name + ' adjusted score (' + direction + ' barcode)')
            all_scale_columns[i].append(len(header) - 1)

            if INCLUDE_START_AND_END:
                header.append(barcode_name + ' read start (' + direction + ' barcode)')
                header.append(barcode_name + ' read end (' + direction + ' barcode)')

            if INCLUDE_BEFORE_AND_AFTER_SEQ:
                header.append(barcode_name + ' before seq (' + direction + ' barcode)')
            if INCLUDE_SEQ:
                header.append(barcode_name + ' read seq (' + direction + ' barcode)')
            if INCLUDE_BEFORE_AND_AFTER_SEQ:
                header.append(barcode_name + ' after seq (' + direction + ' barcode)')

        header.append(barcode_name + ' best start read barcode score')
        header.append(barcode_name + ' best end read barcode score')
        header.append(barcode_name + ' mean best barcode score')

    header += ['Best barcode match', 'Final barcode call', 'Call confidence']
    header_string = '\t'.join(header)

    return header_string, \
           score_columns_forward_2D, score_columns_reverse_2D, \
           score_columns_forward_template, score_columns_reverse_template, \
           score_columns_forward_complement, score_columns_reverse_complement


def max_with_empty_strings(*args):
    return max(x for x in args if x != '')


def mean_with_empty_strings(*args):
    return statistics.mean(x for x in args if x != '')



