"""
Miscellaneous function used by various scripts.

Author: Ryan Wick
email: rrwick@gmail.com
"""

import sys
import os
import subprocess
import random
import math
import gzip


def float_to_str(num, decimals, max_num=0):
    """
    Converts a number to a string. Will add left padding based on the max value to ensure numbers
    align well.
    """
    if decimals == 0:
        return int_to_str(int(round(num)), max_num=max_num)
    if num is None:
        num_str = 'n/a'
    else:
        num_str = '%.' + str(decimals) + 'f'
        num_str = num_str % num
        parts = num_str.split('.')
        before_decimal = parts[0]
        after_decimal = parts[1]
        num_str = int_to_str(int(before_decimal)) + '.' + after_decimal
    if max_num > 0:
        max_str = float_to_str(max_num, decimals)
        num_str = num_str.rjust(len(max_str))
    return num_str


def int_to_str(num, max_num=0):
    """
    Converts a number to a string. Will add left padding based on the max value to ensure numbers
    align well.
    """
    if num is None:
        num_str = 'n/a'
    else:
        num_str = '{:,}'.format(num)
    max_str = '{:,}'.format(int(max_num))
    return num_str.rjust(len(max_str))


def check_files_and_programs(files, spades_path=None, graphmap_path=None,
                             makeblastdb_path=None, tblastn_path=None, gene_db_path=None,
                             pilon_path=None, samtools_path=None, bowtie2_path=None,
                             bowtie2_build_path=None):
    """
    Checks to make sure all files in the list are present and either program, as needed.
    """
    for file in files:
        check_file_exists(file)
    if graphmap_path:
        check_graphmap(graphmap_path)
    if spades_path:
        check_spades(spades_path)
    if makeblastdb_path and tblastn_path:
        check_blast(makeblastdb_path, tblastn_path, gene_db_path)
    if pilon_path:
        check_pilon(pilon_path, samtools_path, bowtie2_path, bowtie2_build_path)


def check_file_exists(filename):  # type: (str) -> bool
    """
    Checks to make sure the single given file exists.
    """
    if not os.path.isfile(filename):
        quit_with_error('could not find ' + filename)


def quit_with_error(message):  # type: (str) -> None
    """
    Displays the given message and ends the program's execution.
    """
    print('Error:', message, file=sys.stderr)
    sys.exit(1)


def check_graphmap(graphmap_path):
    """
    Makes sure the GraphMap executable is available.
    """
    if not find_program_with_which(graphmap_path):
        quit_with_error('could not find GraphMap at ' + graphmap_path +
                        ', either fix path or run with --no_graphmap')


def check_spades(spades_path):
    """
    Makes sure the SPAdes executable is available.
    """
    if not find_program_with_which(spades_path):
        quit_with_error('could not find SPAdes at ' + spades_path)

    command = [spades_path, '-h']
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = process.communicate()

    if not err.decode():
        quit_with_error('SPAdes was found but does not produce output (make sure to use '
                        '"spades.py" location, not "spades")')


def check_blast(makeblastdb_path, tblastn_path, gene_db_path):
    """
    Makes sure the BLAST executables are available.
    """
    if not find_program_with_which(makeblastdb_path):
        quit_with_error('could not find makeblastdb - either specify its location using '
                        '--makeblastdb_path or use --no_rotate to remove BLAST dependency')
    if not find_program_with_which(tblastn_path):
        quit_with_error('could not find tblastn - either specify its location using '
                        '--tblastn_path or use --no_rotate to remove BLAST dependency')
    if not os.path.isfile(gene_db_path):
        quit_with_error('could not find file: ' + gene_db_path +
                        '\neither specify a different start gene database using --start_genes '
                        'or use --no_rotate')


def check_pilon(pilon_path, samtools_path, bowtie2_path, bowtie2_build_path):
    """
    Makes sure the Pilon executable is available.
    """
    if not find_program_with_which('java'):
        quit_with_error('could not find java - either specify its location using '
                        '--pilon_path or use --no_pilon to remove Java dependency')
    if not find_program_with_which(samtools_path):
        quit_with_error('could not find samtools - either specify its location using '
                        '--samtools_path or use --no_pilon to remove Samtools dependency')
    if not find_program_with_which(bowtie2_path):
        quit_with_error('could not find bowtie2 - either specify its location using '
                        '--bowtie2_path or use --no_pilon to remove Bowtie2 dependency')
    if not find_program_with_which(bowtie2_build_path):
        quit_with_error('could not find bowtie2-build - either specify its location using '
                        '--bowtie2_build_path or use --no_pilon to remove Bowtie2 dependency')
    if not get_pilon_jar_path(pilon_path):
        quit_with_error('could not find pilon.jar - either specify its location using --pilon_path '
                        'or use --no_pilon to remove Pilon dependency')


def get_pilon_jar_path(pilon_path):
    """
    Returns the path to pilon.jar. If the given path is correct, it just returns that,
    as an absolute path. Otherwise it tries to use the which command to find it.
    """
    if os.path.isfile(pilon_path):
        return os.path.abspath(pilon_path)
    try:
        pilon_path = subprocess.check_output(['which', 'pilon.jar'],
                                             stderr=subprocess.STDOUT).decode().strip()
    except subprocess.CalledProcessError:
        pass
    if os.path.isfile(pilon_path):
        return pilon_path
    else:
        return None


def find_program_with_which(executable_path):
    """
    Returns whether or not the executable was found.
    """
    process = subprocess.Popen(['which', executable_path], stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE)
    out, err = process.communicate()
    return bool(out) and not bool(err)


def get_mean_and_st_dev(num_list):
    """
    This function returns the mean and standard deviation of the given list of numbers.
    """
    num = len(num_list)
    if num == 0:
        return None, None
    mean = sum(num_list) / num
    if num == 1:
        return mean, None
    sum_squares = sum((x - mean) ** 2 for x in num_list)
    st_dev = (sum_squares / (num - 1)) ** 0.5
    return mean, st_dev


def print_progress_line(completed, total, base_pairs=None, prefix=None, end_newline=False):
    """
    Prints a progress line to the screen using a carriage return to overwrite the previous progress
    line.
    """
    progress_str = ''
    if prefix:
        progress_str += prefix
    progress_str += int_to_str(completed) + ' / ' + int_to_str(total)
    if total > 0:
        percent = 100.0 * completed / total
    else:
        percent = 0.0
    progress_str += ' (' + '%.1f' % percent + '%)'
    if base_pairs is not None:
        progress_str += ' - ' + int_to_str(base_pairs) + ' bp'
    end_char = '\n' if end_newline else ''
    print('\r' + progress_str, end=end_char, flush=True)


def get_nice_header(header):
    """
    For a header with a SPAdes/Velvet format, this function returns a simplified string that is
    just NODE_XX where XX is the contig number.
    For any other format, this function trims off everything following the first whitespace.
    """
    if is_header_spades_format(header):
        return 'NODE_' + header.split('_')[1]
    else:
        return header.split()[0]


def is_header_spades_format(contig_name):
    """
    Returns whether or not the header appears to be in the SPAdes/Velvet format.
    Example: NODE_5_length_150905_cov_4.42519
    """
    contig_name_parts = contig_name.split('_')
    return len(contig_name_parts) > 5 and \
        (contig_name_parts[0] == 'NODE' or contig_name_parts[0] == 'EDGE') and \
        contig_name_parts[2] == 'length' and contig_name_parts[4] == 'cov'


def reverse_complement(seq):
    """
    Given a DNA sequences, this function returns the reverse complement sequence.
    """
    return ''.join([complement_base(seq[i]) for i in range(len(seq) - 1, -1, -1)])


def complement_base(base):
    """
    Given a DNA base, this returns the complement.
    """
    if base == 'A':
        return 'T'
    if base == 'T':
        return 'A'
    if base == 'G':
        return 'C'
    if base == 'C':
        return 'G'
    if base == 'a':
        return 't'
    if base == 't':
        return 'a'
    if base == 'g':
        return 'c'
    if base == 'c':
        return 'g'
    forward = 'RYSWKMryswkmBDHVbdhvNn.-?'
    reverse = 'YRSWMKyrswmkVHDBvhdbNn.-?N'
    return reverse[forward.find(base)]


def get_random_base():
    """
    Returns a random base with 25% probability of each.
    """
    rand_int = random.randint(0, 3)
    if rand_int == 0:
        return 'A'
    elif rand_int == 1:
        return 'C'
    elif rand_int == 2:
        return 'G'
    elif rand_int == 3:
        return 'T'


def get_random_sequence(length):
    """
    Returns a random sequence of the given length.
    """
    sequence = ''
    for _ in range(length):
        sequence += get_random_base()
    return sequence


def get_median(sorted_list):
    """
    Returns the median of a list of numbers. Assumes the list has already been sorted.
    """
    count = len(sorted_list)
    index = (count - 1) // 2
    if count % 2:
        return sorted_list[index]
    else:
        return (sorted_list[index] + sorted_list[index + 1]) / 2.0


def get_percentile(unsorted_list, percentile):
    """
    Returns a percentile of a list of numbers. Doesn't assume the list has already been sorted.
    Implements the nearest rank method:
    https://en.wikipedia.org/wiki/Percentile#The_Nearest_Rank_method
    """
    return get_percentile_sorted(sorted(unsorted_list), percentile)


def get_percentile_sorted(sorted_list, percentile):
    """
    Same as the above function, but assumes the list is already sorted.
    """
    if not sorted_list:
        return 0.0
    fraction = percentile / 100.0
    rank = int(math.ceil(fraction * len(sorted_list)))
    if rank == 0:
        return sorted_list[0]
    return sorted_list[rank - 1]


def weighted_average(num_1, num_2, weight_1, weight_2):
    """
    A simple weighted mean of two numbers.
    """
    weight_sum = weight_1 + weight_2
    return num_1 * (weight_1 / weight_sum) + num_2 * (weight_2 / weight_sum)


def weighted_average_list(nums, weights):
    """
    A simple weighted mean of a list of numbers.
    """
    w_sum = sum(weights)
    if w_sum == 0.0:
        return 0.0
    else:
        return sum(num * (weights[i] / w_sum) for i, num in enumerate(nums))


def print_section_header(message, verbosity, last_newline=True):
    """
    Prints a header for std out, unless verbosity is zero, in which case it does nothing.
    """
    if verbosity > 0:
        print('\n')
        print(message)
        end_char = '\n' if last_newline else ''
        print('-' * len(message), flush=True, end=end_char)


def round_to_nearest_odd(num):
    """
    Rounds a float to an odd integer.
    """
    round_up = int(math.ceil(num))
    if round_up % 2 == 0:
        round_up += 1
    round_down = int(math.floor(num))
    if round_down % 2 == 0:
        round_down -= 1
    up_diff = round_up - num
    down_diff = num - round_down
    if up_diff > down_diff:
        return round_down
    else:
        return round_up


def get_compression_type(filename):
    """
    Attempts to guess the compression (if any) on a file using the first few bytes.
    http://stackoverflow.com/questions/13044562
    """
    magic_dict = {'gz': (b'\x1f', b'\x8b', b'\x08'),
                  'bz2': (b'\x42', b'\x5a', b'\x68'),
                  'zip': (b'\x50', b'\x4b', b'\x03', b'\x04')}
    max_len = max(len(x) for x in magic_dict)

    unknown_file = open(filename, 'rb')
    file_start = unknown_file.read(max_len)
    unknown_file.close()
    compression_type = 'plain'
    for filetype, magic_bytes in magic_dict.items():
        if file_start.startswith(magic_bytes):
            compression_type = filetype
    if compression_type == 'bz2':
        quit_with_error('cannot use bzip2 format - use gzip instead')
    if compression_type == 'zip':
        quit_with_error('cannot use zip format - use gzip instead')
    return compression_type


def get_sequence_file_type(filename):
    """
    Determines whether a file is FASTA or FASTQ.
    """
    if get_compression_type(filename) == 'gz':
        open_func = gzip.open
    else:  # plain text
        open_func = open
    seq_file = open_func(filename, 'rt')
    first_char = seq_file.read(1)

    if first_char == '>':
        return 'FASTA'
    elif first_char == '@':
        return 'FASTQ'
    else:
        raise ValueError('File is neither FASTA or FASTQ')


def get_num_agreement(num_1, num_2):
    """
    Returns a value between 0.0 and 1.0 describing how well the numbers agree.
    1.0 is perfect agreement and 0.0 is the worst.
    """
    if num_1 == 0.0 and num_2 == 0.0:
        return 1.0
    if num_1 < 0.0 and num_2 < 0.0:
        num_1 *= -1
        num_2 *= -1
    if num_1 * num_2 < 0.0:
        return 0.0
    return min(num_1, num_2) / max(num_1, num_2)


def flip_number_order(num_1, num_2):
    """
    Given two segment numbers, this function possibly flips them around. It returns the new numbers
    (either unchanged or flipped) and whether or not a flip took place. The decision is somewhat
    arbitrary, but it needs to be consistent so when we collect bridging read sequences they are
    always in the same direction.
    """
    if num_1 > 0 and num_2 > 0:
        flip = False
    elif num_1 < 0 and num_2 < 0:
        flip = True
    elif num_1 < 0:  # only num_1 is negative
        flip = abs(num_1) > abs(num_2)
    else:  # only num_2 is negative
        flip = abs(num_2) > abs(num_1)
    if flip:
        return (-num_2, -num_1), True
    else:
        return (num_1, num_2), False


def load_fasta_or_fastq(filename):
    """
    Returns a list of tuples (header, seq) for each record in the fasta/fastq file.
    """
    file_type = get_sequence_file_type(filename)
    if file_type == 'FASTA':
        return load_fasta(filename)
    else:  # FASTQ
        return load_fastq(filename)


def load_fasta(filename):
    """
    Returns a list of tuples (header, seq) for each record in the fasta file.
    """
    fasta_seqs = []
    fasta_file = open(filename, 'rt')
    name = ''
    sequence = ''
    for line in fasta_file:
        line = line.strip()
        if not line:
            continue
        if line[0] == '>':  # Header line = start of new contig
            if name:
                fasta_seqs.append((name.split()[0], sequence))
                sequence = ''
            name = line[1:]
        else:
            sequence += line
    if name:
        fasta_seqs.append((name.split()[0], sequence))
    fasta_file.close()
    return fasta_seqs


def load_fastq(fastq_filename):
    """
    Returns a list of tuples (header, seq) for each record in the fastq file.
    """
    if get_compression_type(fastq_filename) == 'gz':
        open_func = gzip.open
    else:  # plain text
        open_func = open
    reads = []
    with open_func(fastq_filename, 'rt') as fastq:
        for line in fastq:
            name = line.strip()[1:].split()[0]
            sequence = next(fastq).strip()
            _ = next(fastq)
            _ = next(fastq)
            reads.append((name, sequence))
    return reads


def score_function(val, half_score_val):
    """
    For inputs of 0.0 and greater, this function returns a value between 0.0 and 1.0, approaching
    1.0 with large values. The half_score_val argument is the point at which the function returns
    0.5. If it's large the function approaches 1.0 more slowly, if it's small the function
    approaches 1.0 more quickly.
    """
    return 1.0 - (half_score_val / (half_score_val + val))


def strip_read_extensions(read_file_name):
    """
    This function removes extensions from a file name.
    """
    base_name = os.path.basename(read_file_name)
    name_parts = base_name.split('.')
    for i in range(2):
        if len(name_parts) > 1 and len(name_parts[-1]) <= 5:
            name_parts = name_parts[:-1]
    return '.'.join(name_parts)
