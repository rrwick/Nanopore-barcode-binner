import random
import math

max_read_length = 1000000.0

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


with open('random_reads2.fasta', 'wt') as random_reads:
    for i in range(5000):
        log_read_length = random.uniform(math.log(24.0, 10.0), math.log(max_read_length, 10.0))
        read_length = int(round(10 ** log_read_length))
        rand_seq = get_random_sequence(read_length)
        random_reads.write('>read')
        random_reads.write(str(i+1))
        random_reads.write('\n')
        random_reads.write(rand_seq)
        random_reads.write('\n')
