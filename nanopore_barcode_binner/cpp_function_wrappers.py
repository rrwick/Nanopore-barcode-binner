"""
This script makes use of several C++ functions which are in cpp_functions.so. They are wrapped in
similarly named Python functions.

Author: Ryan Wick
email: rrwick@gmail.com
"""

import os
from ctypes import CDLL, cast, c_char_p, c_int, c_double, c_void_p, c_bool, POINTER
from .misc import quit_with_error

SO_FILE = 'cpp_functions.so'
SO_FILE_FULL = os.path.join(os.path.dirname(os.path.realpath(__file__)), SO_FILE)
if not os.path.isfile(SO_FILE_FULL):
    quit_with_error('could not find ' + SO_FILE + ' - please reinstall')
C_LIB = CDLL(SO_FILE_FULL)

C_LIB.barcodeAlignment.argtypes = [c_char_p,  # Read sequence
                                   c_char_p,  # Barcode sequence
                                   c_int,     # Match score
                                   c_int,     # Mismatch score
                                   c_int,     # Gap open score
                                   c_int]     # Gap extension score
C_LIB.barcodeAlignment.restype = c_void_p     # String describing alignment


# This function cleans up the heap memory for the C strings returned by the other C functions. It
# must be called after them.
C_LIB.freeCString.argtypes = [c_void_p]
C_LIB.freeCString.restype = None


def barcode_alignment(read_sequence, barcode_sequence,
                      match_score, mismatch_score, gap_open_score, gap_extend_score):
    """
    Python wrapper for barcodeAlignment C++ function.
    """
    ptr = C_LIB.barcodeAlignment(read_sequence.encode('utf-8'), barcode_sequence.encode('utf-8'),
                                 match_score, mismatch_score, gap_open_score, gap_extend_score)
    result_string = c_string_to_python_string(ptr)

    return result_string


def c_string_to_python_string(c_string):
    """
    This function casts a C string to a Python string and then calls a function to delete the C
    string from the heap.
    """
    python_string = cast(c_string, c_char_p).value.decode()
    C_LIB.freeCString(c_string)
    return python_string
