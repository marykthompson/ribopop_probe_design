'''
probe_helpers.py
Functions used for parsing ranges, etc.
'''

from collections.abc import Iterable
import pandas as pd
from collections import defaultdict
import numpy as np

def range_defined(ranges):
    '''
    Check whether ranges is defined. Let's us differentiate between nan/None and
    a potential list of ranges.
    '''
    if isinstance(ranges, Iterable):
        return(not pd.isnull(ranges).all())
    else:
        return(not pd.isnull(ranges))

def get_subregion_ranges(ranges_string):
    '''
    Convert the provided ranges string into an array, 0-based, closed interval.
    The reason for making it closed is to match the alignment-derived intervals.
    e.g. "1-2000,2300-3400" -> [[0, 1999], [2299, 3399]]
    '''
    subregion_ranges = []
    subregion_list = [i for i in ranges_string.replace(', ',',').split(',')]
    for i in subregion_list:
        start, end = [int(j.strip()) for j in i.split('-')]
        #convert to 0-based, half-open
        start -= 1
        end -= 1
        subregion_ranges.append([start, end])
    subregion_ranges = np.array(subregion_ranges)
    return subregion_ranges

def set_snake_args(args, snakemake):
    '''
    A hack to get the argparse defaults overriden by snakemake args, if provided
    '''
    argdict = vars(args)
    toset = defaultdict(list)
    for arg in argdict:
        if hasattr(snakemake.input, arg):
            toset['input'].append(arg)
        if hasattr(snakemake.params, arg):
            toset['params'].append(arg)
        if hasattr(snakemake.output, arg):
            toset['output'].append(arg)

    for block in toset:
        for arg in toset[block]:
            setattr(args, arg, getattr(getattr(snakemake, block), arg))
    return args
