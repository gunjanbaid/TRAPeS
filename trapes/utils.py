"""
This module contains utility functions for various tasks.

"""

import os
import subprocess


def format_files(full_path, bam, unmapped, output):
    found = True
    nbam = full_path + bam
    if bam.startswith('/'):
        nbam = full_path + bam[1:]
    elif bam.startswith('./'):
        nbam = full_path + bam[2:]
    nunmapped = full_path + unmapped
    if unmapped.startswith('/'):
        nunmapped = full_path + unmapped[1:]
    if unmapped.startswith('./'):
        nunmapped = full_path + unmapped[2:]
    if ((os.path.isfile(nunmapped)) & (os.path.isfile(nbam))):
        noutput = make_output_dir(output, full_path)
    else:
        noutput = output
        found = False
    return (found, nbam, nunmapped, noutput)


def make_output_dir(output, full_path):
    noutput = output
    if output.startswith('/'):
        noutput = output[1:]
    if output.startswith('./'):
        noutput = output[2:]
    if output.endswith('/'):
        noutput = output[:-1]
    if output.find('/') != -1:
        out_arr = noutput.split('/')
        curr_path = full_path
        for i in range(0, len(out_arr) - 1):
            curr_path = curr_path + out_arr[i] + '/'
            if not os.path.exists(curr_path):
                os.makedirs(curr_path)
    noutput = full_path + noutput
    return noutput


def find_seq_and_length_of_aa(aa_seq):
    f_len = 0
    f_seq = ''
    start_arr = []
    stop_arr = []
    start_m = aa_seq.find('M')
    while start_m != -1:
        start_arr.append(start_m)
        start_m = aa_seq.find('M', start_m + 1)
    stop_pos = aa_seq.find('*')
    while stop_pos != -1:
        stop_arr.append(stop_pos)
        stop_pos = aa_seq.find('*', stop_pos + 1)

    if ((len(start_arr) == 0) | (len(stop_arr) == 0)):
        return (f_seq, f_len)
    for st_p in start_arr:
        curr_stop = find_stop(st_p, stop_arr)
        if curr_stop == -1:
            return (f_seq, f_len)
        else:
            curr_len = curr_stop - st_p
            if curr_len >= f_len:
                f_len = curr_len
                f_seq = aa_seq[st_p:curr_stop]
    return f_seq


def find_stop(st_p, stop_arr):
    for x in stop_arr:
        if x > st_p:
            return x
    return -1


def get_c_info(bed_entry, id_name_dict, fasta_dict):
    bed_arr = bed_entry.strip('\n').split('\t')
    c_id = bed_arr[3]
    c_name = id_name_dict[c_id]
    c_seq = fasta_dict[c_id]
    return (c_seq, c_name, c_id)


def is_paired_end(filename):
    """
    Parameters
    ----------
    filename : str
        Name of BAM/SAM file.

    Returns
    -------
    boolean
        True if reads are paired-end, False otherwise. 

    """
    return int(subprocess.check_output(["samtools", "view", "-c", "-f", "1", filename], universal_newlines=True).strip())


def remove_file(*args):
    for filename in args:
        subprocess.call(["rm", filename])

