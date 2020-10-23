#!/usr/bin/env python3
# fastcov.py
# SK

import os
import sys
import time
import argparse
import numpy as np
import pysam
import ctypes
from multiprocessing import Pool
import matplotlib as mpl
mpl.use('Agg') # do not require X window
from matplotlib import pyplot as plt


#####

def error(string, error_type=1):
    sys.stderr.write(f'ERROR: {string}\n')
    sys.exit(error_type)


def log(string, newline_before=False):
    if newline_before:
        sys.stderr.write('\n')
    sys.stderr.write(f'LOG: {string}\n')

#####


def worker_get_coverage_chromo(bamfile, chromo, len_chromo):

    log(f'{os.getpid()} working on {bamfile} ...')
    start_time = time.time()
    with pysam.AlignmentFile(bamfile, 'rb') as pysam_bamfile:
        # fetch alignments
        alignments = pysam_bamfile.fetch(chromo, 0, len_chromo)
        num_aln = 0

        cov = np.zeros(len_chromo, dtype=int)

        for aln in alignments:
            if aln.is_secondary:
                continue
            num_aln += 1
            cov[aln.positions] += 1
        
    log(f'{bamfile}: {num_aln} alignments done in {time.time()-start_time:.2f} sec.')

    return num_aln, cov
    


if __name__ == '__main__':
    main_start_time = time.time()
    log('Started fastcov.py ...')


    # parse arguments
    bamfiles = sys.argv[1:]
    num_bamfiles = len(bamfiles)

    parser = argparse.ArgumentParser(description='Plot the coverage based on some bam files.')
    parser.add_argument("bamfile", nargs='+', help="Alignment files to include in the coverage plot.")
    args = parser.parse_args()


    # input checking
    bamfiles = args.bamfile
    for bamfile in bamfiles:
        if not os.path.isfile(bamfile):
            error(f'Not a file: {bamfile}')
        if not os.path.isfile(bamfile + '.bai'):
            log(f'Bam index missing for file: {bamfile}. Creating ...')
            ret = os.system(f'samtools index {bamfile}')
            if ret != 0:
                log(f'Warning: samtools index returned exit code {ret} - everything might crash and burn.')



    # find length of sequence to get coverage for
    chromo = None
    len_chromo = None
    
    with pysam.AlignmentFile(bamfiles[0], 'rb') as pysam_bamfile:
        chromo = pysam_bamfile.references[0]
        len_chromo = pysam_bamfile.lengths[0]

    log(f'Got reference name: {chromo} - length: {len_chromo}')


    # find out how many CPUs we can use
    num_cpus = 2
    max_cpus = os.cpu_count()
    if not max_cpus:
        log('Could not determine number of available CPUs, defaulting to 2.')
    else:
        if max_cpus > 2:
            num_cpus = max_cpus
        log(f'Using {num_cpus} of {max_cpus} available CPUs.')


    # make a numpy array to hold cov data
    cov_data = np.zeros([len(bamfiles), len_chromo], dtype=int)
    process_results = []


    # make worker pool
    pool = Pool(num_cpus)

    log(f'Made worker pool with {num_cpus} processes.')


    log(f'Parsing alignments ...')
    parse_start_time = time.time()

    # send work to pool
    for bamfile in bamfiles:
        res = pool.apply_async(worker_get_coverage_chromo, (bamfile, chromo, len_chromo))
        process_results.append(res)

    # wait for results
    total_num_aln = 0
    for i, res in enumerate(process_results):
        # write to array
        num_aln, cov_data[i] = res.get()
        total_num_aln += num_aln

    parse_end_time = time.time()
    parse_time = parse_end_time - parse_start_time
    log(f'Parsed {total_num_aln} alignments in {parse_time:.2f} seconds')


    print(cov_data)


    # plot stuff
    # TODO

