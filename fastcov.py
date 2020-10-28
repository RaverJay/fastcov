#!/usr/bin/env python3
# fastcov.py
# SK

import os
import sys
import time
import argparse
import numpy as np
import pandas as pd
import seaborn as sns
import pysam
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

def worker_get_coverage_ref_name(bamfile, ref_name, pos_start, pos_end):

    log(f'{os.getpid()} working on {bamfile} ...')
    start_time = time.time()
    with pysam.AlignmentFile(bamfile, 'rb') as pysam_bamfile:

        # fetch alignments
        alignments = pysam_bamfile.fetch(ref_name, pos_start, pos_start + pos_len)
        num_aln = 0

        coverage = np.zeros(pos_len, dtype=int)

        for aln in alignments:
            num_aln += 1
            if aln.is_secondary:
                continue

            # record coverage
            for pos in aln.positions:
                pos_shift = pos - pos_start
                # early stop
                if pos_shift >= pos_len:
                    break             

                if 0 <= pos_shift < pos_len:
                    coverage[pos_shift] += 1
        
    log(f'{bamfile}: {num_aln} alignments done in {time.time()-start_time:.2f} sec.')

    return num_aln, coverage
    

if __name__ == '__main__':
    main_start_time = time.time()
    log('Started fastcov.py ...')


    ###
    # parse arguments
    parser = argparse.ArgumentParser(description='Plot the coverage based on some bam files.')
    parser.add_argument("bamfile", nargs='+', help="Alignment files to include in the coverage plot.")
    parser.add_argument("-p", "--position", help="Specify a genomic position to plot exclusively. Format: <ref_name>[:<start>-<stop>] (i.e. coordinates are optional and must be 1-based and inclusive)")
    parser.add_argument("-l", "--logscale", action='store_true', help="Use logarithmic scale on y-axis.")
    args = parser.parse_args()


    ###
    # input checking
    bamfiles = args.bamfile
    num_bamfiles = len(bamfiles)
    for bamfile in bamfiles:
        if not os.path.isfile(bamfile):
            error(f'Not a file: {bamfile}')
        if not os.path.isfile(bamfile + '.bai'):
            log(f'Bam index missing for file: {bamfile}. Trying "samtools index {bamfile}" ...')
            ret = os.system(f'samtools index {bamfile}')
            if ret != 0:
                log(f'Warning: samtools index returned exit code {ret} - everything might crash and burn.')
    log(f'Number of .bam files: {num_bamfiles}')

    # position to get coverage for
    ref_name = None
    pos_start = None
    pos_end = None

    with pysam.AlignmentFile(bamfiles[0], 'rb') as pysam_bamfile:
        # positions are 1-based (inclusive) in this block

        if args.position:
            if ':' in args.position:
                # with coords
                ref_name, startstop_str = args.position.split(':')
                pos_start_str, pos_end_str = startstop_str.split('-')
                pos_start = int(pos_start_str)
                pos_end = int(pos_end_str)

                if pos_start < 1:
                    error(f'Illegal start position: {pos_start}')
                if pos_start > pos_end:
                    error(f'Start position is greater than end position: {pos_start} vs {pos_end}')

            else:
                # no coords
                ref_name = args.position
                log(f'No coordinates given for reference {ref_name}, assuming whole reference. Inferring length from first alignment: {bamfiles[0]}')
                if ref_name not in pysam_bamfile.references:
                    error(f'Reference {ref_name} not found in alignment {bamfiles[0]}')
                # get start/stop
                pos_start = 1
                pos_end = pysam_bamfile.lengths[pysam_bamfile.references.index(ref_name)]

        else:
            # no position
            log(f'No position given, assuming whole reference. Taking first reference name from first alignment: {bamfiles[0]}')
            ref_name = pysam_bamfile.references[0]
            pos_start = 1
            pos_end = pysam_bamfile.lengths[pysam_bamfile.references.index(ref_name)]


    log(f'Position: {ref_name}:{pos_start}-{pos_end} - length: {pos_end-pos_start+1}')

    # convert from 1-based inclusive (genomic) to 0-based half open interval (pythonic)
    pos_start -= 1
    pos_len = pos_end - pos_start


    ###
    # find out how many processes we should use, make worker pool
    num_cpus = 2
    max_cpus = os.cpu_count()
    if not max_cpus:
        log('Could not determine number of available CPUs, choosing default maximum of 2.')
    else:
        if max_cpus > 2:
            num_cpus = max_cpus
        if num_bamfiles < num_cpus:
            num_cpus = num_bamfiles
        log(f'Using {num_cpus} of {max_cpus} available CPUs.')

    # make a numpy array to hold cov data
    cov_data = np.zeros([len(bamfiles), pos_len], dtype=int)
    process_results = []

    # make worker pool
    pool = Pool(num_cpus)

    log(f'Made worker pool with {num_cpus} processes.')
    log(f'Parsing alignments ...')
    parse_start_time = time.time()

    # send work to pool
    for bamfile in bamfiles:
        res = pool.apply_async(worker_get_coverage_ref_name, (bamfile, ref_name, pos_start, pos_len))
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


    # get into pandas dataframe
    data = pd.DataFrame(cov_data, index=bamfiles, columns=range(pos_start+1, pos_end+1)).T
    # print(data)


    ###
    # plotting
    log('Plotting ...')
    output_name = 'fastcov_output.pdf'

    # plot it
    sns.set(style="whitegrid")
    fig, ax = plt.subplots(figsize=(16, 8))

    p = sns.lineplot(data=data, linewidth=1.5, dashes=False)

    # formatting
    if args.logscale:
        plt.yscale('log')
    else:
        yup = plt.ylim()[1]
        plt.ylim(-yup/40, yup+(yup/40))

    
    # save figure
    plt.savefig(output_name, bbox_inches='tight')
    log(f'Wrote {output_name}')
    log(f'All done. Took {time.time() - main_start_time:.2f} seconds total. Have a nice day!')

