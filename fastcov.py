#!/usr/bin/env python3
# fastcov.py
# SK

from matplotlib import pyplot as plt
import os
import sys
import time
from multiprocessing import Pool
import numpy as np
import pandas as pd
import seaborn as sns
import pysam
import matplotlib as mpl

import warnings

from input_interpretation import parse_args, check_input, parse_or_infer_reference_and_position, BamFileChunk
from custom_logging import log

mpl.use('Agg')  # do not require X window
# suppress matplotlib warning, e.g. about legend
warnings.filterwarnings("ignore")


def fastcov_main():
    main_start_time = time.time()
    log('Started fastcov.py ...')

    ###
    # parse arguments
    args = parse_args()

    ###
    # params
    do_plots = True

    ###
    # input checking
    bam_files, num_bam_files = check_input(args)

    # position (chunk) to get coverage for
    with pysam.AlignmentFile(bam_files[0], 'rb') as pysam_bam_file:
        chunk = parse_or_infer_reference_and_position(
            args, bam_files, pysam_bam_file)

    log(f'Position: {chunk.reference}:{chunk.start}-{chunk.end} - length: {chunk.end-chunk.start+1}')

    ###
    num_cpus = determine_number_of_cpus_to_use(num_bam_files)

    # make a numpy array to hold cov data
    cov_data = np.zeros([len(bam_files), chunk.length], dtype=int)
    process_results = []

    # make worker pool
    pool = Pool(num_cpus)

    log(f'Made worker pool with {num_cpus} processes.')
    log(f'Parsing alignments ...')
    parse_start_time = time.time()
    total_num_aln = parse_bam_files(
        bam_files, chunk, cov_data, pool, process_results)
    parse_end_time = time.time()
    parse_time = parse_end_time - parse_start_time
    log(f'Parsed {total_num_aln} alignments in {parse_time:.2f} seconds')

    # get coverage data into pandas dataframe
    data = pd.DataFrame(cov_data, index=bam_files,
                        columns=range(chunk.start+1, chunk.end+1)).T

    # data output
    if args.csv_out:
        if args.csv_out == '-':
            log('Writing csv to stdout ...')
            data.to_csv(sys.stdout, header=(not args.csv_no_header))
            log(f'Wrote csv data to stdout.')
        else:
            log('Writing csv ...')
            data.to_csv(args.csv_out, header=(not args.csv_no_header))
            log(f'Wrote csv data to {args.csv_out}.')
        if not args.output_file:
            do_plots = False

    # plotting
    if do_plots:
        log('Plotting ...')
        output_name = plot(args, chunk, data)
        log(f'Wrote plot to {output_name}')

    log(f'All done. Took {time.time() - main_start_time:.2f} seconds total. Have a nice day!')


def determine_number_of_cpus_to_use(num_bam_files):
    # find out how many processes we should use, make worker pool
    num_cpus = 2
    max_cpus = os.cpu_count()
    if not max_cpus:
        log('Could not determine number of available CPUs, choosing default maximum of 2.')
    else:
        if max_cpus > 2:
            num_cpus = max_cpus
        if num_bam_files < num_cpus:
            num_cpus = num_bam_files
        log(f'Using {num_cpus} of {max_cpus} available CPUs.')
    return num_cpus


def parse_bam_files(bam_files, chunk, cov_data, pool, process_results):
    # send work to pool
    for bam_file in bam_files:
        res = pool.apply_async(worker_get_coverage_ref_name, (bam_file, chunk))
        process_results.append(res)
    # wait for results
    total_num_aln = 0
    for i, res in enumerate(process_results):
        # write to array
        num_aln, cov_data[i] = res.get()
        total_num_aln += num_aln
    return total_num_aln


def worker_get_coverage_ref_name(bam_file: str, chunk: BamFileChunk):

    log(f'Process {os.getpid()} working on {bam_file} ...')
    start_time = time.time()
    with pysam.AlignmentFile(bam_file, 'rb') as pysam_bamfile:

        # fetch alignments
        alignments = pysam_bamfile.fetch(
            chunk.reference, chunk.start, chunk.start+chunk.length)
        num_aln = 0

        coverage = np.zeros(chunk.length, dtype=int)

        for aln in alignments:
            num_aln += 1
            if aln.is_secondary:
                continue

            # record coverage
            for pos in aln.positions:
                pos_shift = pos - chunk.start
                # early stop
                if pos_shift >= chunk.length:
                    break

                if 0 <= pos_shift < chunk.length:
                    coverage[pos_shift] += 1

    log(f'{bam_file}: {num_aln} alignments done in {time.time() - start_time:.2f} sec.')

    return num_aln, coverage


def plot(args, chunk, data) -> str:

    output_name = 'fastcov_output.pdf'
    if args.output_file:
        output_name = args.output_file

    # plot it
    sns.set(style="whitegrid")
    fig, ax = plt.subplots(figsize=(20, 9))
    sns.lineplot(data=data, linewidth=1.5, dashes=False, alpha=0.8, legend=False)

    # formatting and annotation
    ydn, yup = plt.ylim()
    if args.logscale:
        data_max = data.max().max()
        next_pow_10 = 10
        while next_pow_10 <= data_max:
            next_pow_10 *= 10
        plt.ylim(0.5, next_pow_10 * 2)
        plt.yscale('log')
    else:
        plt.ylim(-yup / 40, yup + (yup / 40))
    plt.vlines([chunk.start + .5, chunk.end + .5], *plt.ylim(),
               'grey', linestyles='dashed', linewidth=.5)
    plt.title(f'Coverage at {chunk.reference}:{chunk.start + 1}-{chunk.end}')
    plt.ylabel('Coverage')
    plt.xlabel('Reference position')
    plt.text(0.045, 0.065, f'{(chunk.start + 1):,}', rotation=90,
             ha='right', va='bottom', transform=ax.transAxes)
    plt.text(0.956, 0.065, f'{(chunk.end):,}', rotation=90,
             ha='left', va='bottom', transform=ax.transAxes)
    
    xlims = plt.xlim()
    locs, labels = plt.xticks()
    plt.xticks(locs, [f'{int(pos):,}' for pos in locs])
    plt.xlim(xlims)

    short_names = [name.rsplit('/', 1)[-1][:-4] if name.endswith('.bam') else name.rsplit('/', 1)[-1]  for name in data.columns]
    plt.legend(labels=short_names)

    # save figure
    plt.savefig(output_name, bbox_inches='tight')
    return output_name


if __name__ == '__main__':
    fastcov_main()
