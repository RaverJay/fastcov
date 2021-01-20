import argparse
import os
from typing import NamedTuple
from custom_logging import error, log


class BamFileChunk(NamedTuple):
    reference: str
    start: int
    end: int
    length: int


def parse_args():
    parser = argparse.ArgumentParser(
        description='Plot the coverage based on some bam files.')
    parser.add_argument("bamfile",
                        nargs='+',
                        help="Alignment files to include in the coverage plot.")
    parser.add_argument("-p", "--position",
                        help="Specify a genomic position to plot exclusively. Format: <ref_name>[:<start>-<stop>]\n"
                             "Coordinates are 1-based and inclusive. "
                             "Start and/or stop are optional with fallbacks 1 and <length_of_ref> respectively "
                             "(i.e. 'chr1', 'chr1:-200', 'chr1:100-' and 'chr1:100-200 are legal)")
    parser.add_argument("-l", "--logscale", action='store_true',
                        help="Use logarithmic scale on y-axis.")
    parser.add_argument("-o", "--output_file",
                        help="Specify plot output filename. File extension defines the format "
                             "(default: fastcov_output.pdf)")
    parser.add_argument("-c", "--csv_out",
                        help="Specify csv data output filename. Will disable plot output by default, "
                             "specify --output_file to re-enable plot output.")
    args = parser.parse_args()
    return args


def check_input(args):
    bam_files = []
    for bam_file in args.bamfile:
        if not os.path.isfile(bam_file):
            error(f'Not a file: {bam_file}')
        if bam_file in bam_files:
            log(f'WARNING: Skipping duplicate bam file input: {bam_file}')
            continue
        if not os.path.isfile(bam_file + '.bai'):
            log(
                f'Bam index missing for file: {bam_file}. Trying "samtools index {bam_file}" ...')
            ret = os.system(f'samtools index {bam_file}')
            if ret != 0:
                error(
                    f'ERROR: samtools index returned exit code {ret}')
        # checked
        bam_files.append(bam_file)
        
    num_bam_files = len(bam_files)
    log(f'Number of .bam files: {num_bam_files}')
    return bam_files, num_bam_files


def parse_or_infer_reference_and_position(args, bam_files, pysam_bam_file) -> BamFileChunk:
    # positions are 1-based (inclusive) in this block
    if args.position:
        if ':' in args.position:
            # with coords
            ref_name, start_stop_str = args.position.split(':')
            if '-' not in start_stop_str:
                error('Please provide a start and/or stop position. '
                      'When providing only a start or stop position, '
                      'indicate the intervals side to use by pre- or postfixing the value with "-" '
                      '(e.g. "100-" or "-200"). '
                      'The other side will be inferred from the given reference')
            pos_start_str, pos_end_str = start_stop_str.split('-')
            pos_start = convert_to_int_or_fallback(pos_start_str,
                                                   fallback=1)
            pos_end = convert_to_int_or_fallback(pos_end_str,
                                                 fallback=pysam_bam_file.lengths[pysam_bam_file.references.index(ref_name)])

            if pos_start < 1:
                error(f'Illegal start position: {pos_start}')
            if pos_start > pos_end:
                error(
                    f'Start position is greater than end position: {pos_start} vs {pos_end}')

        else:
            # no coords
            ref_name = args.position
            log(f'No coordinates given for reference {ref_name}, assuming whole reference. '
                f'Inferring length from first alignment: {bam_files[0]}')
            if ref_name not in pysam_bam_file.references:
                error(
                    f'Reference {ref_name} not found in alignment {bam_files[0]}')
            # get start/stop
            pos_start = 1
            pos_end = pysam_bam_file.lengths[pysam_bam_file.references.index(
                ref_name)]
    else:
        # no position
        log(
            f'No position given, assuming whole reference. Taking first reference name from first alignment: {bam_files[0]}')
        ref_name = pysam_bam_file.references[0]
        pos_start = 1
        pos_end = pysam_bam_file.lengths[pysam_bam_file.references.index(
            ref_name)]

    # convert from 1-based inclusive (genomic) to 0-based half open interval (pythonic)
    pos_start -= 1
    pos_len = pos_end - pos_start
    return BamFileChunk(end=pos_end, start=pos_start, reference=ref_name, length=pos_len)


def convert_to_int_or_fallback(string: str, fallback: int):
    try:
        return int(string)
    except ValueError:
        return fallback
