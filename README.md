# fastcov.py

For when you just need a coverage plot.

Just throw in your indexed bamfiles, done. Will process the files in parallel. Uses the first reference of the first file if you don't give a specific genomic position.

Feature requests welcome!

## Quick start
```
fastcov.py my_alignment.bam another_alignment.bam
```
generates a coverage plot called `fastcov_output.pdf`

## Example

Get a logplot of some SARS-CoV-2 direct RNA sequencing data

`./fastcov.py -l -o example.png flowcell*.bam`

![logplot of SARS-CoV-2 coverages](/images/example.png)

## Dependencies
makes use of
* `pysam`
* `numpy`
* `pandas`
* `seaborn`

## Full usage
```
usage: fastcov.py [-h] [-p POSITION] [-l] [-o OUTPUT_FILE] [-c CSV_OUT]
                  [--csv_no_header]
                  bamfile [bamfile ...]

Plot the coverage based on some bam files.

positional arguments:
  bamfile               Alignment files to include in the coverage plot.

optional arguments:
  -h, --help            show this help message and exit
  -p POSITION, --position POSITION
                        Specify a genomic position to plot exclusively.
                        Format: <ref_name>[:<start>-<stop>] Coordinates are
                        1-based and inclusive. Start and/or stop are optional
                        with fallbacks 1 and <length_of_ref> respectively
                        (i.e. 'chr1', 'chr1:-200', 'chr1:100-' and
                        'chr1:100-200 are legal)
  -l, --logscale        Use logarithmic scale on y-axis.
  -o OUTPUT_FILE, --output_file OUTPUT_FILE
                        Specify plot output filename. File extension defines
                        the format (default: fastcov_output.pdf)
  -c CSV_OUT, --csv_out CSV_OUT
                        Specify csv data output filename. Use '-' to write to
                        stdout. Will disable plot output by default, specify
                        --output_file to re-enable plot output.
  --csv_no_header       Suppress column names in csv output.
```
