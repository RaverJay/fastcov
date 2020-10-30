# fastcov
Just get me some coverage plots, dude

# Quick start
`fastcov.py my_alignment.bam another_alignment.bam`

generates a coverage plot called `fastcov_output.pdf`

# Dependencies
makes use of
* `numpy`
* `pandas`
* `pysam`
* `seaborn`

# Full usage for now
```
usage: fastcov.py [-h] [-p POSITION] [-l] [-o OUTFILE] bamfile [bamfile ...]

Plot the coverage based on some bam files.

positional arguments:
  bamfile               Alignment files to include in the coverage plot.

optional arguments:
  -h, --help            show this help message and exit
  -p POSITION, --position POSITION
                        Specify a genomic position to plot exclusively.
                        Format: <ref_name>[:<start>-<stop>] (i.e. coordinates
                        are optional and must be 1-based and inclusive)
  -l, --logscale        Use logarithmic scale on y-axis.
  -o OUTFILE, --outfile OUTFILE
                        Specify output filename. File extension defines the
                        format (default: fastcov_output.pdf)
```
