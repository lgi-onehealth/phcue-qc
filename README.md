# phcue-qc: a fast and simple summary of FASTQ files

## Background

`phcue-qc` (pronounced `fq qc`) is a simple tool for summarizing FASTQ files for quality control
purposes.

## Running

```bash
phcue-qc [-t int] [-m int] <fastq-file> > <summary_text>
```

The required input:
* `fastq-file`: a FASTQ file, typically the R1 for paired-end sequencing where we would have fewer errors.

The options include:
* `-t` or `--threads`: the number of threads to use for parallelization.
* `-m` or `--min-count`: the minimum number of times a kmer must be observed to be included in the summary.

Get help:

```bash
phcue-qc -h

phcue-qc 0.2.0
A fast and simple summary of FASTQ.

USAGE:
    phcue-qc [OPTIONS] <input>

FLAGS:
    -h, --help       Prints help information
    -V, --version    Prints version information

OPTIONS:
    -m, --min-count <min-count>     [default: 3]
    -t, --threads <threads>         [default: 1]

ARGS:
    <input>    Input file
```

Get the current version:
```bash
phcue-qc --version
```

## The output

An example output:

```
Total reads: 636,311
Total bases: 160,520,099
Min read len: 35
Max read len: 301
Mean Q score: 36.09
Total kmers: 10,279,035
Total count: 147,771,738
Total unique: 4,807,539
Elapsed time: 6.555 seconds
Processing rate: 97072.62 reads/second
```
