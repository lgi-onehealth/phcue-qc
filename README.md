# ph-cue: a fast and simple summary of FASTQ files

## Background

`ph-cue` (pronounced `fq`) is a simple tool for summarizing FASTQ files for quality control
purposes.

## Running

```bash
ph-cue [-t int] <fastq-file> > <summary_text>
```

Get help:

```bash
ph-cue -h
```

Get the current version:
```bash
ph-cue --version
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
