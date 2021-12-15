extern crate bloom;
extern crate needletail;
use std::path::PathBuf;
use std::collections::HashMap;
use std::time::{Instant};
use structopt::StructOpt;
use needletail::{parse_fastx_file, Sequence}; //, FastxReader};
use bloom::{ASMS, BloomFilter};

#[derive(Debug, StructOpt)]
#[structopt(name = "example", about = "An example of StructOpt usage.")]
struct Opt {
    /// Activate debug mode
    // short and long flags (-d, --debug) will be deduced from the field's name
    #[structopt(short, long)]
    debug: bool,

    /// Set speed
    // we don't want to name it "speed", need to look smart
    #[structopt(short = "v", long = "velocity", default_value = "42")]
    speed: f64,

    /// Input file
    #[structopt(parse(from_os_str))]
    input: PathBuf,

    /// Output file, stdout if not present
    #[structopt(parse(from_os_str))]
    output: Option<PathBuf>,

    /// File name: only required when `out-type` is set to `file`
    #[structopt(name = "FILE", required_if("out-type", "file"))]
    file_name: Option<String>,
}

fn main() {
    let opt = Opt::from_args();
    println!("Working on: {:?}", opt.input);
    let mut n_bases: u64 = 0;
    let mut total_qual: u64 = 0;
    let mut n_seqs: u64 = 0;
    let mut max_read_len: u16 = 0;
    let mut min_read_len: u16 = std::u16::MAX;
    let mut read_len: u16;
    let mut bf: BloomFilter = BloomFilter::with_rate(0.01, 10_000_000);
    // let mut n_valid_kmers = 0;
    let mut kmer_counter: HashMap<String, u64> = HashMap::with_capacity(10_000_000);
    let mut reader = parse_fastx_file(&opt.input).expect("valid path/file");
    eprintln!("Starting to parse file...");
    let start = Instant::now();
    while let Some(record) = reader.next() {
        let seqrec = record.expect("invalid record");
        // keep track of the total number of bases
        n_bases += seqrec.num_bases() as u64;
        // keep track of Q score
        total_qual += seqrec.qual().unwrap().iter().map(|&q| q as u64).sum::<u64>();
        // keep track of the number of sequences
        n_seqs += 1;
        // keep track of the max and min read lengths
        read_len = seqrec.num_bases() as u16;
        if read_len > max_read_len {
            max_read_len = read_len;
        } else if read_len < min_read_len {
            min_read_len = read_len;
        }
        // normalize to make sure all the bases are consistently capitalized and
        // that we remove the newlines since this is FASTA
        let norm_seq = seqrec.normalize(false);
        // we make a reverse complemented copy of the sequence first for
        // `canonical_kmers` to draw the complemented sequences from.
        let rc = norm_seq.reverse_complement();
        // now we keep track of the number of AAAAs (or TTTTs via
        // canonicalization) in the file; note we also get the position (i.0;
        // in the event there were `N`-containing kmers that were skipped)
        // and whether the sequence was complemented (i.2) in addition to
        // the canonical kmer (i.1)
        for (_, kmer, _) in norm_seq.canonical_kmers(31, &rc) {
            let kmer = String::from_utf8_lossy(kmer).into_owned();
            // if !bf.contains(&kmer) {
            //     bf.insert(&kmer);
            //     n_valid_kmers += 1;
            // }
            if !bf.contains(&kmer) {
                bf.insert(&kmer);
            } else {
                *kmer_counter.entry(kmer).or_insert(1) += 1;
            }
        }
        if n_seqs % 10_000 == 0 {
            let time_spent = start.elapsed().as_millis();
            eprintln!("Total reads processed: {} ({} reads/msecond - total time {} seconds)", n_seqs, n_seqs as f64/time_spent as f64, time_spent as f64/1000.0);
        }
        if n_seqs > 20_000 {
            break;
        }
    }
    println!("Total reads: {}", n_seqs);
    println!("Total nucleotides: {}", n_bases);
    println!("Mean Q score: {}", (total_qual as f64/n_bases as f64)-33.0);
    println!("Max read length: {}", max_read_len);
    println!("Min read length: {}", min_read_len);
    for i in 2..11 {
        let total_kmers = kmer_counter.values().filter(|&x| *x >= i).count();
        println!("Total kmers witn count of at least {}: {}", i, total_kmers);
    }
}