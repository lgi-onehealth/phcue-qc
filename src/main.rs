extern crate bio;
extern crate bloom;
extern crate dashmap;
extern crate fastq;
extern crate num_format;

use std::path::PathBuf;
use std::sync::{Arc, Mutex};
// use std::sync::atomic::{AtomicU32, AtomicU64};
// use std::collections::HashMap;
use std::time::{Instant};
use structopt::StructOpt;
use bio::alphabets::dna;
use dashmap::DashMap;
// use bloom::{ASMS, BloomFilter};
use fastq::{parse_path, Record, RefRecord};
use num_format::{Locale, ToFormattedString};

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

#[derive(Debug, Copy, Clone)]
struct FastqStats {
    num_reads: u32,
    num_bases: u64,
    min_read_len: u16,
    max_read_len: u16,
    total_q_score: u64,
}

impl FastqStats {
    fn new() -> Self {
        FastqStats {
            num_reads: 0,
            num_bases: 0,
            min_read_len: std::u16::MAX,
            max_read_len: 0,
            total_q_score: 0,
        }
    }

    fn update(&mut self, local_num_reads: u32, local_num_bases: u64, 
        local_min_read_len: u16, local_max_read_len: u16, local_total_q_score: u64) {
        self.num_reads += local_num_reads;
        self.num_bases += local_num_bases;
        self.min_read_len = std::cmp::min(self.min_read_len, local_min_read_len);
        self.max_read_len = std::cmp::max(self.max_read_len, local_max_read_len);
        self.total_q_score += local_total_q_score;
    }

    fn print_stats(self, offset: Option<f64>) {
        println!("Total reads: {}", self.num_reads.to_formatted_string(&Locale::en));
        println!("Total bases: {}", self.num_bases.to_formatted_string(&Locale::en));
        println!("Min read len: {}", self.min_read_len);
        println!("Max read len: {}", self.max_read_len);
        println!("Mean Q score: {:.2}", self.mean_q_score(offset));
    }

    fn _get(self) -> (u32, u64, u16, u16, u64, u64) {
        (self.num_reads, self.num_bases, self.min_read_len, self.max_read_len, self.total_q_score, self.total_q_score)
    }

    fn mean_q_score(self, offset: Option<f64>) -> f64 {
        (self.total_q_score as f64 / self.num_bases as f64) - offset.unwrap_or(33.0)
    }    
}


trait Kmers {
    fn get_kmers<'a> (&self, k: usize) ->  Vec<Vec<u8>>;
}

impl Kmers for RefRecord<'_> {
    
    fn get_kmers<'a> (&self, k: usize) ->  Vec<Vec<u8>> {
        let seq = self.seq().iter().map(|b| *b).collect::<Vec<u8>>();
        let comp = seq.iter().map(|b| dna::complement(*b)).collect::<Vec<u8>>();
        let iter = seq.windows(k).zip(comp.windows(k));
        // Box::new(iter.map(|ks| std::borrow::Cow::from(std::cmp::min(ks.0, ks.1).copy())))
        iter.map(|ks| {
            let a = ks.0.iter().map(|n| *n).collect::<Vec<u8>>();
            let b = ks.1.iter().map(|n| *n).collect::<Vec<u8>>();
            std::cmp::min(a, b)}
        ).collect::<Vec<Vec<u8>>>()
        
    }
}

fn main() {
    let opt = Opt::from_args();
    eprintln!("Working on: {:?}", opt.input);
    let fastq_stats = Arc::new(Mutex::new(FastqStats::new()));
    let kmers = Arc::new(DashMap::with_capacity(10_000_000));
    // let mut bf: BloomFilter = BloomFilter::with_rate(0.01, 10_000_000);
    // let mut kmer_counter: HashMap<String, u64> = HashMap::with_capacity(10_000_000);
    eprintln!("Starting to parse file...");
    let start = Instant::now();
    parse_path(Some(&opt.input), |parser| {
        let n_threads = 8;
        let local_stats = Arc::clone(&fastq_stats);
        let local_kmers = Arc::clone(&kmers);
        let _results: Vec<u32> = parser.parallel_each(n_threads, move |record_sets| {
            let msg = String::from("Starting to process new batch of reads...");
            eprintln!("{}", msg);
            let kmer_len = 21;
            let mut total_reads: u32 = 0;
            let mut total_bases: u64 = 0;
            let mut total_qual: u64 = 0;
            let mut min_read_len: u16 = std::u16::MAX;
            let mut max_read_len: u16 = 0;
            for record_set in record_sets {
                for _record in record_set.iter() {
                    total_reads += 1;
                    total_bases += _record.seq().len() as u64;
                    total_qual += _record.qual().iter().map(|x| *x as u64).sum::<u64>();
                    min_read_len = std::cmp::min(min_read_len, _record.seq().len() as u16);
                    max_read_len = std::cmp::max(max_read_len, _record.seq().len() as u16);
                    let rc_win = dna::revcomp(_record.seq()).into_iter().rev().collect::<Vec<u8>>();
                    let seq_win = _record.seq().windows(kmer_len);
                    seq_win.zip(rc_win.windows(kmer_len)).for_each(|ks| {
                        let k1 = ks.0.clone().to_vec();
                        let mut k2 = ks.1.clone().to_vec();
                        k2.reverse();
                        *local_kmers.entry(std::cmp::min(k1, k2)).or_insert(0) += 1});
                }
            }
            let mut ls = local_stats.lock().unwrap();
            ls.update(total_reads, total_bases, min_read_len, max_read_len, total_qual);
            println!("Finalizing batch...");
            total_reads
        }).expect("Invalid FASTQ file.");
        println!("Finished processing file.");
        let qual_offset = 33.0;
        let elapsed = start.elapsed().as_millis() as f64 * 0.001;
        let fq_lock = fastq_stats.lock().unwrap();
        fq_lock.print_stats(Some(qual_offset));
        println!("Total kmers: {}", kmers.len().to_formatted_string(&Locale::en));
        for i in 1..50 {
            let num_kmers = kmers.iter().filter(|x| *x.value() == i).count().to_formatted_string(&Locale::en);
            println!("Unique kmers with count of {}: {}", i, num_kmers);    
        }
        let total_count = kmers.iter().map(|x| *x.value()).sum::<u64>();
        let total_unique_kmers = kmers.iter().filter(|x| *x.value() > 2).count().to_formatted_string(&Locale::en);
        println!("Total count: {}", total_count.to_formatted_string(&Locale::en));
        println!("Total unique: {}", total_unique_kmers);
        println!("Elapsed time: {} seconds", elapsed);
        let proc_rate = fq_lock.num_reads as f64 / elapsed;

        println!("Processing rate: {:.2} reads/second", proc_rate);
        // for k in kmers.iter() {
        //     let kmer = k.key();
        //     let count = k.value();
        //     let kmer_str = String::from_utf8(kmer.to_vec()).unwrap();
        //     let kmer_str = kmer_str.replace("\0", "");
        //     println!("{}\t{}", kmer_str, count.to_formatted_string(&Locale::en));
        // }
    }).expect("Invalid compression");
}