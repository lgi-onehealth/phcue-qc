extern crate bloom;
extern crate fastq;
extern crate num_format;

use std::path::PathBuf;
use std::sync::{Arc, Mutex};
// use std::sync::atomic::{AtomicU32, AtomicU64};
// use std::collections::HashMap;
use std::time::{Instant};
use structopt::StructOpt;
// use bloom::{ASMS, BloomFilter};
use fastq::{parse_path, Record};
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

fn main() {
    let opt = Opt::from_args();
    eprintln!("Working on: {:?}", opt.input);
    let fastq_stats = Arc::new(Mutex::new(FastqStats::new()));
    // let mut bf: BloomFilter = BloomFilter::with_rate(0.01, 10_000_000);
    // let mut kmer_counter: HashMap<String, u64> = HashMap::with_capacity(10_000_000);
    eprintln!("Starting to parse file...");
    let start = Instant::now();
    parse_path(Some(&opt.input), |parser| {
        let n_threads = 4;
        let local_stats = Arc::clone(&fastq_stats);
        let _results: Vec<u32> = parser.parallel_each(n_threads, move |record_sets| {
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
                }
            }
            let mut ls = local_stats.lock().unwrap();
            ls.update(total_reads, total_bases, min_read_len, max_read_len, total_qual);
            total_reads
        }).expect("Invalid FASTQ file.");
        let qual_offset = 33.0;
        let elapsed = start.elapsed().as_secs() as f64;
        let fq_lock = fastq_stats.lock().unwrap();
        fq_lock.print_stats(Some(qual_offset));
        println!("Elapsed time: {} seconds", elapsed);
        println!("Processing rate: {} reads/second", fq_lock.num_reads as f64 / elapsed);
    }).expect("Invalid compression");
}