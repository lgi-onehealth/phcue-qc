extern crate dashmap;
extern crate fastq;
extern crate hdrhistogram;
extern crate num_format;

use clap::Parser;
use dashmap::DashMap;
use fastq::{parse_path, Record};
use hdrhistogram::Histogram;
use num_format::{Locale, ToFormattedString};
use serde::{Deserialize, Serialize};
// use serde_json::Result;
use find_peaks::PeakFinder;
use std::ops::Deref;
use std::path::PathBuf;
use std::sync::{Arc, Mutex};
use std::time::Instant;

#[derive(Debug, Parser)]
#[clap(
    name = "ph-cue",
    version,
    about = "A fast and simple summary of FASTQ."
)]
struct Opt {
    /// Input file
    #[clap(parse(from_os_str))]
    input: PathBuf,

    // Number of threads
    #[clap(short = 't', long = "threads", default_value = "1")]
    threads: usize,

    // Minimum kmer count to consider a kmer
    #[clap(short = 'm', long = "min-count", default_value = "3")]
    min_count: u64,

    // Whether the experiment is single_end
    #[clap(short = 's', long = "single-end")]
    single_end: bool,

    // Whether to output the results in JSON
    #[clap(short = 'j', long = "json")]
    json: bool,
}

#[derive(Debug, Copy, Clone, Serialize, Deserialize)]
struct FastqStats {
    num_reads: u32,
    num_bases: u64,
    min_read_len: u16,
    max_read_len: u16,
    #[serde(skip)]
    total_q_score: u64,
    depth_coverage: u64,
    genome_size: u64,
    mean_q_score: f64,
    #[serde(skip)]
    multiplier: u8,
}

impl FastqStats {
    fn new() -> Self {
        FastqStats {
            num_reads: 0,
            num_bases: 0,
            min_read_len: std::u16::MAX,
            max_read_len: 0,
            total_q_score: 0,
            depth_coverage: 0,
            genome_size: 0,
            multiplier: 2,
            mean_q_score: 0.0,
        }
    }

    fn update(
        &mut self,
        local_num_reads: u32,
        local_num_bases: u64,
        local_min_read_len: u16,
        local_max_read_len: u16,
        local_total_q_score: u64,
    ) {
        self.num_reads += local_num_reads;
        self.num_bases += local_num_bases;
        self.min_read_len = std::cmp::min(self.min_read_len, local_min_read_len);
        self.max_read_len = std::cmp::max(self.max_read_len, local_max_read_len);
        self.total_q_score += local_total_q_score;
    }

    fn update_multiplier(&mut self, single_end: bool) {
        self.multiplier = if single_end { 1 } else { 2 };
    }

    fn get_depth_coverage(&mut self, hist: &Histogram<u64>, min_count: &u64) {
        // let mut depth_coverage: u64 = 0;
        // let mut count_value: u64 = 0;
        let counts: Vec<i64> = hist
            .iter_recorded()
            .map(|x| x.count_at_value() as i64)
            .collect();
        let values: Vec<u64> = hist
            .iter_recorded()
            .map(|x| x.value_iterated_to())
            .collect();
        let mut pf = PeakFinder::new(&counts);
        pf.with_min_prominence(10);
        pf.with_min_height(10);
        let peaks = pf.find_peaks();
        // for v in hist.iter_recorded() {
        //     if v.value_iterated_to() < *min_count {
        //         continue;
        //     }
        //     if v.count_at_value() > count_value {
        //         depth_coverage = v.value_iterated_to();
        //         count_value = v.count_at_value();
        //         println!(
        //             "depth_coverage: {}\t count_value: {}",
        //             depth_coverage, count_value
        //         );
        //     }
        // }
        self.depth_coverage = values[peaks[0].middle_position() as usize];
    }

    fn get_genome_size(&mut self, kmers: &Arc<DashMap<u64, u64>>, min_count: &u64) {
        self.genome_size = kmers.iter().filter(|x| *x.value() > *min_count).count() as u64;
    }

    fn print_stats(self) {
        println!(
            "Total reads: {}",
            (self.multiplier as u32 * self.num_reads).to_formatted_string(&Locale::en)
        );
        println!(
            "Total bases: {}",
            (self.multiplier as u64 * self.num_bases).to_formatted_string(&Locale::en)
        );
        println!(
            "Estimated genome size: {}",
            self.genome_size.to_formatted_string(&Locale::en)
        );
        println!("Min read len: {}", self.min_read_len);
        println!("Max read len: {}", self.max_read_len);
        println!("Mean Q score: {:.2}", self.mean_q_score);
        println!(
            "Depth of coverage: {}",
            self.multiplier as u64 * self.depth_coverage
        );
    }

    fn _get(self) -> (u32, u64, u16, u16, u64, u64, u64) {
        (
            self.num_reads,
            self.num_bases,
            self.min_read_len,
            self.max_read_len,
            self.total_q_score,
            self.total_q_score,
            self.depth_coverage,
        )
    }

    fn get_mean_q_score(&mut self, offset: Option<f64>) -> f64 {
        if self.mean_q_score == 0.0 {
            self.mean_q_score =
                (self.total_q_score as f64 / self.num_bases as f64) - offset.unwrap_or(33.0)
        }
        return self.mean_q_score;
    }
}

fn transform(base: char) -> u64 {
    match base {
        'A' | 'a' => 0,
        'C' | 'c' => 1,
        'G' | 'g' => 2,
        'T' | 't' => 3,
        _ => 4,
    }
}

fn main() {
    let opt = Opt::from_args();
    eprintln!("Working on: {:?}", opt.input);
    let _filename = opt.input.clone().file_name().unwrap().to_str().unwrap();
    let fastq_stats = Arc::new(Mutex::new(FastqStats::new()));
    let kmers = Arc::new(DashMap::with_capacity(10_000_000));
    eprintln!("Starting to parse file...");
    let start = Instant::now();
    parse_path(Some(&opt.input), |parser| {
        let n_threads = opt.threads;
        let local_stats = Arc::clone(&fastq_stats);
        let local_kmers = Arc::clone(&kmers);
        let _results: Vec<u32> = parser
            .parallel_each(n_threads, move |record_sets| {
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
                        let mut r = 0;
                        let mut a = 0;
                        let mask: u64 = (1 << 2 * kmer_len) - 1;
                        let shift = (kmer_len - 1) * 2;
                        for (i, c) in _record.seq().to_vec().iter().enumerate() {
                            let c = *c as char;
                            let c = transform(c);
                            if c < 4 {
                                r = (r << 2 | c) & mask;
                                a = a >> 2 | (3 - c) << shift;
                                if i >= kmer_len - 1 {
                                    let canonical_kmer = std::cmp::min(r, a);
                                    *local_kmers.entry(canonical_kmer).or_insert(0) += 1;
                                }
                            } else {
                                r = 0;
                                a = 0;
                            }
                        }
                    }
                }
                let mut ls = local_stats.lock().unwrap();
                ls.update(
                    total_reads,
                    total_bases,
                    min_read_len,
                    max_read_len,
                    total_qual,
                );
                eprintln!("Finalizing batch...");
                total_reads
            })
            .expect("Invalid FASTQ file.");
        eprintln!("Finished processing file.");
        let qual_offset = 33.0;
        let elapsed = start.elapsed().as_millis() as f64 * 0.001;
        let mut hist = Histogram::<u64>::new(2).unwrap();
        kmers.iter().filter(|c| *c.value() > 1).for_each(|c| {
            hist += *c.value() as u64;
        });
        let mut fq_lock = fastq_stats.lock().unwrap();
        fq_lock.update_multiplier(opt.single_end);
        fq_lock.get_depth_coverage(&hist, &opt.min_count);
        fq_lock.get_genome_size(&kmers, &opt.min_count);
        fq_lock.get_mean_q_score(Some(qual_offset));
        if opt.json {
            let serialized = serde_json::to_string_pretty(fq_lock.deref()).unwrap();
            println!("{}", serialized);
        } else {
            fq_lock.print_stats();
        }
        eprintln!("Elapsed time: {} seconds", elapsed);
        let proc_rate = fq_lock.num_reads as f64 / elapsed;

        eprintln!("Processing rate: {:.2} reads/second", proc_rate);
    })
    .expect("Invalid compression");
}
