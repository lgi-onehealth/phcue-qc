extern crate bloom;
extern crate fastq;
use std::path::PathBuf;
// use std::collections::HashMap;
use std::time::{Instant};
use structopt::StructOpt;
// use bloom::{ASMS, BloomFilter};
use fastq::{parse_path, Record};

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
    eprintln!("Working on: {:?}", opt.input);
    // let mut bf: BloomFilter = BloomFilter::with_rate(0.01, 10_000_000);
    // let mut kmer_counter: HashMap<String, u64> = HashMap::with_capacity(10_000_000);
    eprintln!("Starting to parse file...");
    let start = Instant::now();
    parse_path(Some(&opt.input), |parser| {
        let n_threads = 4;
        let results: Vec<(u32,u64,u64)> = parser.parallel_each(n_threads, move |record_sets| {
            let mut total_reads: u32 = 0;
            let mut total_bases: u64 = 0;
            let mut total_qual: u64 = 0;
            for record_set in record_sets {
                for _record in record_set.iter() {
                    total_reads += 1;
                    total_bases += _record.seq().len() as u64;
                    total_qual += _record.qual().iter().map(|x| *x as u64).sum::<u64>();
        }
            }
            (total_reads, total_bases, total_qual)
        }).expect("Invalid FASTQ file.");
        let qual_offset = 33.0;
        let elapsed = start.elapsed().as_secs() as f64;
        let total_reads: u32 = results.iter().map(|x| x.0).sum();
        let total_bases: u64 = results.iter().map(|x| x.1).sum();
        let mean_q_score: f64 = (results.iter().map(|x| x.2).sum::<u64>() as f64 / total_bases as f64) - qual_offset;
        let reads_per_second = total_reads as f64 / elapsed;
        println!("Total reads: {} ({} reads/second)", total_reads, reads_per_second);
        println!("Total bases: {}", total_bases);
        println!("Total bases (parallel): {}", par_total_nucleotides.load(std::sync::atomic::Ordering::Relaxed));
        println!("Mean quality score: {}", mean_q_score);
        println!("Total time spent: {}s", elapsed);
        println!("Total readsets processed: {}", results.len());
        println!("The results vector: {:?}", results);
    }).expect("Invalid compression");
}