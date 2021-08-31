use in_place_fastx;

use rayon::iter::ParallelBridge;
use rayon::iter::ParallelIterator;

const K: usize = 13;
const KMER_SPACE: usize = 4_usize.pow(K as u32);
const BLOCK_SIZE: u64 = 2u64.pow(16);

in_place_fastx::fastq_sharedstate!(
    Parser,
    std::collections::HashMap<u64, std::sync::atomic::AtomicU64>,
    |record: in_place_fastx::block::Record, counter: &std::collections::HashMap<u64, std::sync::atomic::AtomicU64>| {
    if record.sequence.len() < K {
        return ;
    }

    for kmer in cocktail::tokenizer::Tokenizer::new(record.sequence, K as u8) {
        counter.get(&kmer).unwrap().fetch_add(1, std::sync::atomic::Ordering::SeqCst);
    }
    }
);

fn main() -> in_place_fastx::error::Result<()> {
    let mut counter = std::collections::HashMap::with_capacity(KMER_SPACE);
    let parser = Parser::new();

    for kmer in 0..KMER_SPACE {
        counter.insert(kmer as u64, std::sync::atomic::AtomicU64::new(0));
    }

    let mut args = std::env::args();
    let _ = args.next();
    for input in args {
        parser.with_blocksize(BLOCK_SIZE, input, &counter)?;
    }

    for (kmer, count) in counter.iter() {
        println!(
            "{},{}",
            cocktail::kmer::kmer2seq(*kmer, K as u8),
            count.load(std::sync::atomic::Ordering::SeqCst)
        );
    }

    Ok(())
}
