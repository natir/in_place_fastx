use fastmap::fastq::Parser as fastmapParser;

type Counter<T, const N: usize> = [T; N];
trait AbsCounter {
    fn new() -> Self;
}

impl AbsCounter for Counter<std::sync::atomic::AtomicU64, 1024> {
    fn new() -> Self {
        unsafe { std::mem::transmute([0_u64; 1024]) }
    }
}

const K: u32 = 5;
const KMER_SPACE: usize = 2_usize.pow(K * 2);

struct ParserParallel {
    pub counter: std::collections::HashMap<u64, std::sync::atomic::AtomicU64>,
}

impl fastmap::fastq::Parser for ParserParallel {}

fn worker(
    record: fastmap::fastq::Record,
    data: &std::collections::HashMap<u64, std::sync::atomic::AtomicU64>,
) {
    for kmer in cocktail::tokenizer::Tokenizer::new(record.1, K as u8) {
        data.get(&kmer)
            .unwrap()
            .fetch_add(1, std::sync::atomic::Ordering::SeqCst);
    }
}

fn main() -> fastmap::error::Result<()> {
    let mut parser = ParserParallel {
        counter: std::collections::HashMap::new(),
    };

    for kmer in 0..KMER_SPACE {
        parser
            .counter
            .insert(kmer as u64, std::sync::atomic::AtomicU64::new(0));
    }

    let mut args = std::env::args();
    let _ = args.next();
    for input in args {
        parser.multithread_by_block_with_blocksize(131072, input, &parser.counter, worker)?;
    }

    for (kmer, count) in parser.counter.iter() {
        println!(
            "{},{}",
            cocktail::kmer::kmer2seq(*kmer, K as u8),
            count.load(std::sync::atomic::Ordering::SeqCst)
        );
    }

    Ok(())
}
