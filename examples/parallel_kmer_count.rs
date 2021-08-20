use in_place_fastx::fastq::Parser as in_place_fastxParser;

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

impl in_place_fastx::fastq::Parser for ParserParallel {}

fn worker(
    record: in_place_fastx::fastq::Record,
    data: &std::collections::HashMap<u64, std::sync::atomic::AtomicU64>,
) {
    if record.1.len() < K as usize{
	return
    }
    for kmer in cocktail::tokenizer::Tokenizer::new(record.1, K as u8) {
        data.get(&kmer)
            .unwrap()
            .fetch_add(1, std::sync::atomic::Ordering::SeqCst);
    }
}

fn main() -> in_place_fastx::error::Result<()> {
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
        parser.multithread_by_block_with_blocksize(1048576, input, &parser.counter, worker)?;
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
