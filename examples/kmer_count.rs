use in_place_fastx::block;
use in_place_fastx::fastq;
use in_place_fastx::parser::Sequential;

type Counter<T, const N: usize> = [T; N];
trait AbsCounter {
    fn new() -> Self;
}

impl<const N: usize> AbsCounter for Counter<u64, N> {
    fn new() -> Self {
        [0; N]
    }
}

const K: u32 = 5;
const KMER_SPACE: usize = 2_usize.pow(K * 2);

struct Parser {
    pub counter: std::collections::HashMap<u64, u64>,
}

impl Sequential<fastq::Producer, fastq::Reader> for Parser {
    fn record(&mut self, record: block::Record) {
        if record.sequence.len() < K as usize {
            return;
        }
        for kmer in cocktail::tokenizer::Tokenizer::new(record.sequence, K as u8) {
            (*self.counter.get_mut(&kmer).unwrap()) += 1;
        }
    }
}

fn main() -> in_place_fastx::error::Result<()> {
    let mut parser = Parser {
        counter: std::collections::HashMap::new(),
    };

    for kmer in 0..KMER_SPACE {
        parser.counter.insert(kmer as u64, 0);
    }

    let mut args = std::env::args();
    let _ = args.next();
    for input in args {
        parser.with_blocksize(1048576, input)?;
    }

    for (kmer, count) in parser.counter.iter() {
        println!("{},{}", cocktail::kmer::kmer2seq(*kmer, K as u8), count);
    }

    Ok(())
}
