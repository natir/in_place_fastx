use fastmap::fastq::Parser as fastmapParser;

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

impl fastmap::fastq::Parser for Parser {
    fn record(&mut self, record: fastmap::fastq::Record) {
        for kmer in cocktail::tokenizer::Tokenizer::new(record.1, K as u8) {
            (*self.counter.get_mut(&kmer).unwrap()) += 1;
        }
    }
}

fn main() -> fastmap::error::Result<()> {
    let mut parser = Parser {
        counter: std::collections::HashMap::new(),
    };

    for kmer in 0..KMER_SPACE {
        parser.counter.insert(kmer as u64, 0);
    }

    let mut args = std::env::args();
    let _ = args.next();
    for input in args {
        parser.file_with_blocksize(131072, input)?;
    }

    for (kmer, count) in parser.counter.iter() {
        println!("{},{}", cocktail::kmer::kmer2seq(*kmer, K as u8), count);
    }

    Ok(())
}
