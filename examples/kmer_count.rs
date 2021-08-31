use in_place_fastx;

const K: usize = 13;
const KMER_SPACE: usize = 4_usize.pow(K as u32);
const BLOCK_SIZE: u64 = 2u64.pow(16);

in_place_fastx::fastq_sequential!(
    Parser,
    std::collections::HashMap<u64, u64>,
    |record: in_place_fastx::block::Record, counter: &mut std::collections::HashMap<u64, u64>| {
    if record.sequence.len() < K {
        return ;
    }

    for kmer in cocktail::tokenizer::Tokenizer::new(record.sequence, K as u8) {
        *counter.entry(kmer).or_insert(0) += 1;
    }
    }
);

fn main() -> in_place_fastx::error::Result<()> {
    let mut counter = std::collections::HashMap::with_capacity(KMER_SPACE);
    let mut parser = Parser::new();

    let mut args = std::env::args();
    let _ = args.next();
    for input in args {
        parser.with_blocksize(BLOCK_SIZE, input, &mut counter)?;
    }

    for (kmer, count) in counter.iter() {
        println!("{},{}", cocktail::kmer::kmer2seq(*kmer, K as u8), count);
    }

    Ok(())
}
