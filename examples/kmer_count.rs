use in_place_fastx;

const K: usize = 13;
const KMER_SPACE: usize = 4_usize.pow(K as u32);
const BLOCK_SIZE: u64 = 2u64.pow(18);

in_place_fastx::fastq_sequential!(
    Parser,
    Vec<u8>,
    |record: in_place_fastx::block::Record, counter: &mut Vec<u8>| {
        if record.sequence.len() < K {
            return;
        }

        for kmer in cocktail::tokenizer::Canonical::new(record.sequence, K as u8) {
            counter[(kmer as usize) >> 1] += 1;
        }
    }
);

fn main() -> in_place_fastx::error::Result<()> {
    let mut counter: Vec<u8> = vec![0; cocktail::kmer::get_hash_space_size(K as u8) as usize];

    let mut parser = Parser::new();

    let mut args = std::env::args();
    let _ = args.next();
    for input in args {
        parser.with_blocksize(BLOCK_SIZE, input, &mut counter)?;
    }

    for (kmer, count) in counter.iter().enumerate() {
        println!(
            "{},{}",
            cocktail::kmer::kmer2seq(kmer as u64, K as u8),
            count
        );
    }

    Ok(())
}
