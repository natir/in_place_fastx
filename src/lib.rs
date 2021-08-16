pub mod error;
pub mod fastq;

#[cfg(test)]
mod tests {
    //use super::*;

    use std::io::Write;

    use rand::Rng;
    use rand::SeedableRng;

    fn generate_fastq(seed: u64) -> tempfile::NamedTempFile {
        let mut rng = rand::rngs::StdRng::seed_from_u64(seed);

        let mut file = tempfile::NamedTempFile::new().unwrap();

        let dna = [b'A', b'C', b'T', b'G'];
        let qual = (0..94).collect::<Vec<u8>>();

        for i in 0..10_000 {
            let dna_seq = (0..150)
                .map(|_| dna[rng.gen_range(0..4)] as char)
                .collect::<String>();
            let qual_seq = (0..150)
                .map(|_| (qual[rng.gen_range(0..94)] + 33) as char)
                .collect::<String>();

            writeln!(file, "@{}\n{}\n+{}\n{}", i, dna_seq, i, qual_seq).unwrap();
        }

        file
    }
}
