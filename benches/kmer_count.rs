/* std use */
use std::io::Write;

use needletail::FastxReader;
/* crate use */
use rand::Rng;
use rand::SeedableRng;

/* project use */
use in_place_fastx;
use rayon::iter::ParallelBridge;
use rayon::iter::ParallelIterator;

/* utils function */
fn generate_fastq(seed: u64, nb_seq: usize, length: usize) -> tempfile::NamedTempFile {
    let mut rng = rand::rngs::StdRng::seed_from_u64(seed);

    let mut file = tempfile::NamedTempFile::new().unwrap();

    let dna = [b'A', b'C', b'T', b'G'];
    let qual = (0..94).collect::<Vec<u8>>();

    for i in 0..nb_seq {
        let dna_seq = (0..length)
            .map(|_| dna[rng.gen_range(0..4)] as char)
            .collect::<String>();
        let qual_seq = (0..length)
            .map(|_| (qual[rng.gen_range(0..94)] + 33) as char)
            .collect::<String>();

        writeln!(file, "@{}\n{}\n+{}\n{}", i, dna_seq, i, qual_seq).unwrap();
    }

    file
}

/* Counter definition */
type Counter<T, const N: usize> = [T; N];
trait AbsCounter {
    fn new() -> Self;
}

#[cfg(not(tarpaulin_include))]
impl<const N: usize> AbsCounter for Counter<u64, N> {
    fn new() -> Self {
        [0; N]
    }
}

#[cfg(not(tarpaulin_include))]
impl AbsCounter for Counter<std::sync::atomic::AtomicU64, 1024> {
    fn new() -> Self {
        unsafe { std::mem::transmute([0_u64; 1024]) }
    }
}

const K: u32 = 5;
const KMER_SPACE: usize = 2_usize.pow(K * 2);

/* in_place_fastx parser definition */
in_place_fastx::fastq_sequential!(
    Parser,
    Counter<u64, KMER_SPACE>,
    |record: in_place_fastx::block::Record, counter: &mut Counter<u64, KMER_SPACE>| {
    for kmer in cocktail::tokenizer::Tokenizer::new(record.sequence, K as u8) {
            counter[kmer as usize] += 1;
    }
    }
);

fn in_place_fastx_kmer_count<P>(path: P, block_length: u64) -> Counter<u64, KMER_SPACE>
where
    P: AsRef<std::path::Path>,
{
    let mut counter = Counter::new();
    let mut parser = Parser::new();
    parser
        .with_blocksize(block_length, path, &mut counter)
        .unwrap();

    counter
}

/* in_place_fastx parallel parser definition */
in_place_fastx::fastq_sharedstate!(
    ParserParallel, Counter<std::sync::atomic::AtomicU64, KMER_SPACE>, |record: in_place_fastx::block::Record, counter: &Counter<std::sync::atomic::AtomicU64, KMER_SPACE>| {
    for kmer in cocktail::tokenizer::Tokenizer::new(record.sequence, K as u8) {
        counter[kmer as usize].fetch_add(1, std::sync::atomic::Ordering::SeqCst);
    }
});

fn in_place_fastx_kmer_count_parallel<P>(
    path: P,
    block_length: u64,
) -> Counter<std::sync::atomic::AtomicU64, KMER_SPACE>
where
    P: AsRef<std::path::Path>,
{
    let mut counter = Counter::new();
    let mut parser = ParserParallel::new();

    parser.with_blocksize(block_length, path, &counter).unwrap();

    counter
}

/* rust bio parser */
fn bio_kmer_count<P>(path: P, capacity: u64) -> Counter<u64, KMER_SPACE>
where
    P: AsRef<std::path::Path>,
{
    let mut records = bio::io::fastq::Reader::new(std::io::BufReader::with_capacity(
        capacity as usize,
        std::fs::File::open(path).unwrap(),
    ))
    .records();
    let mut counter = Counter::new();

    while let Some(Ok(record)) = records.next() {
        for kmer in cocktail::tokenizer::Tokenizer::new(record.seq(), K as u8) {
            counter[kmer as usize] += 1;
        }
    }

    counter
}

/* needletail parser */
fn needletail_kmer_count<P>(path: P, capacity: u64) -> Counter<u64, KMER_SPACE>
where
    P: AsRef<std::path::Path>,
{
    let mut reader = needletail::parser::FastqReader::with_capacity(
        std::fs::File::open(path).unwrap(),
        capacity as usize,
    );

    let mut counter = Counter::new();

    while let Some(Ok(record)) = reader.next() {
        for kmer in cocktail::tokenizer::Tokenizer::new(record.raw_seq(), K as u8) {
            counter[kmer as usize] += 1;
        }
    }

    counter
}

fn blocksize(c: &mut criterion::Criterion) {
    let file = generate_fastq(42, 100_000, 150);

    let mut g = c.benchmark_group("blocksize");

    g.throughput(criterion::Throughput::Elements(
        file.path().metadata().unwrap().len(),
    ));

    rayon::ThreadPoolBuilder::new()
        .num_threads(4)
        .build_global()
        .unwrap();

    for power2 in 11..24 {
        g.bench_with_input(
            criterion::BenchmarkId::new("in_place_fastx", 2_u64.pow(power2)),
            &2_u64.pow(power2),
            |b, &length| b.iter(|| criterion::black_box(in_place_fastx_kmer_count(&file, length))),
        );
        g.bench_with_input(
            criterion::BenchmarkId::new("in_place_fastx_parallel", 2_u64.pow(power2)),
            &2_u64.pow(power2),
            |b, &length| {
                b.iter(|| criterion::black_box(in_place_fastx_kmer_count_parallel(&file, length)))
            },
        );
        g.bench_with_input(
            criterion::BenchmarkId::new("bio", 2_u64.pow(power2)),
            &2_u64.pow(power2),
            |b, &length| b.iter(|| criterion::black_box(bio_kmer_count(&file, length))),
        );
        g.bench_with_input(
            criterion::BenchmarkId::new("needletail", 2_u64.pow(power2)),
            &2_u64.pow(power2),
            |b, &length| b.iter(|| criterion::black_box(needletail_kmer_count(&file, length))),
        );
    }
}

fn setup(c: &mut criterion::Criterion) {
    let _ = env_logger::builder().is_test(true).try_init();

    blocksize(c);
}

criterion::criterion_group!(benches, setup);

criterion::criterion_main!(benches);
