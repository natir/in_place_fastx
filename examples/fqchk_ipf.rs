use std::fmt::Display;

/* crate use */
use clap::Clap;

/* project use */
use in_place_fastx;
use lazy_static::lazy_static;

#[derive(clap::Clap, Debug)]
#[clap(
    name = "fqchk",
    version = "0.1",
    author = "Pierre Marijon <pierre.marijon-ext@aphp.fr>",
    about = "A clone of seqtk fqchk"
)]
struct Command {
    #[clap(short = 'i', long = "input", about = "Fastq input")]
    pub inputs: Vec<String>,

    #[clap(
        short = 'q',
        long = "quality-threshold",
        about = "Quality threshold",
        default_value = "20"
    )]
    pub qual_t: u8,

    #[clap(
        short = 'b',
        long = "blocksize",
        about = "Control default blocksize",
        default_value = "16384"
    )]
    pub blocksize: u64,
}

const NUC_SPACE: usize = 4;
const QUAL_SPACE: usize = 93;
lazy_static! {
    static ref PERR: [f64; QUAL_SPACE] = {
        let mut tmp: [f64; QUAL_SPACE] = [0.0; QUAL_SPACE];
        let mut k = 0;
        while k != QUAL_SPACE {
            if k < 4 {
                tmp[k] = 0.5
            } else {
                tmp[k] = 10_f64.powf(-1.0 * k as f64)
            }
            k += 1;
        }
        tmp
    };
}

#[derive(Default, Clone)]
pub struct PosInfo {
    threshold: u8,
    quality_counts: Vec<u64>,
    bases_count: [u64; NUC_SPACE],
}

impl PosInfo {
    pub fn new(threshold: u8) -> Self {
        Self {
            threshold,
            quality_counts: vec![0; QUAL_SPACE],
            bases_count: [0; NUC_SPACE],
        }
    }
}

impl Display for PosInfo {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut sum = 0.0;
        let mut sum_low = 0.0;
        let mut q_sum = 0.0;
        let mut p_sum: f64 = 0.0;
        for k in 0..QUAL_SPACE {
            sum += self.quality_counts[k] as f64;
            if k < self.threshold as usize {
                sum_low += self.quality_counts[k] as f64;
            }
            q_sum += k as f64 * self.quality_counts[k] as f64;
            p_sum += k as f64 * PERR[k];
        }

        self.quality_counts
            .iter()
            .enumerate()
            .map(|(qual, count)| qual as u64 * *count)
            .sum::<u64>() as f64;

        write!(f, "{}\t", sum)?; // number of bases
        for nuc in 0..NUC_SPACE {
            // nucleotide order is A, C, G and T
            write!(f, "{:.1}\t", 100.0 * self.bases_count[nuc] as f64 / sum)?;
        }

        write!(f, "{:.1}\t", q_sum / sum)?; // average quality
        write!(f, "{:.1}\t", -4.343 * ((p_sum + 1e-6) / (sum + 1e-6)).ln())?; // average error

        write!(f, "{:.1}\t", 100.0 * sum_low / sum)?;
        write!(f, "{:.1}\t", 100.0 * (sum - sum_low) / sum)?;

        Ok(())
    }
}

#[derive(Default)]
pub struct Data {
    pub phred_offset: u8,
    pub threshold: u8,

    pub min_len: u64,
    pub max_len: u64,
    pub average_len: f64,
    pub seq_count: u64,

    pub pos_infos: Vec<PosInfo>,
}

impl Data {
    pub fn new(threshold: u8) -> Self {
        Self {
            phred_offset: 33,
            threshold,

            min_len: u64::MAX,
            max_len: u64::MIN,
            average_len: 0.0,
            seq_count: 0,

            pos_infos: vec![PosInfo::new(threshold); 150],
        }
    }
}

impl Display for Data {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(
            f,
            "min_len: {}; max_len {}; avg_len {:.2}; {} distinct quality values",
            self.min_len, self.max_len, self.average_len, 0
        )?;

        writeln!(f, "POS\t#bases\t%A\t%C\t%G\t%T\tavgQ\terrQ\t%low\t%high")?;

        for (i, pos) in self.pos_infos.iter().enumerate() {
            writeln!(f, "{}\t{}", i, pos)?;
        }

        Ok(())
    }
}

in_place_fastx::fastq_sequential!(
    Sequential,
    Data,
    |record: in_place_fastx::block::Record, data: &mut Data| {
        // Update min and max length
        if data.min_len > record.quality.len() as u64 {
            data.min_len = record.quality.len() as u64;
        }

        if data.max_len < record.quality.len() as u64 {
            data.max_len = record.quality.len() as u64;
        }

        data.seq_count += 1;
        data.average_len = (1.0 / data.seq_count as f64) * record.sequence.len() as f64
            + (1.0 - (1.0 / data.seq_count as f64)) * data.average_len as f64;

        // Update count table length
        if record.quality.len() > data.pos_infos.len() {
            data.pos_infos.extend(vec![
                PosInfo::new(data.threshold);
                record.quality.len() - data.pos_infos.len()
            ]);
        }

        for (indice, (nuc, qual)) in record.sequence.iter().zip(record.quality).enumerate() {
            data.pos_infos[indice].quality_counts[(qual - data.phred_offset) as usize] += 1;
            data.pos_infos[indice].bases_count[(nuc >> 1 & 0b11) as usize] += 1;
        }
    }
);

fn main() -> in_place_fastx::error::Result<()> {
    let params = Command::parse();

    let mut data = Data::new(params.qual_t);
    let mut parser = Sequential::new();
    for input in params.inputs {
        parser.parse(input, &mut data)?;
    }

    println!("{}", data);
    Ok(())
}
