const NUC_SPACE: usize = 8;
const QUAL_SPACE: usize = 93;
lazy_static::lazy_static! {
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

#[derive(clap::Clap, Debug)]
#[clap(
    name = "fqchk",
    version = "0.1",
    author = "Pierre Marijon <pierre.marijon-ext@aphp.fr>",
    about = "A clone of seqtk fqchk"
)]
pub struct Command {
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

#[derive(Default, Clone)]
pub struct PosInfo {
    pub threshold: u8,
    pub quality_counts: Vec<u64>,
    pub bases_count: [u64; NUC_SPACE], // A -> 0, C -> 1, G -> 3, T -> 2, N -> 7
}

impl PosInfo {
    pub fn new(threshold: u8) -> Self {
        Self {
            threshold,
            quality_counts: vec![0; QUAL_SPACE],
            bases_count: [0; NUC_SPACE],
        }
    }

    pub fn add(&mut self, rhs: &PosInfo) {
        for i in 0..QUAL_SPACE {
            self.quality_counts[i] += rhs.quality_counts[i];
        }

        for i in 0..NUC_SPACE {
            self.bases_count[i] += rhs.bases_count[i];
        }
    }
}

impl std::fmt::Display for PosInfo {
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
        write!(f, "{:.1}\t", 100.0 * self.bases_count[0] as f64 / sum)?; // A
        write!(f, "{:.1}\t", 100.0 * self.bases_count[1] as f64 / sum)?; // C
        write!(f, "{:.1}\t", 100.0 * self.bases_count[3] as f64 / sum)?; // G
        write!(f, "{:.1}\t", 100.0 * self.bases_count[2] as f64 / sum)?; // T
        write!(f, "{:.1}\t", 100.0 * self.bases_count[7] as f64 / sum)?; // N

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

            pos_infos: Vec::new(),
        }
    }
}

impl std::fmt::Display for Data {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(
            f,
            "min_len: {}; max_len {}; avg_len {:.2}; {} distinct quality values",
            self.min_len, self.max_len, self.average_len, 0
        )?;

        writeln!(
            f,
            "POS\t#bases\t%A\t%C\t%G\t%T\t%N\tavgQ\terrQ\t%low\t%high"
        )?;

        let mut all = PosInfo::new(self.threshold);
        for (i, pos) in self.pos_infos.iter().enumerate() {
            writeln!(f, "{}\t{}", i, pos)?;
            all.add(pos);
        }
        writeln!(f, "ALL\t{}", all)?;

        Ok(())
    }
}
