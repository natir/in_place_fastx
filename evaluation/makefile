../target/release/examples/fqchk_ipf: ../examples/fqchk_ipf.rs
	cargo build --release --example fqchk_ipf


../target/release/examples/fqchk_needletail: ../examples/fqchk_needletail.rs
	cargo build --release --example fqchk_needletail


../target/release/examples/fqchk_bio: ../examples/fqchk_bio.rs
	cargo build --release --example fqchk_bio


fqchk.json: ../target/release/examples/fqchk_ipf reads.fastq
	hyperfine --warmup 3 --export-json $@ -n seqtk -n ipf_seq -n needletail -n rust-bio "seqtk fqchk reads.fastq" "../target/release/examples/fqchk_ipf -i reads.fastq" "../target/release/examples/fqchk_needletail -i reads.fastq" "../target/release/examples/fqchk_bio -i reads.fastq"


reads.fastq:
	curl ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR022/ERR022075/ERR022075_1.fastq.gz | seqtk seq - > $@
