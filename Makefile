

build:
	cargo build --release


index: build
	/usr/bin/time -v target/release/isomatch index \
		--reffa test/GRCh38.p14.allChr.fa \
		test/isoseq_transcripts.sorted.filtered_lite.clean.gtf