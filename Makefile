

build:
	cargo fmt && cargo build --release


index: build
	/usr/bin/time -v target/release/isomatch index \
		--reffa test/GRCh38.p14.allChr.fa \
		test/isoseq_transcripts.sorted.filtered_lite.clean.gtf

index2: build
	/usr/bin/time -v target/release/isomatch index \
		--reffa test/GRCh38.p14.allChr.fa \
		test/gencode.v49.basic.annotation.sorted.gtf.gz 




merge: build
	/usr/bin/time -v target/release/isomatch merge \
		-o test/merge.gtf.gz --splice-policy major --tss-wob 400 --tes-wob 400 \
		test/isoseq_transcripts.sorted.filtered_lite.clean.gtf \
# 		test/isoseq_transcripts.sorted.filtered_lite.clean.perturbed.smoke.gtf.gz \
# 		test/isoseq_transcripts.sorted.filtered_lite.clean.perturbed.smoke.gtf.gz  \
# 		test/gencode.v49.basic.annotation.sorted.gtf.gz \

INPUT := test/isoseq_transcripts.sorted.filtered_lite.clean.gtf
N := 100

merge2: build
	/usr/bin/time -v target/release/isomatch merge \
		-o test/merge.gtf \
		$(foreach i,$(shell seq 1 $(N)),$(INPUT))
		