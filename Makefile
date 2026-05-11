

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

index3: build
	/usr/bin/time -v target/release/isomatch index \
		--reffa test/hg38.fa  --skip-missing-ref-seqids \
		test/hg002_ont_drna.isoquant.gtf.gz 


merge: build
	/usr/bin/time -v target/release/isomatch merge \
		-o test/merge.gtf.gz -d 3 -a 3 -s 200 -e 200 \
		--guide-tss /ssd2/projects/isomatch-dev/evidence/human.guide.tss.bed \
		--guide-tes /ssd2/projects/isomatch-dev/evidence/human.guide.tes.bed \
		test/isoseq_transcripts.sorted.filtered_lite.clean.gtf \
# 		test/isoseq_transcripts.sorted.filtered_lite.clean.perturbed.smoke.gtf.gz \
# 		test/isoseq_transcripts.sorted.filtered_lite.clean.perturbed.smoke.gtf.gz  \
# 		test/gencode.v49.basic.annotation.sorted.gtf.gz \

INPUT := test/isoseq_transcripts.sorted.filtered_lite.clean.gtf
N := 1000

merge2: build
	/usr/bin/time -v target/release/isomatch merge \
		-o test/merge.gtf \
		$(foreach i,$(shell seq 1 $(N)),$(INPUT))
	
merge3: build
	/usr/bin/time -v target/release/isomatch merge \
		-o test/merge.gtf.gz --splice-policy major -d 3 -a 3 -s 200 -e 200 \
		test/hg002_ont_drna.isoquant.gtf.gz \
		test/isoseq_transcripts.sorted.filtered_lite.clean.gtf