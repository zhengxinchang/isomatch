

build:
	cargo fmt && cargo build --release


index2: build
	/usr/bin/time -v target/release/isomatch index \
		--ref-fa test/GRCh38.p14.allChr.fa \
		test/isoseq_transcripts.sorted.filtered_lite.clean.gtf

index: build
	/usr/bin/time -v target/release/isomatch index -q \
		--ref-fa test/GRCh38.p14.allChr.fa \
		test/gencode.v49.basic.annotation.sorted.gtf.gz 

index3: build
	/usr/bin/time -v target/release/isomatch index \
		--ref-fa test/hg38.fa  --skip-missing-ref-chr \
		test/hg002_ont_drna.isoquant.gtf.gz 


merge: build
	/usr/bin/time -v target/release/isomatch merge \
		-o test/merge.gtf.gz --terminal-refine both \
		--guide-tss /ssd2/projects/isomatch-dev/evidence/human.grch38.tss.bed \
		--guide-tes /ssd2/projects/isomatch-dev/evidence/human.grch38.tes.bed \
		test/isoseq_transcripts.sorted.filtered_lite.clean.gtf 
# 		test/isoseq_transcripts.sorted.filtered_lite.clean.perturbed.smoke.gtf.gz \
# 		test/isoseq_transcripts.sorted.filtered_lite.clean.perturbed.smoke.gtf.gz  \
# 		test/gencode.v49.basic.annotation.sorted.gtf.gz \

INPUT := test/isoseq_transcripts.sorted.filtered_lite.clean.gtf
N := 1000

merge_1k: build
	/usr/bin/time -v target/release/isomatch merge \
		-o test/merge.large.gtf.gz \
		$(foreach i,$(shell seq 1 $(N)),$(INPUT))

merge2gffcompare: build
	/usr/bin/time -v test/gffcompare \
		-o test/gffcompare.large-merge.gtf \
		$(foreach i,$(shell seq 1 $(N)),$(INPUT))

	
merge3: build
	/usr/bin/time -v target/release/isomatch merge \
		-t none \
		-o test/drna.merge.gtf.gz --splice-policy major -d 3 -a 3 -s 200 -e 200 \
		test/hg002_ont_drna.isoquant.gtf.gz \

make_tss_tes_plot_drna:
	python3 stage/sonia-test/plot_merge_tss_tes.py   test/drna.merge.gtf.gz.merged.gtf.gz  --min-count 2   --top 50   -o test/test/drna.merge.gtf.gz.merged_tss_tes
	

merge_no_ends_cmp: build
	/usr/bin/time -v target/release/isomatch merge \
		-o test/merge.noends -t none -d 3 -a 3 -s 200 -e 200 \
		test/isoseq_transcripts.sorted.filtered_lite.clean.gtf

make_tss_tes_plot:
	python3 stage/sonia-test/plot_merge_tss_tes.py   test/merge.noends.merged.gtf.gz     --min-count 2   --top 50   -o test/test/merge.gtf.gz.merged_tss_tes


merge_no_ends_cmp_gencode: build
	/usr/bin/time -v target/release/isomatch merge \
		-o test/gencode.merge.noends -t none -d 3 -a 3 -s 200 -e 200 \
		test/gencode.v49.basic.annotation.sorted.gtf.gz

make_tss_tes_plot_gencode:
	python3 stage/sonia-test/plot_merge_tss_tes.py   test/gencode.merge.noends.merged.gtf.gz     --min-count 2   --top 50   -o test/test/gencode.merge.gtf.gz.merged_tss_tes


merge4 : build
	/usr/bin/time -v target/release/isomatch merge \
		-o test/merge-single.gtf.gz  -s 0 -e 0 -u 0 -S 0 -E 0 -U 0 --mono-ovlp 1.0\
		test/isoseq_transcripts.sorted.filtered_lite.clean.gtf

classify1index:build
	/usr/bin/time -v target/release/isomatch index \
		--ref-fa test/GRCh38.p14.allChr.fa \
		test/isoseq_transcripts.sorted.filtered_lite.clean.gtf 

classify1:build 

	/usr/bin/time -v target/release/isomatch classify \
		-s test/GRCh38.p14.allChr.fa \
		-r test/gencode.v49.basic.annotation.sorted.gtf.gz \
		-o test/classify.test test/isoseq_transcripts.sorted.filtered_lite.clean.gtf