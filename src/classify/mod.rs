use std::{
    collections::BTreeMap,
    fs::File,
    io::{self, BufWriter, Write},
};

use anyhow::Result as AnyResult;
use flate2::{Compression, write::GzEncoder};
use log::info;
use num_format::{Locale, ToFormattedString};
use serde::Serialize;

use crate::{
    ClassifyArgs, IndexArgs,
    classify::{
        classify_policy::{ClassifyRecord, update_group3_seq_context, update_group4_regions},
        query_ptir::QueryPTIRManager,
        ref_ptir_manager::RefPTIRManager,
    },
    index::{
        fasta::{FaType, FastaReader},
        run_index,
    },
    merge::guide::GuideDb,
    traits::ArgValidate,
    utils::{check_index_ready, greetings2, print_json_block},
};

pub mod class_code;
pub mod classify_error;
pub mod classify_policy;
pub mod compare;
pub mod query_ptir;
pub mod ref_ptir;
pub mod ref_ptir_manager;

#[derive(Debug, Default, Serialize)]
pub struct ClassifyStats {
    pub reference_tx_cnt: u64,
    pub query_tx_cnt: u64,
    pub classified_tx_cnt: u64,
    pub mono_exon_tx_cnt: u64,
    pub multi_exon_tx_cnt: u64,
    pub all_canonical_tx_cnt: u64,
    pub cage_supported_tx_cnt: u64,
    pub cage_supported_pct: f64,
    pub polya_supported_tx_cnt: u64,
    pub polya_supported_pct: f64,
    pub structural_category_cnt: BTreeMap<String, u64>,
    pub structural_category_pct: BTreeMap<String, f64>,
    pub subcategory_counts: BTreeMap<String, BTreeMap<String, u64>>,
}

impl ClassifyStats {
    pub fn new(query_tx_count: usize, reference_tx_count: usize) -> Self {
        Self {
            query_tx_cnt: query_tx_count as u64,
            reference_tx_cnt: reference_tx_count as u64,
            ..Self::default()
        }
    }

    pub fn observe_record(
        &mut self,
        main_category: &str,
        subcategory: &str,
        exon_count: u16,
        all_canonical: bool,
        within_cage_peak: Option<bool>,
        within_poly_a_peak: Option<bool>,
    ) {
        self.classified_tx_cnt += 1;

        if exon_count <= 1 {
            self.mono_exon_tx_cnt += 1;
        } else {
            self.multi_exon_tx_cnt += 1;
        }

        if all_canonical {
            self.all_canonical_tx_cnt += 1;
        }
        if matches!(within_cage_peak, Some(true)) {
            self.cage_supported_tx_cnt += 1;
        }
        if matches!(within_poly_a_peak, Some(true)) {
            self.polya_supported_tx_cnt += 1;
        }

        *self
            .structural_category_cnt
            .entry(main_category.to_string())
            .or_insert(0) += 1;
        *self
            .subcategory_counts
            .entry(main_category.to_string())
            .or_default()
            .entry(subcategory.to_string())
            .or_insert(0) += 1;
    }

    pub fn finalize(&mut self) {
        if self.classified_tx_cnt == 0 {
            self.cage_supported_pct = 0.0;
            self.polya_supported_pct = 0.0;
            return;
        }

        self.cage_supported_pct =
            (self.cage_supported_tx_cnt as f64 * 100.0 / self.classified_tx_cnt as f64 * 10000.0)
                .round()
                / 10000.0;
        self.polya_supported_pct =
            (self.polya_supported_tx_cnt as f64 * 100.0 / self.classified_tx_cnt as f64 * 10000.0)
                .round()
                / 10000.0;
        self.structural_category_pct = self
            .structural_category_cnt
            .iter()
            .map(|(category, count)| {
                (
                    category.clone(),
                    (*count as f64 * 100.0 / self.classified_tx_cnt as f64 * 10000.0).round()
                        / 10000.0,
                )
            })
            .collect();
    }
}

impl ArgValidate for ClassifyArgs {
    fn validate(&self) {}
}

pub fn run_classify(args: ClassifyArgs) -> AnyResult<()> {
    greetings2(&args);
    args.validate();

    info!("Loading 1 gtf(s)");
    // auto indexing
    info!(
        "Checking input GTF indexes; missing/corrupted/outdated indexes will be created automatically."
    );
    for gtf in [&args.input, &args.ref_gtf] {
        if !check_index_ready(gtf) {
            info!("Indexing {}", gtf.display());
            let mut index_args = IndexArgs {
                input: gtf.clone(),
                ref_fa: args.ref_fa.clone(),
                seqfa: None,
                out: None,
                skip_missing_ref_chr: args.skip_missing_ref_chr,
                quiet: true,
            };
            run_index(&mut index_args)?;
        }
    }

    // open output files
    let mut class_table_path = args.out.clone();
    class_table_path.add_extension("classification.txt.gz");

    let mut class_table_writer = File::create(&class_table_path)
        .map(|f| BufWriter::new(GzEncoder::new(f, Compression::default())))?;

    class_table_header(&mut class_table_writer)?;

    let mut gtf_path = args.out.clone();
    gtf_path.add_extension("annotated.gtf.gz");

    let mut gtf_writer = File::create(&gtf_path)
        .map(|f| BufWriter::new(GzEncoder::new(f, Compression::default())))?;

    update_gtf_header(&mut gtf_writer, &args)?;

    let mut classify_info_path = args.out.clone();
    classify_info_path.add_extension("classify_info.json");

    info!("Loading Reference GTF");

    let ref_ptir_manager = RefPTIRManager::open(&args.ref_gtf)?;

    info!("Loading Reference FASTA");
    let mut fa_reader = FastaReader::open(&args.ref_fa, FaType::Ref)?;

    let reftss_db = match &args.guide_tss {
        Some(path) => Some(GuideDb::from_bed_path(
            path,
            crate::merge::guide::GuideBEDType::Tss,
            &args.chrmap.as_ref(),
        )?),
        None => None,
    };

    let reftes_db = match &args.guide_tes {
        Some(path) => Some(GuideDb::from_bed_path(
            path,
            crate::merge::guide::GuideBEDType::Tes,
            &args.chrmap.as_ref(),
        )?),
        None => None,
    };

    info!("Loading Query GTF");

    let mut query_ptir_manager = QueryPTIRManager::open(&args.input)?;
    let mut stats = ClassifyStats::new(
        query_ptir_manager.total_tx_n(),
        ref_ptir_manager.ptirs.len(),
    );

    info!(
        "{} query transcripts loaded",
        query_ptir_manager
            .total_tx_n()
            .to_formatted_string(&Locale::en)
    );

    info!("Start classification");
    let mut processed_tx = 0u64;
    loop {
        {
            processed_tx += 1;
            if processed_tx % 10000 == 0 {
                info!(
                    "Classified {} transcripts",
                    processed_tx.to_formatted_string(&Locale::en)
                )
            }
        }

        let Some(query_ptir) = query_ptir_manager.next_record() else {
            break;
        };

        let mut best_class_record = ClassifyRecord::new(&query_ptir, &ref_ptir_manager, &args);

        update_group3_seq_context(&mut best_class_record, &query_ptir, &mut fa_reader, &args);

        update_group4_regions(
            &mut best_class_record,
            &query_ptir,
            &reftss_db,
            &reftes_db,
            &args,
        );

        best_class_record.observe_stats(&mut stats);
        best_class_record.write_to_file(&mut class_table_writer, &mut gtf_writer)?;
    }

    stats.finalize();
    class_table_writer.flush()?;
    gtf_writer.flush()?;
    print_json_block("Classify summary", &stats);

    let mut classify_info_writer = File::create(&classify_info_path)?;
    let classify_info_json = serde_json::to_string_pretty(&stats)?;
    classify_info_writer.write(classify_info_json.as_bytes())?;
    classify_info_writer.flush()?;

    info!(
        "Classification file has been saved to: {}",
        class_table_path.display()
    );
    info!("Annotated GTF has been saved to: {}", gtf_path.display());
    info!(
        "Classify summary has been saved to: {}",
        classify_info_path.display()
    );

    info!("Finished!");
    Ok(())
}

pub fn class_table_header(writer: &mut dyn Write) -> Result<(), io::Error> {
    writeln!(
        writer,
        "isoform_id\tchrom\tstrand\tquery_length\tquery_exon_n\tstructural_category\tsubcategory\tref_gene_id\tref_gene_name\tref_tx_id\tref_length\tref_exon_n\tref_strand\tdiff_to_tss\tdiff_to_tes\tmatched_junctions\tmatched_exons\tbite\tall_canonical\tperc_a_downstream_tts\tseq_a_downstream_tts\twithin_cage_peak\tdist_to_cage_peak\twithin_poly_a_peak\tdist_to_poly_a_site"
    )
}

pub fn update_gtf_header(gtf_writer: &mut dyn Write, args: &ClassifyArgs) -> Result<(), io::Error> {
    // read args.input and get headers starts with #
    let escape = |value: &str| value.replace('\\', "\\\\").replace('"', "\\\"");
    let mut has_isom_version = false;
    let mut has_ref_tx_id = false;
    let mut has_ref_gene_id = false;
    let mut has_ref_gene_name = false;
    let mut has_category = false;
    let mut has_subcategory = false;

    let mut reader = crate::utils::open_file_bufread(&args.input)?;
    let mut line = String::new();
    loop {
        line.clear();
        if reader.read_line(&mut line)? == 0 {
            break;
        }
        if !line.starts_with('#') {
            break;
        }

        has_isom_version |= line.starts_with("##ISOM <VERSION>");
        has_ref_tx_id |= line.contains("ID=\"ISOM_REF_TX_ID\"");
        has_ref_gene_id |= line.contains("ID=\"ISOM_REF_GENE_ID\"");
        has_ref_gene_name |= line.contains("ID=\"ISOM_REF_GENE_NAME\"");
        has_category |= line.contains("ID=\"ISOM_CATEGORY\"");
        has_subcategory |= line.contains("ID=\"ISOM_SUBCATEGORY\"");

        gtf_writer.write_all(line.as_bytes())?;
        if !line.ends_with('\n') {
            gtf_writer.write_all(b"\n")?;
        }
    }

    // check if ISOM version exists if not, add
    if !has_isom_version {
        writeln!(
            gtf_writer,
            "##ISOM <VERSION> version=\"{}\"; program=\"{}\"; schema=\"{}\"",
            env!("CARGO_PKG_VERSION"),
            env!("CARGO_PKG_NAME"),
            crate::constants::ISOM_GTF_SCHEMA
        )?;
    }

    // check if headers for ISOM_REF_TX_ID, ISOM_REF_GENE_ID, ISOM_REF_GENE_NAME, ISOM_CATEGORY ISOM_SUBCATEGORY, exists, if not , add.
    if !has_ref_tx_id {
        writeln!(
            gtf_writer,
            "##ISOM <FORMAT> ID=\"ISOM_REF_TX_ID\"; Description=\"best matching reference transcript ID for this classified transcript\";"
        )?;
    }
    if !has_ref_gene_id {
        writeln!(
            gtf_writer,
            "##ISOM <FORMAT> ID=\"ISOM_REF_GENE_ID\"; Description=\"reference gene ID associated with this classified transcript\";"
        )?;
    }
    if !has_ref_gene_name {
        writeln!(
            gtf_writer,
            "##ISOM <FORMAT> ID=\"ISOM_REF_GENE_NAME\"; Description=\"reference gene name associated with this classified transcript\";"
        )?;
    }
    if !has_category {
        writeln!(
            gtf_writer,
            "##ISOM <FORMAT> ID=\"ISOM_CATEGORY\"; Description=\"SQANTI3-style structural category assigned by isomatch classify\";"
        )?;
    }
    if !has_subcategory {
        writeln!(
            gtf_writer,
            "##ISOM <FORMAT> ID=\"ISOM_SUBCATEGORY\"; Description=\"SQANTI3-style structural subcategory assigned by isomatch classify\";"
        )?;
    }

    // append a line about the classify command
    let command = std::env::args_os()
        .map(|arg| escape(&arg.to_string_lossy()))
        .collect::<Vec<_>>()
        .join(" ");
    writeln!(gtf_writer, "##ISOM <COMMAND> cmd={}", command)?;

    Ok(())
}
