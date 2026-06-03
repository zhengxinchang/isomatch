use std::{
    fs::File,
    io::{self, BufWriter, Write},
};

use anyhow::Result as AnyResult;
use flate2::{Compression, write::GzEncoder};
use log::{error, info};
use num_format::{Locale, ToFormattedString};

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
    utils::{check_index_ready, greetings2},
};

pub mod class_code;
pub mod classify_error;
pub mod classify_policy;
pub mod compare;
pub mod query_ptir;
pub mod ref_ptir;
pub mod ref_ptir_manager;

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
            info!("Re-indexing {}", gtf.display());
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

        best_class_record.write_to_file(&mut class_table_writer)?;
    }

    info!(
        "classification file saved to: {}",
        class_table_path.display()
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
