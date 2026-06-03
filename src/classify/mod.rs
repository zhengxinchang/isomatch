use std::{
    fs::File,
    io::{self, BufWriter, Write},
};

use anyhow::Result as AnyResult;
use flate2::{Compression, write::GzEncoder};
use log::info;
use num_format::{Locale, ToFormattedString};

use crate::{
    ClassifyArgs,
    classify::{
        classify_policy::{ClassifyRecord, update_group3_seq_context, update_group4_3rd_party},
        query_ptir::QueryPTIRManager,
        ref_ptir_manager::RefPTIRManager,
    },
    index::fasta::{FaType, FastaReader},
    traits::ArgValidate,
    utils::greetings2,
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

    // open output files
    let mut class_table_table = args.out.clone();
    class_table_table.add_extension("classification.txt.gz");

    let mut class_table_table_writer = File::create(class_table_table)
        .map(|f| BufWriter::new(GzEncoder::new(f, Compression::default())))?;

    class_table_header(&mut class_table_table_writer)?;

    info!("Loading Reference GTF");

    let ref_ptir_manager = RefPTIRManager::open(&args.ref_gtf)?;

    info!("Loading Reference FASTA");
    let mut fa_reader = FastaReader::open(&args.ref_fa, FaType::Ref)?;

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

        update_group4_3rd_party(&mut best_class_record, &query_ptir, None, None, &args);

        best_class_record.write_to_file(&mut class_table_table_writer)?;
    }

    info!("Finished!");
    Ok(())
}

pub fn class_table_header(writer: &mut dyn Write) -> Result<(), io::Error> {
    writeln!(
        writer,
        "isoform_id\tchrom\tstrand\tquery_length\tquery_exon_n\tstructural_category\tsubcategory\tref_gene_id\tref_gene_name\tref_tx_id\tref_length\tref_exon_n\tref_strand\tdiff_to_tss\tdiff_to_tes\tmatched_junctions\tmatched_exons\tbite\tall_canonical\tperc_a_downstream_tts\tseq_a_downstream_tts\twithin_cage_peak\tdist_to_cage_peak\twithin_poly_a_peak\tdist_to_poly_a_site"
    )
}
