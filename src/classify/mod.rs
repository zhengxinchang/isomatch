use std::{
    fs::File,
    io::{self, BufWriter, Write},
};

use ahash::HashSet;
use anyhow::{Result as AnyResult, anyhow};
use flate2::{Compression, write::GzEncoder};
use log::info;
use noodles_fasta::record;
use num_format::{Locale, ToFormattedString};

use crate::{
    ClassifyArgs,
    classify::{
        class_code::ClassCode,
        classify_policy::{ClassifyRecord, classify, update_group3, update_group4},
        query_ptir::{QueryPTIR, QueryPTIRManager},
        ref_ptir::{RefPTIR, RefPTIRManager},
    },
    index::fasta::{FaType, FastaReader},
    traits::ArgValidate,
};

pub mod class_code;
pub mod classify_error;
pub mod classify_policy;
pub mod query_ptir;
pub mod ref_ptir;

impl ArgValidate for ClassifyArgs {
    fn validate(&self) {}
}

pub fn run_classify(args: ClassifyArgs) -> AnyResult<()> {
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
    let mut classes = Vec::new();
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
        let ref_candidates =
            ref_ptir_manager.find_ovlp(&query_ptir.chr_name, query_ptir.start(), query_ptir.end());

        let best_class_record = match ref_candidates {
            None => ClassifyRecord::new_intergenic(&query_ptir),
            Some(ref_candidate_vec) => {
                classes.clear();
                for candidate in &ref_candidate_vec {
                    classes.push(classify(&query_ptir, candidate))
                }

                let best = find_best_class(&mut classes);

                let mut best = if matches!(best.cc(), ClassCode::NNC(_))
                    || matches!(best.cc(), ClassCode::NIC(_))
                {
                    if let Some(fusion_class_record) = find_fusion(&classes, &query_ptir) {
                        fusion_class_record
                    } else {
                        best
                    }
                } else {
                    best
                };

                update_group3(&mut best, &query_ptir, &mut fa_reader);
                update_group4(&mut best, &query_ptir, None, None);
                // update_group3(&mut best,&query_ptir,&mut fa_reader);
                best
            }
        };
        best_class_record.write_to_file(&mut class_table_table_writer)?;
    }

    info!("Finished!");
    Ok(())
}

/// find the best class when all candidate reference ptir within the same gene.
pub fn find_best_class(mut classes: &mut Vec<ClassifyRecord>) -> ClassifyRecord {
    let best_rank = classes.iter().map(|c| priority_rank(c.cc())).min().unwrap();
    classes.retain(|c| priority_rank(c.cc()) == best_rank);

    let winner = match classes[0].cc() {
        ClassCode::FSM(_) | ClassCode::ISM(_) => classes
            .into_iter()
            .reduce(|best, c| {
                if c.endpoint_total_diff() < best.endpoint_total_diff() {
                    c
                } else {
                    best
                }
            })
            .unwrap(),
        ClassCode::NIC(_) | ClassCode::NNC(_) => classes
            .into_iter()
            .reduce(|best, c| {
                if c.matched_junctions() > best.matched_junctions() {
                    c
                } else if c.matched_junctions() == best.matched_junctions()
                    && c.exon_count_diff() < best.exon_count_diff()
                {
                    c
                } else {
                    best
                }
            })
            .unwrap(),
        _ => classes.into_iter().next().unwrap(),
    };

    winner.clone()
}

/// check if the query matched with fusion or morejunction.
pub fn find_fusion(
    class_records: &Vec<ClassifyRecord>,
    query_ptir: &QueryPTIR,
) -> Option<ClassifyRecord> {
    let mut geneset = HashSet::default();
    class_records.iter().for_each(|record| {
        if matches!(record.cc(), ClassCode::NIC(_)) || matches!(record.cc(), ClassCode::NNC(_)) {
            geneset.insert(record.ref_gene_name());
        }
    });

    if geneset.len() > 1 {
        Some(ClassifyRecord::new_fusion(query_ptir, &geneset))
    } else {
        None
    }
}

fn priority_rank(cc: ClassCode) -> u8 {
    match cc {
        ClassCode::FSM(_) => 0,
        ClassCode::ISM(_) => 1,
        ClassCode::NIC(_) => 2,
        ClassCode::NNC(_) => 3,
        ClassCode::Fusion => 4,        // determine before any FSM etc
        ClassCode::MoreJunctions => 5, // pigeon extended, same priority as fusion
        ClassCode::Antisense => 6,     // associationOverlapping 第一优先
        ClassCode::GenicGenomic => 7,  // genic in squanti3
        ClassCode::GenicIntron => 8,   // genic but in intron
        ClassCode::Intergenic => 9,    // rest
    }
}

pub fn class_table_header(writer: &mut dyn Write) -> Result<(), io::Error> {
    writeln!(
        writer,
        "isoform_id\tchrom\tstrand\tquery_length\tquery_exon_n\tstructural_category\tsubcategory\tref_gene_id\tref_gene_name\tref_tx_id\tref_length\tref_exon_n\tref_strand\tdiff_to_tss\tdiff_to_tes\tmatched_junctions\tmatched_exons\tbite\tall_canonical\tperc_a_downstream_tts\tseq_a_downstream_tts\twithin_cage_peak\tdist_to_cage_peak\twithin_poly_a_peak\tdist_to_poly_a_site"
    )
}
