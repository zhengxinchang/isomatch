use std::{
    collections::{HashMap, HashSet},
    fs::File,
    path::PathBuf,
};

use crate::{
    error::{self, Error},
    fasta::{FaType, FastaReader},
    gtf::{ChromID, GtfReader, Strand},
    index::{AttrIndexBuilder, ChromBlockBuilder, IndexBuilder},
};

#[derive(Debug, Clone, Default)]
pub struct BuildReport {
    pub transcript_count: u64,
    pub gene_count: u64,
    pub skipped_transcript_count: u64,
    pub skipped_gene_count: u64,
    pub missing_seqid_count: u64,
    pub missing_seqids: Vec<String>,

    pub plus_strand_transcript_count: u64,
    pub minus_strand_transcript_count: u64,
    pub unknown_strand_transcript_count: u64,

    pub mono_exon_transcript_count: u64,
    pub multi_exon_transcript_count: u64,
    pub all_canonical_transcript_count: u64,
    pub partial_canonical_transcript_count: u64,
    pub non_canonical_transcript_count: u64,

    pub junction_count: u64,
    pub canonical_junction_count: u64,
    pub non_canonical_junction_count: u64,
    pub canonical_junction_ratio: f64,
}

#[derive(Debug, Clone)]
pub struct BuildConfig {
    pub gtf_path: PathBuf,
    pub reference_fasta: PathBuf,
    pub transcript_fasta: Option<PathBuf>,
    pub isomx_output: PathBuf,
    pub isoms_output: PathBuf,
    pub skip_missing_reference_seqids: bool,
    pub temp_dir: Option<PathBuf>,
}

#[non_exhaustive]
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum BuildEvent {
    LoadingFasta,
    IndexingGtf,
    InitializingBuilders,
    MissingReferenceSequence { seqid: String },
    ProcessingChromosome { name: String },
    SkippingChromosome { name: String },
    Finalizing,
}

pub fn build_index(config: &BuildConfig) -> error::Result<BuildReport> {
    build_index_with_events(config, |_| {})
}

pub fn build_index_with_events<F>(config: &BuildConfig, mut emit: F) -> error::Result<BuildReport>
where
    F: FnMut(BuildEvent),
{
    emit(BuildEvent::LoadingFasta);
    let mut reference = FastaReader::open(&config.reference_fasta, FaType::Ref)?;
    let mut transcript_sequences = match &config.transcript_fasta {
        Some(path) => Some(FastaReader::open(path, FaType::Seq)?),
        None => None,
    };

    emit(BuildEvent::IndexingGtf);
    let mut gtf_reader = match &config.temp_dir {
        Some(temp_dir) => GtfReader::new_in(&config.gtf_path, temp_dir)?,
        None => GtfReader::new(&config.gtf_path)?,
    };
    let profile = gtf_reader.profile().clone();

    let mut missing_seqids: Vec<String> = profile
        .chrom_names
        .iter()
        .filter(|chrom| !reference.contains(chrom))
        .cloned()
        .collect();
    missing_seqids.sort();

    if !missing_seqids.is_empty() {
        if !config.skip_missing_reference_seqids {
            return Err(Error::MissingReferenceSequences {
                seqids: missing_seqids,
            });
        }

        for seqid in &missing_seqids {
            emit(BuildEvent::MissingReferenceSequence {
                seqid: seqid.clone(),
            });
        }
    }

    let missing_seqid_set: HashSet<String> = missing_seqids.iter().cloned().collect();
    let missing_chrom_ids: HashSet<ChromID> = missing_seqid_set
        .iter()
        .filter_map(|chrom| profile.chrom_name_to_id.get(chrom).copied())
        .collect();
    let chrom_names: Vec<String> = profile
        .chrom_names
        .iter()
        .filter(|chrom| !missing_seqid_set.contains(*chrom))
        .cloned()
        .collect();
    let output_chrom_ids: HashMap<ChromID, u16> = profile
        .chrom_names
        .iter()
        .filter(|chrom| !missing_seqid_set.contains(*chrom))
        .enumerate()
        .filter_map(|(output_idx, chrom)| {
            profile
                .chrom_name_to_id
                .get(chrom)
                .map(|profile_id| (*profile_id, (output_idx + 1) as u16))
        })
        .collect();

    if chrom_names.is_empty() {
        return Err(Error::NoIndexableSequences);
    }

    emit(BuildEvent::InitializingBuilders);
    let isomx_file = File::create(&config.isomx_output)?;
    let mut index_builder = IndexBuilder::new(
        isomx_file,
        chrom_names,
        profile.file_size,
        profile.md5,
        true,
        config.transcript_fasta.is_some(),
        missing_seqids.clone(),
    )?;

    let total_indexable_transcripts = gtf_reader.transcript_count_excluding(&missing_chrom_ids);
    let isoms_file = File::create(&config.isoms_output)?;
    let mut attribute_builder =
        AttrIndexBuilder::new(isoms_file, total_indexable_transcripts, &profile.md5)?;

    let mut report = BuildReport {
        missing_seqid_count: missing_seqids.len() as u64,
        missing_seqids,
        ..BuildReport::default()
    };
    let mut gene_ids = HashSet::new();
    let mut skipped_gene_ids = HashSet::new();
    let mut current_chrom_id = 0u16;
    let mut chrom_block: Option<ChromBlockBuilder> = None;
    let mut next_transcript_idx = 0u64;

    while let Some(mut transcript) = gtf_reader.next()? {
        let chrom_name = gtf_reader
            .chrom_name(transcript.chrom_id)
            .ok_or_else(|| Error::IndexBuild {
                reason: format!("invalid chromosome id {}", transcript.chrom_id),
            })?
            .to_string();

        if current_chrom_id != transcript.chrom_id {
            if let Some(block) = chrom_block.take() {
                index_builder.add_chrom(block)?;
            }

            current_chrom_id = transcript.chrom_id;
            if missing_chrom_ids.contains(&current_chrom_id) {
                emit(BuildEvent::SkippingChromosome {
                    name: chrom_name.clone(),
                });
                chrom_block = None;
            } else {
                emit(BuildEvent::ProcessingChromosome {
                    name: chrom_name.clone(),
                });
                let output_chrom_id = output_chrom_ids
                    .get(&current_chrom_id)
                    .copied()
                    .ok_or_else(|| Error::IndexBuild {
                        reason: format!("missing output chromosome id for {chrom_name}"),
                    })?;
                chrom_block = Some(ChromBlockBuilder::init(output_chrom_id));
            }
        }

        if missing_chrom_ids.contains(&transcript.chrom_id) {
            report.skipped_transcript_count += 1;
            skipped_gene_ids.insert(transcript.gene_id);
            continue;
        }

        transcript.set_gidx(next_transcript_idx);
        let attribute = transcript.attr_string.clone();
        let gene_id = transcript.gene_id.clone();
        let summary = chrom_block
            .as_mut()
            .ok_or_else(|| Error::IndexBuild {
                reason: "cannot access chromosome block".to_string(),
            })?
            .add_tx(
                transcript,
                &chrom_name,
                &mut reference,
                &mut transcript_sequences,
            )?;

        report.transcript_count += 1;
        gene_ids.insert(gene_id);
        match summary.strand {
            Strand::Plus => report.plus_strand_transcript_count += 1,
            Strand::Minus => report.minus_strand_transcript_count += 1,
            Strand::Unknown => report.unknown_strand_transcript_count += 1,
        }

        if summary.exon_count <= 1 {
            report.mono_exon_transcript_count += 1;
        } else {
            report.multi_exon_transcript_count += 1;
            let junction_count = (summary.exon_count - 1) as u64;
            let canonical_count = summary.canonical_junction_count as u64;

            if canonical_count == junction_count {
                report.all_canonical_transcript_count += 1;
            } else if canonical_count == 0 {
                report.non_canonical_transcript_count += 1;
            } else {
                report.partial_canonical_transcript_count += 1;
            }

            report.junction_count += junction_count;
            report.canonical_junction_count += canonical_count;
            report.non_canonical_junction_count += junction_count - canonical_count;
        }

        if let Some(attribute) = attribute {
            attribute_builder.dump_attr(attribute, next_transcript_idx)?;
        }

        next_transcript_idx =
            next_transcript_idx
                .checked_add(1)
                .ok_or_else(|| Error::IndexBuild {
                    reason: "written transcript index exceeded u64".to_string(),
                })?;
    }

    if let Some(block) = chrom_block.take() {
        index_builder.add_chrom(block)?;
    }

    emit(BuildEvent::Finalizing);
    index_builder.finalize()?;
    attribute_builder.finish()?;

    report.gene_count = gene_ids.len() as u64;
    report.skipped_gene_count = skipped_gene_ids.len() as u64;
    report.canonical_junction_ratio = if report.junction_count == 0 {
        0.0
    } else {
        report.canonical_junction_count as f64 / report.junction_count as f64
    };

    Ok(report)
}
