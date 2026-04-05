use anyhow::Context;
use log::{error, info};

use crate::{
    IndexArgs,
    fasta::{self, FastaReader},
    gtf::{self, profile_gtf},
    index::format::ChromBlockBuilder,
    traits::ArgValidate,
};
pub use anyhow::Result;
pub mod builder;
pub mod format;
pub mod index_error;
pub mod reader;

impl ArgValidate for IndexArgs {
    fn validate(&self) {
        let mut error_msg = "".to_string();
        let mut has_error = false;

        if !self.input.exists() {
            error_msg.push_str(&format!(
                "\nInput GTF file does not exist: {:?}",
                self.input
            ));
            has_error = true;
        }

        if !self.reffa.exists() {
            error_msg.push_str(&format!(
                "\nReference FASTA file does not exist: {:?}",
                self.reffa
            ));
            has_error = true;
        }

        let mut fai1 = self.reffa.clone();
        fai1.add_extension("fai");
        if !fai1.exists() {
            error_msg.push_str(&format!(
                "\nReference FASTA index file does not exist: {:?}, use ' samtools faidx {} ' to create one.",
                fai1,
                self.reffa.display()
            ));
            has_error = true;
        }

        if let Some(seqfa) = &self.seqfa {
            if !seqfa.exists() {
                error_msg.push_str(&format!(
                    "\nSequence FASTA file does not exist: {:?}",
                    seqfa
                ));
                has_error = true;
            }

            let mut seqfai1 = seqfa.clone();
            seqfai1.add_extension("fai");
            if !seqfai1.exists() {
                error_msg.push_str(&format!(
                    "\nSequence FASTA index for {:?} does not exist, use ' samtools faidx {} ' to create one.",
                    seqfai1,
                    seqfa.display()
                ));
                has_error = true;
            }
        }

        if has_error {
            error!("Error validating arguments: {}", error_msg);
            std::process::exit(1);
        }
    }
}

pub fn run_index(args: &mut IndexArgs) -> Result<()> {
    args.validate();

    info!("Creating isomatch index for {}", args.input.display());

    info!("Loading Reference and/or Sequence FASTA...");

    let mut ref_far = FastaReader::open(args.reffa.clone(), fasta::FaType::Ref)
        .with_context(|| format!("Can not load reference sequence: {}", args.reffa.display()))?;

    let mut seq_far = if let Some(seqfa) = &args.seqfa {
        Some(
            FastaReader::open(seqfa.clone(), fasta::FaType::Seq).with_context(|| {
                format!(
                    "Can not load sequence from reference genome: {}",
                    seqfa.display()
                )
            })?,
        )
    } else {
        None
    };

    info!("Profiling GTF...");
    let gtf_size = std::fs::metadata(&args.input)
        .expect("Can not read gtf file metadata")
        .len();
    let md5 = crate::utils::checksum_file(&args.input)
        .expect("Can not read gtf file for checksum")
        .0;
    let chrom_names = profile_gtf(&args.input).expect("Can not profile GTF file");

    let output_path = if let Some(out) = &args.out {
        out.clone()
    } else {
        let mut default_out = args.input.clone();
        default_out.set_extension("isomx");
        default_out
    };

    info!("Init Builder...");
    let mut builder = builder::IndexBuilder::new(
        std::fs::File::create(&output_path).expect("Can not create output file"),
        chrom_names,
        gtf_size,
        md5,
        true,
        args.seqfa.is_some(),
    )
    .expect("Can not init index builder");

    info!("Loading GTF...");
    let gtf_reader = gtf::MyGTFReader::new(&args.input)
        .with_context(|| format!("Can not open GTF file: {}", args.input.display()))?;

    let mut current_chrom = String::new();
    let mut chrom_id = 0u16;
    let mut chrom_block: Option<ChromBlockBuilder> = None;

    for tx_record in gtf_reader {
        if current_chrom != tx_record.chrom {
            if let Some(cb) = chrom_block.take() {
                builder.add_chrom(cb)?;
            }
            current_chrom = tx_record.chrom.clone();
            chrom_id += 1;
            chrom_block = Some(ChromBlockBuilder::init(chrom_id));
            info!("Processing chromosome {}", &current_chrom);
        }
        chrom_block
            .as_mut()
            .expect("Can not access chromblock")
            .add_tx(tx_record, &mut ref_far, &mut seq_far)?;
    }

    if let Some(cb) = chrom_block.take() {
        builder.add_chrom(cb)?;
    }
    builder.finalize()?;

    info!("Index written to {:?}", output_path);
    Ok(())
}
