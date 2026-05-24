use anyhow::{Result as AnyResult, anyhow};
use log::info;

use crate::{
    ClassifyArgs,
    classify::{classify_policy::classify, query_ptir::QueryPTIRManager, ref_ptir::RefPTIRManager},
    index::fasta::{FaType, FastaReader},
    traits::ArgValidate,
};

pub mod class_code;
pub mod classify_error;
pub mod classify_policy;
pub mod query_ptir;
pub mod ref_ptir;

impl ArgValidate for ClassifyArgs {
    fn validate(&self) {
        todo!()
    }
}

pub fn run_classify(args: ClassifyArgs) -> AnyResult<()> {
    args.validate();

    info!("Loading Reference GTF");

    let mut ref_ptir_manager = RefPTIRManager::open(&args.ref_gtf)?;

    info!("Loading Reference FASTA");
    let mut fa_reader = FastaReader::open(&args.ref_fa, FaType::Ref)?;

    info!("Load Query GTF");

    let mut query_ptir_manager = QueryPTIRManager::open(&args.cmp_gtf)?;

    info!("Start classification");

    loop {
        let Some(query_ptir) = query_ptir_manager.next_record() else {
            break;
        };
        let mut ref_candidates =
            ref_ptir_manager.find_ovlp(&query_ptir.chr_name, query_ptir.start(), query_ptir.end());

        let class = match ref_candidates {
            Some(candidates) => {
                let mut classes = Vec::new();
                for candidate in candidates {
                    classes.push(classify(&query_ptir, candidate, None, None, &mut fa_reader))
                }
            }
            None => todo!(),
        };
    }

    Ok(())
}
