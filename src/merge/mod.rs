use std::{fs::File, path};

use crate::{MergeArgs, index::reader::IndexReader, traits::ArgValidate};

pub mod policy;
use anyhow::{Ok, Result};

impl ArgValidate for MergeArgs {
    fn validate(&self) {
        // place holder
    }
}



pub fn run_merge(args:MergeArgs) ->Result<()> {

    // open all files (isomx) into a vec


    let fhs: Vec<IndexReader> = args.inputs.iter().map(|pathb| {
        let f = File::open(pathb)?;
        IndexReader::open(f)
    }).collect::<std::io::Result<Vec<_>>>()?;


    // collect all chromsome from all files and build a unique list

    // get chromsome names and get chromblockreader from all indexxreader, for each chromsome do:

    // k-way merge, build super cluster

    // build the junction-level custer 

    // merge canonical tx

    // non-canoniical tx to canonical tx

    // report to unified GTF

    Ok(())
}
