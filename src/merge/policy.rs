use rustc_hash::FxHashMap;

use crate::{
    MergeArgs,
    core::{ptir::PTIR, status::TxType},
    merge::{merge_error::MergeError, mptir::MPTIR},
};

pub fn merge_cluster(
    n_exon: u16,
    strand: u8,
    cluster_idx: &Vec<usize>,
    scluster: &Vec<PTIR>,
    args: &MergeArgs,
) -> Result<Vec<MPTIR>, MergeError> {
    if n_exon != 1 {
        // multiple exon tx
        // further splice ptirs into canonical and non-canonical
        let mut canonical_vec = Vec::new();
        let mut noncanonical_vec = Vec::new();
        for ptir_idx in cluster_idx {
            match scluster[*ptir_idx].tx_type {
                TxType::ALLC => {
                    canonical_vec.push(ptir_idx);
                }
                TxType::NOTC | TxType::PRTC => {
                    noncanonical_vec.push(ptir_idx);
                }
                TxType::MONO => {
                    return Err(MergeError::TxType {
                        reason: "Mono exon transcript shouldn't have exon number > 1".to_string(),
                    });
                }
            }
        }
        // merge canonical

        let mut mptirs = merge_canonical(canonical_vec, args);

        let rest_ptir_idxs = noncannonical_to_canonical(&mut mptirs, noncanonical_vec, args);

        let rest_mptirs = merge_rest_noncanonical(rest_ptir_idxs, args);

        mptirs.extend(rest_mptirs.into_iter());
    } else {
        // mono-exon
    }

    todo!()
}

// Merge canonical transcripts into at least one MPTIR, each one is a 
// backbone
pub fn merge_canonical(canonical_vec: Vec<&usize>,scluster:&Vec<PTIR>, args: &MergeArgs) -> Vec<MPTIR> {
    // the ptir has same strand, same exons and all canonical
    // only consider the TSS and TES distance is ok
    // because all canonical should be able to correctly splice

    for anchor_idx in 0..canonical_vec.len() {
        for cmp_idx in anchor_idx+1 .. canonical_vec.len() {
            let anchor = canonical_vec[anchor_idx];
            let cmp  = canonical_vec[anchor_idx];
            let anchor_ptir = &scluster[*anchor];
            let cmp_ptir = &scluster[*cmp];

            for exon_idx in 0..(anchor_ptir.n_exons -1) as usize {
                let anchor_junc  = &(anchor_ptir.junction_vec.unwrap()[exon_idx]);
                let cmp_junc = &cmp_ptir.junction_vec.unwrap()[exon_idx];

            }
        }
    }
    todo!()
}

pub fn noncannonical_to_canonical<'a>(
    mptirs: &mut Vec<MPTIR>,
    noncanonical_vec: Vec<&'a usize>,
    args: &MergeArgs,
) -> Vec<&'a usize> {
    todo!()
}

pub fn merge_rest_noncanonical(rest_non_canon: Vec<&usize>,scluster:&Vec<PTIR>, args: &MergeArgs) -> Vec<MPTIR> {
    todo!()
}

pub fn merge_mono_exon(mono_vec: Vec<&usize>, scluster:&Vec<PTIR>,args: &MergeArgs) -> Vec<MPTIR> {
    todo!()
}

/// find the refhash equivalent in the super cluster
pub fn find_refhash_equivalent(scluster: &Vec<PTIR>) -> FxHashMap<usize, Vec<usize>> {
    todo!()
}
