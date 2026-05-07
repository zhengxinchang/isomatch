use log::error;
use rustc_hash::FxHashMap;

use crate::{
    MergeArgs,
    core::{ptir::PTIR, status::TxType, tx_strand::ISOMSTRAND},
    merge::{merge_error::MergeError, mptir::MPTIR},
};

pub fn merge_cluster(
    n_exon: u16,
    strand: ISOMSTRAND,
    cluster_idx: &Vec<usize>,
    scluster: &[PTIR],
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

        let mut mptirs = merge_canonical(canonical_vec, scluster, args);

        let rest_ptir_idxs = noncannonical_to_canonical(&mut mptirs, noncanonical_vec, args);

        let rest_mptirs = merge_rest_noncanonical(rest_ptir_idxs, scluster, args);

        mptirs.extend(rest_mptirs.into_iter());
    } else {
        // mono-exon
    }

    todo!()
}

// Merge canonical transcripts into at least one MPTIR, each one is a 
// backbone
// for canonical transcirpt, no wobble are allowed by default.
pub fn merge_canonical(canonical_vec: Vec<&usize>,super_cluster:&[PTIR], args: &MergeArgs) -> Vec<MPTIR> {
    // the ptir has same strand, same exons and all canonical
    // only consider the TSS and TES distance is ok
    // because all canonical should be able to correctly splice

    for anchor_idx in 0..canonical_vec.len() {
        for cmp_idx in anchor_idx+1 .. canonical_vec.len() {
            let anchor = canonical_vec[anchor_idx];
            let cmp  = canonical_vec[anchor_idx];
            let anchor_ptir = &super_cluster[*anchor];
            let cmp_ptir = &super_cluster[*cmp];

            // let anchor_junctions =  anchor_ptir.junction_vec.as_ref();
            // let cmp_junctions = cmp_ptir.junction_vec.as_ref();

            // for exon_idx in 0..(anchor_ptir.n_exons -1) as usize {
            //     let anchor_junc  = &(anchor_ptir.junction_vec.unwrap()[exon_idx]);
            //     let cmp_junc = &cmp_ptir.junction_vec.unwrap()[exon_idx];

            // }

            match (anchor_ptir.junction_vec.as_ref(), cmp_ptir.junction_vec.as_ref()) {
                (Some(a),Some(b)) => {
                    let mut sj_identical = true;
                    for exon_idx in 0..a.len() -1 {
                        
                        if a[exon_idx].0 != b[exon_idx].0 || a[exon_idx].1 != b[exon_idx].1 {
                            sj_identical = false;
                        }
                    }

                 


                }
                _ => {
                    error!("Junction is missing for transcript...");
                    std::process::exit(1);
                }
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

pub fn merge_rest_noncanonical(rest_non_canon: Vec<&usize>, scluster:&[PTIR], args: &MergeArgs) -> Vec<MPTIR> {
    todo!()
}

pub fn merge_mono_exon(mono_vec: Vec<&usize>, scluster:&[PTIR],args: &MergeArgs) -> Vec<MPTIR> {
    todo!()
}

/// find the refhash equivalent in the super cluster
pub fn find_refhash_equivalent(scluster: &[PTIR]) -> FxHashMap<usize, Vec<usize>> {
    todo!()
}
