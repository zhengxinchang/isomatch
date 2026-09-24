use crate::core::core_error::TxBaseError;
use crate::core::junction_pool::JunctionPool;
use crate::core::splice_site_pair::SpliceSitePair;
use crate::core::splice_site_pool::SpliceSitePool;
use crate::core::splice_site_span::SpliceSiteSpan;
use crate::core::string_pool::StringPool;
use crate::core::tx_base::TxBase;
use crate::index::IndexStats;
use crate::index::fasta::FastaReader;
use crate::index::gtf::Transcript;
use crate::index::index_error::IndexError;
use crate::utils;

use libgtf::gtf::Strand;

// // Flag for index status
// // bit 0, sequence from reference genome (0) or tx sequence (1)
// // bit 1, gtf format, 0: plan text, 1: bgzipped
// pub struct Flags {
//     pub bits: u64,
// }

// impl Flags {
//     // bit 0: 0 = gtf format, 1 = bgzipped
//     const GTF_FORMAT_BIT: u64 = 1 << 0;
//     // bit 1: index has ref hash
//     const REF_HASH_BIT: u64 = 1 << 1;
//     // bit 2: index has seq length
//     const SEQ_HASH_BIT: u64 = 1 << 2;

//     pub fn new() -> Self {
//         Self { bits: 0 }
//     }

//     /// true = ref genome sequence hash is valid
//     pub fn set_ref_hash(&mut self, has_ref_hash: bool) {
//         if has_ref_hash {
//             self.bits &= !Self::REF_HASH_BIT; // clear → ref
//         } else {
//             self.bits |= Self::REF_HASH_BIT; // set → tx seq
//         }
//     }

//     /// Returns true if sequence is from reference genome
//     pub fn get_ref_hash(&self) -> bool {
//         self.bits & Self::REF_HASH_BIT == 0
//     }

//     /// true = seq hash is valid
//     pub fn set_seq_hash(&mut self, has_seq_hash: bool) {
//         if has_seq_hash {
//             self.bits &= !Self::SEQ_HASH_BIT; //
//         } else {
//             self.bits |= Self::SEQ_HASH_BIT; //
//         }
//     }

//     /// Returns true if sequence hash is valid
//     /// true = valid, false = invalid
//     pub fn get_seq_hash(&self) -> bool {
//         self.bits & Self::SEQ_HASH_BIT == 0
//     }

//     /// true = gtf format, false = bgzipped
//     pub fn set_gtf_format(&mut self, is_bgzipped: bool) {
//         if is_bgzipped {
//             self.bits |= Self::GTF_FORMAT_BIT; // set → bgzipped
//         } else {
//             self.bits &= !Self::GTF_FORMAT_BIT; // clear → plain text
//         }
//     }

//     /// Returns true if GTF is bgzipped
//     /// true = bgzipped, false = plain text
//     pub fn get_gtf_format(&self) -> bool {
//         self.bits & Self::GTF_FORMAT_BIT != 0
//     }
// }

/// Builder for constructing a single chrom's data block.
pub struct ChromBlockBuilder {
    pub chrom_id: u16,
    pub txs: Vec<TxBase>,
    pub junction_pool: JunctionPool,
    pub splice_site_pool: SpliceSitePool,
    pub string_pool: StringPool,
}

impl ChromBlockBuilder {
    pub fn init(chrom_id: u16) -> Self {
        Self {
            chrom_id,
            txs: Vec::new(),
            junction_pool: JunctionPool::new(),
            splice_site_pool: SpliceSitePool::new(),
            string_pool: StringPool::new(),
        }
    }

    pub fn add_tx(
        &mut self,
        gtf_tx: Transcript,
        chrom_name: &str,
        refr: &mut FastaReader,
        seqr: &mut Option<FastaReader>,
        stats: &mut IndexStats,
    ) -> Result<(), IndexError> {
        // let intron:Vec<u32> = gtf_tx.exons.iter().flat_map(|(e1,e2)|[*e1,*e2]).collect::<Vec<_>>();
        let intron: Vec<u32> = gtf_tx
            .exons
            .windows(2)
            .flat_map(|w| [w[0].1, w[1].0]) // [exon_n.end, exon_{n+1}.start, ...]
            .collect();

        let junction_span =
            self.junction_pool
                .add(&intron)
                .map_err(|e| IndexError::JunctionPoolAdd {
                    id: gtf_tx.tx_id.clone(),
                    reason: e.to_string(),
                })?;

        let tx_id_span =
            self.string_pool
                .add(gtf_tx.tx_id.as_str())
                .map_err(|e| IndexError::StringPoolAdd {
                    id: gtf_tx.tx_id.clone(),
                    reason: e.to_string(),
                })?;

        let gene_id_span = self.string_pool.add(gtf_tx.gene_id.as_str()).map_err(|e| {
            IndexError::StringPoolAdd {
                id: gtf_tx.gene_id.clone(),
                reason: e.to_string(),
            }
        })?;

        let (refhash, splice_site_pairs) = if gtf_tx.exons.len() == 1 {
            // if mono-exon, dont need to do this calculation
            (0, Vec::new())
        } else {
            let reference_seq = refr
                .fetch(
                    chrom_name,
                    gtf_tx.get_0based_start() as usize,
                    gtf_tx.end as usize,
                    true,
                )
                .map_err(|e| IndexError::FetchSeqFailed {
                    reason: format!(
                        "Can not fetch reference sequence for transcript {} on {}:{}-{}: {}",
                        gtf_tx.tx_id, chrom_name, gtf_tx.start, gtf_tx.end, e
                    ),
                })?;

            let mut exon_offsets: Vec<(u32, u32)> = gtf_tx.get_0based_exon_relative_offset();

            let splice_sites_offsets: Vec<(usize, usize, usize, usize)> = exon_offsets
                .windows(2)
                .map(|e| {
                    (
                        e[0].1 as usize,
                        e[0].1 as usize + 2,
                        (e[1].0 as usize).saturating_sub(2),
                        e[1].0 as usize,
                    )
                })
                .collect();
            // println!("{:?}",gtf_tx.exons);
            // println!("{:?}",exon_offsets);
            // println!("{}",reference_seq.len());

            // shift the first and the last position into 3bp close to another side of the exon.
            // only for first and last exon
            // this ensure that transcripts that share same isoform structure (including the small exon shift).
            // will be correctly assinged same hash vlaue.

            let left_exon = exon_offsets
                .first_mut()
                .ok_or_else(|| IndexError::FetchSeqFailed {
                    reason: "Can not obtain the frst exon".to_string(),
                })?;
            if (left_exon.1 - left_exon.0) > 3 {
                left_exon.0 = left_exon.1 - 3;
            }

            let right_exon = exon_offsets
                .last_mut()
                .ok_or_else(|| IndexError::FetchSeqFailed {
                    reason: "Can not obain the last exon".to_string(),
                })?;

            if (right_exon.1 - right_exon.0) > 3 {
                right_exon.1 = right_exon.0 + 3;
            }

            let mut tx_sequence = Vec::new();

            for region in exon_offsets.into_iter() {
                let bases = &reference_seq[region.0 as usize..region.1 as usize];
                tx_sequence.extend_from_slice(&bases);
            }

            let mut splice_site_pairs = Vec::new();
            for (lstart, lend, rstart, rend) in splice_sites_offsets.into_iter() {
                let left = reference_seq[lstart..lend].to_vec();
                let right = reference_seq[rstart..rend].to_vec();
                splice_site_pairs.push((left, right));
            }

            // refhash
            (utils::hash_u8_vec(&tx_sequence), splice_site_pairs)
        };

        // seqhash
        let seqhash = match seqr {
            Some(reader) => {
                if gtf_tx.exons.len() == 1 {
                    // if mono-exon, dont need to do this calculation
                    0
                } else {
                    let sequence = reader.fetch_all(&gtf_tx.tx_id, false)?;

                    let tx_seq_len: usize = gtf_tx
                        .exons
                        .iter()
                        .map(|(s, e)| (*e - *s + 1) as usize)
                        .sum();

                    if sequence.len() != tx_seq_len {
                        dbg!(&gtf_tx);
                        return Err(IndexError::FetchSeqFailed {
                            reason: format!("Actual sequence length ({}) is not equal to GTF derived sequence length ({}). Affected transcript {}",sequence.len(),tx_seq_len,gtf_tx.tx_id).to_string()
                        });
                    }

                    let first_exon = gtf_tx.exons.first().ok_or(IndexError::FetchSeqFailed {
                        reason: "No exons found".to_string(),
                    })?;

                    let last_exon = gtf_tx.exons.last().ok_or(IndexError::FetchSeqFailed {
                        reason: "No exons found".to_string(),
                    })?;

                    let first_exon_len = (first_exon.1 - first_exon.0 + 1) as usize;
                    let last_exon_len = (last_exon.1 - last_exon.0 + 1) as usize;

                    // exons are sorted by genomic position.
                    // For minus strand the RNA 5'-terminal exon is genomically last,
                    // and the RNA 3'-terminal exon is genomically first.
                    let (tss_exon_len, tes_exon_len) = if gtf_tx.strand == Strand::Minus {
                        (last_exon_len, first_exon_len)
                    } else {
                        (first_exon_len, last_exon_len)
                    };

                    // trim from the outer (TSS) side of the 5'-terminal exon
                    let left_trim = if tss_exon_len > 3 {
                        tss_exon_len - 3
                    } else {
                        0
                    };

                    // trim from the outer (TES) side of the 3'-terminal exon
                    let right_trim = if tes_exon_len > 3 {
                        sequence.len() - (tes_exon_len - 3)
                    } else {
                        sequence.len()
                    };

                    let sliced_seq = &sequence[left_trim..right_trim];

                    utils::hash_u8_slice(&sliced_seq)
                }
            }
            None => 0u128,
        };

        let canonical_junction_count = splice_site_pairs
            .iter()
            .map(|(left, right)| SpliceSitePair::pack(left, right, gtf_tx.strand))
            .try_fold(0usize, |count, pair| {
                let pair: SpliceSitePair = pair?;
                Ok::<usize, TxBaseError>(if pair.is_canonical() {
                    count + 1
                } else {
                    count
                })
            })
            .map_err(|e| IndexError::AddGTFTx {
                id: gtf_tx.tx_id.clone(),
                reason: e.to_string(),
            })?;

        // build splice site pairs from intron boundaries
        // let splice_site_pairs: Vec<(&str, &str)> = Vec::new(); // TODO: extract donor/acceptor dinucleotides from reference
        let splice_site_span = if splice_site_pairs.is_empty() {
            SpliceSiteSpan {
                offset: 0,
                count: 0,
            }
        } else {
            self.splice_site_pool
                .add_pairs(&splice_site_pairs, gtf_tx.strand)
                .map_err(|e| IndexError::AddGTFTx {
                    id: gtf_tx.tx_id.clone(),
                    reason: e.to_string(),
                })?
        };

        let tx_base = TxBase::new(
            gtf_tx.gidx,
            self.chrom_id,
            gtf_tx.start,
            gtf_tx.end,
            gtf_tx.strand,
            seqhash,
            refhash,
            gtf_tx.exons.len() as u16,
            splice_site_span,
            junction_span,
            tx_id_span,
            gene_id_span,
        )
        .map_err(|e| IndexError::AddGTFTx {
            id: gtf_tx.tx_id.clone(),
            reason: e.to_string(),
        })?;

        stats.observe_tx(
            gtf_tx.strand,
            gtf_tx.exons.len(),
            canonical_junction_count,
            &gtf_tx.gene_id,
        );
        self.txs.push(tx_base);

        Ok(())
    }

    // pub fn finalize(&mut self) {
    //     // self.txs.sort_unstable();
    // }
}
