//! merged ptir
//! one object for one GTF transcript

use crate::core::ptir::PTIR;

/// MergedPTIR is the representation of a merged transcript
/// it has the prepresentive
pub struct MPTIR {
    repr_start: u32,
    repr_end: u32,
    strand: u8, // 0== plus, 1 == minus
    n_exon: u16,
    junctions: Vec<(u32, u32)>,
    ptir_count: u32, // how many ptir been merged
    ptir_idx: Vec<usize>,
    terminal_vec:Vec<(u32,u32)>
}

impl MPTIR {

    pub fn from_ptir(ptir:&PTIR,ptir_idx:usize) -> MPTIR {
        Self {
            repr_start:0,
            repr_end:0,
            strand: ptir.strand,
            n_exon: ptir.n_exons,
            junctions:ptir.junction_vec.clone().unwrap(),
            ptir_count:1,
            ptir_idx:vec![ptir_idx],
            terminal_vec:vec![(ptir.start,ptir.end)]
        }
    }

    pub fn merge_ptir(&mut self, other_ptir:&PTIR) {

        

    }

    pub fn is_same_junctions(&self, ptir2:&PTIR) -> bool {

    }

    pub fn finalize (&mut self) {

    }
    
}