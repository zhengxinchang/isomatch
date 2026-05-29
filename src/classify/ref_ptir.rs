use crate::core::{ptir::PTIR, string_pool::StringSpan, tx_strand::ISOMSTRAND, tx_type::TxType};

#[derive(Debug, Clone)]
pub struct RefPTIR {
    pub base: PTIR,
    pub attrs: Vec<(StringSpan, StringSpan)>, // key,value
    pub gene_name: String,
    pub chr_name: String,
}

impl RefPTIR {
    pub fn chr_name(&self) -> &str {
        self.chr_name.as_str()
    }

    pub fn start(&self) -> u32 {
        self.base.start
    }

    pub fn end(&self) -> u32 {
        self.base.end
    }

    pub fn standard(&self) -> &ISOMSTRAND {
        &self.base.strand
    }

    pub fn n_exons(&self) -> u16 {
        self.base.n_exons
    }

    pub fn junction_vec(&self) -> &Option<Vec<(u32, u32)>> {
        &self.base.junction_vec
    }

    pub fn junction_vec_ref(&self) -> &[(u32, u32)] {
        self.base.junction_vec.as_deref().unwrap_or(&[])
    }

    pub fn tx_type(&self) -> &TxType {
        &self.base.tx_type
    }

    pub fn exons_vec(&self) -> Vec<(u32, u32)> {
        let Some(junctions) = self.base.junction_vec.as_ref() else {
            return vec![(self.start(), self.end())];
        };

        if junctions.is_empty() {
            return vec![(self.start(), self.end())];
        }

        let mut boundaries = Vec::with_capacity(junctions.len() * 2 + 2);
        boundaries.push(self.start());
        boundaries.extend(junctions.iter().flat_map(|&(left, right)| [left, right]));
        boundaries.push(self.end());

        boundaries
            .chunks_exact(2)
            .map(|chunk| (chunk[0], chunk[1]))
            .collect()
    }

    pub fn reference_gene_name<'a>(&'a self) -> &'a str {
        if self.gene_name.is_empty() {
            self.base.source_geneid.as_str()
        } else {
            self.gene_name.as_str()
        }
    }

    pub fn transcript_len(&self) -> u32 {
        self.exons_vec()
            .iter()
            .map(|(start, end)| end - start)
            .sum()
    }
}
