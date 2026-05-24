use std::fmt;

/// Structural classification of a query transcript against a reference.
///
/// Variants map directly to SQANTI3 classification file values.
/// Priority order: FSM > ISM > NIC > NNC > Fusion > Antisense
///                 > GenicIntron > GenicGenomic > Intergenic
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum ClassCode {
    /// All internal junctions match; same exon count
    FSM(SubFSM),
    /// Consecutive junction subset matches (fewer exons)
    ISM(SubISM),
    /// All donor/acceptor sites known, novel combination
    NIC(SubNIC),
    /// ≥1 novel donor or acceptor site
    NNC(SubNNC),
    /// Spans multiple annotated genes
    Fusion,
    /// Antisense to an annotated gene
    Antisense,
    /// Completely contained within a reference intron
    GenicIntron,
    /// Overlaps both intronic and exonic regions (same strand)
    GenicGenomic,
    /// No overlap with any annotated gene
    Intergenic,

    // ── Pigeon extensions (not in standard SQANTI3) ──
    /// Like Fusion but ≥1 junction shared between multiple genes
    /// (pigeon-specific, not in SQANTI3)
    MoreJunctions,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum SubFSM {
    /// Both ends within 50bp of reference TSS/TTS
    ReferenceMatch,
    /// 3' end differs >50bp
    Alternative3End,
    /// 5' end differs >50bp
    Alternative5End,
    /// Both ends differ >50bp
    Alternative5And3End,
    /// Single-exon query overlapping a mono-exonic reference
    MonoExon,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum SubISM {
    /// Missing junctions at 5' → retains 3' portion
    ThreePrimeFragment,
    /// Missing junctions at 3' → retains 5' portion
    FivePrimeFragment,
    /// Missing junctions at both ends
    InternalFragment,
    /// Junction loss due to intron retention
    IntronRetention,
    /// Single-exon query within a multi-exonic reference exon
    MonoExon,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum SubNIC {
    /// All junction pairs (donor-acceptor) individually known
    CombinationOfKnownJunctions,
    /// Individual sites known, but ≥1 novel pair combination
    CombinationOfKnownSpliceSites,
    /// Intron retention event
    IntronRetention,
    /// Mono-exon spanning a complete annotated intron
    MonoExonIntronRetention,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum SubNNC {
    /// ≥1 novel splice site (multi-exon)
    AtLeastOneNovelSpliceSite,
    /// Mono-exon NNC
    MonoExon,
}

impl ClassCode {
    /// SQANTI3 `structural_category` column value
    pub fn main_category(&self) -> &'static str {
        match self {
            Self::FSM(_) => "full-splice_match",
            Self::ISM(_) => "incomplete-splice_match",
            Self::NIC(_) => "novel_in_catalog",
            Self::NNC(_) => "novel_not_in_catalog",
            Self::Fusion => "fusion",
            Self::Antisense => "antisense",
            Self::GenicIntron => "genic_intron",
            Self::GenicGenomic => "genic",
            Self::Intergenic => "intergenic",
            Self::MoreJunctions => "moreJunctions",
        }
    }

    /// SQANTI3 `subcategory` column value (semicolon-delimited)
    pub fn sub_category(&self) -> &'static str {
        match self {
            // FSM
            Self::FSM(SubFSM::ReferenceMatch) => "multi-exon;reference_match",
            Self::FSM(SubFSM::Alternative3End) => "multi-exon;alternative_3end",
            Self::FSM(SubFSM::Alternative5End) => "multi-exon;alternative_5end",
            Self::FSM(SubFSM::Alternative5And3End) => "multi-exon;alternative_5end_3end",
            Self::FSM(SubFSM::MonoExon) => "mono-exon",
            // ISM
            Self::ISM(SubISM::ThreePrimeFragment) => "multi-exon;3prime_fragment",
            Self::ISM(SubISM::FivePrimeFragment) => "multi-exon;5prime_fragment",
            Self::ISM(SubISM::InternalFragment) => "multi-exon;internal_fragment",
            Self::ISM(SubISM::IntronRetention) => "multi-exon;intron_retention",
            Self::ISM(SubISM::MonoExon) => "mono-exon",
            // NIC
            Self::NIC(SubNIC::CombinationOfKnownJunctions) => {
                "multi-exon;combination_of_known_junctions"
            }
            Self::NIC(SubNIC::CombinationOfKnownSpliceSites) => {
                "multi-exon;combination_of_known_splicesites"
            }
            Self::NIC(SubNIC::IntronRetention) => "multi-exon;intron_retention",
            Self::NIC(SubNIC::MonoExonIntronRetention) => "mono-exon;intron_retention",
            // NNC
            Self::NNC(SubNNC::AtLeastOneNovelSpliceSite) => {
                "multi-exon;at_least_one_novel_splicesite"
            }
            Self::NNC(SubNNC::MonoExon) => "mono-exon",
            // No subcategory for these
            Self::Fusion => "multi-exon",
            Self::Antisense => "multi-exon",
            Self::GenicIntron => "genic_intron",
            Self::GenicGenomic => "genic_genomic",
            Self::Intergenic => "intergenic",
            Self::MoreJunctions => "multi-exon",
        }
    }

    /// Whether this classification involves junction-level matching
    /// (useful for deciding whether associated_transcript is meaningful)
    pub fn is_junction_based(&self) -> bool {
        matches!(
            self,
            Self::FSM(_) | Self::ISM(_) | Self::NIC(_) | Self::NNC(_) | Self::Fusion
        )
    }

    /// Whether the query is mono-exonic
    pub fn is_mono_exon(&self) -> bool {
        matches!(
            self,
            Self::FSM(SubFSM::MonoExon)
                | Self::ISM(SubISM::MonoExon)
                | Self::NIC(SubNIC::MonoExonIntronRetention)
                | Self::NNC(SubNNC::MonoExon)
        )
    }
}

impl fmt::Display for ClassCode {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.main_category())
    }
}
