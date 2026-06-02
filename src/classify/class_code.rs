use std::fmt;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum ClassCode {
    FSM(SubFSM),

    ISM(SubISM),

    NIC(SubNIC),

    NNC(SubNNC),

    Fusion,

    Antisense,

    GenicIntron,

    Genic,

    Intergenic,

    /// Produced when a query is associated with multiple genes, but at least
    /// one query junction is shared by more than one of those reference genes,
    /// so SQANTI3 does not call it `fusion`.
    MoreJunctions,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum SubFSM {
    ReferenceMatch,
    Alternative3End,
    Alternative5End,
    Alternative3And5End,
    MonoExon,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum SubISM {
    Complete, // left most and right most are match, but internal doesn't
    ThreePrimeFragment,
    FivePrimeFragment,
    InternalFragment, // ThreePrimeFragment + FivePrimeFragment
    IntronRetention,  // all splice site match, not all junction match
    MonoExon,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum SubNIC {
    CombinationOfKnownJunctions,
    CombinationOfKnownSpliceSites,
    IntronRetention,
    MonoExonByIntronRetention,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum SubNNC {
    AtLeastOneNovelSpliceSite,
}

impl ClassCode {
    pub fn main_category(&self) -> &'static str {
        match self {
            Self::FSM(_) => "full-splice_match",
            Self::ISM(_) => "incomplete-splice_match",
            Self::NIC(_) => "novel_in_catalog",
            Self::NNC(_) => "novel_not_in_catalog",
            Self::Fusion => "fusion",
            Self::Antisense => "antisense",
            Self::GenicIntron => "genic_intron",
            Self::Genic => "genic",
            Self::Intergenic => "intergenic",
            Self::MoreJunctions => "moreJunctions",
        }
    }

    pub fn sub_category(&self, exon_count: u16) -> &'static str {
        let exon_shape = if exon_count == 1 {
            "mono-exon"
        } else {
            "multi-exon"
        };

        match self {
            Self::FSM(SubFSM::ReferenceMatch) => "reference_match",
            Self::FSM(SubFSM::Alternative3End) => "alternative_3end",
            Self::FSM(SubFSM::Alternative5End) => "alternative_5end",
            Self::FSM(SubFSM::Alternative3And5End) => "alternative_3end5end",
            Self::FSM(SubFSM::MonoExon) => "mono-exon",

            Self::ISM(SubISM::Complete) => "complete",
            Self::ISM(SubISM::ThreePrimeFragment) => "3prime_fragment",
            Self::ISM(SubISM::FivePrimeFragment) => "5prime_fragment",
            Self::ISM(SubISM::InternalFragment) => "internal_fragment",
            Self::ISM(SubISM::IntronRetention) => "intron_retention",
            Self::ISM(SubISM::MonoExon) => "mono-exon",

            Self::NIC(SubNIC::CombinationOfKnownJunctions) => "combination_of_known_junctions",
            Self::NIC(SubNIC::CombinationOfKnownSpliceSites) => "combination_of_known_splicesites",
            Self::NIC(SubNIC::IntronRetention) => "intron_retention",
            Self::NIC(SubNIC::MonoExonByIntronRetention) => "mono-exon_by_intron_retention",

            Self::NNC(SubNNC::AtLeastOneNovelSpliceSite) => "at_least_one_novel_splicesite",

            Self::Fusion
            | Self::Antisense
            | Self::GenicIntron
            | Self::Genic
            | Self::Intergenic
            | Self::MoreJunctions => exon_shape,
        }
    }

    pub fn is_junction_based(&self) -> bool {
        matches!(
            self,
            Self::FSM(_)
                | Self::ISM(_)
                | Self::NIC(_)
                | Self::NNC(_)
                | Self::Fusion
                | Self::MoreJunctions
        )
    }

    pub fn is_mono_exon(&self, exon_count: u16) -> bool {
        exon_count == 1
            || matches!(
                self,
                Self::FSM(SubFSM::MonoExon)
                    | Self::ISM(SubISM::MonoExon)
                    | Self::NIC(SubNIC::MonoExonByIntronRetention)
            )
    }
}

impl fmt::Display for ClassCode {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.main_category())
    }
}
