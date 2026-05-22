pub enum ClassCode {
    FSM(SubLevelFSM),
    ISM(SubLevelISM),
    NIC(SubLevelNIC),
    NNC(SubLevelNNC),
    Antisense(SubLevelAntisense),
    Genic(SubLevelGenic),
    Intergenic(SubLevelIntergenic),
    Fusion(SubLevelFusion),
}

pub enum SubLevelFSM {}

pub enum SubLevelISM {}

pub enum SubLevelNIC {}

pub enum SubLevelNNC {}

pub enum SubLevelAntisense {}

pub enum SubLevelGenic {}

pub enum SubLevelIntergenic {}

pub enum SubLevelFusion {}
