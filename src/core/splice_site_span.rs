#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Hash)]
pub struct SpliceSiteSpan {
    pub offset: u32,
    pub count: u16,
}

impl SpliceSiteSpan {
    pub fn is_empty(self) -> bool {
        self.count == 0
    }

    pub fn end_offset(self) -> u32 {
        self.offset + u32::from(self.count)
    }
}
