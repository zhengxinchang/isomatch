pub trait ArgValidate {
    fn validate(&self);
}

pub trait LogMemSize {
    fn get_mem_size(&self) -> usize;
}
