#[derive(Debug, Default)]
pub(crate) enum GroupIterState {
    #[default]
    Precursor,
    Product(usize),
    Done,
}