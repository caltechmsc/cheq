pub trait AtomView {
    fn atomic_number(&self) -> u8;

    fn position(&self) -> [f64; 3];
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Atom {
    pub atomic_number: u8,
    pub position: [f64; 3],
}

impl AtomView for Atom {
    #[inline(always)]
    fn atomic_number(&self) -> u8 {
        self.atomic_number
    }

    #[inline(always)]
    fn position(&self) -> [f64; 3] {
        self.position
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct CalculationResult {
    pub charges: Vec<f64>,
    pub equilibrated_potential: f64,
    pub iterations: u32,
}
