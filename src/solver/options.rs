#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SolverOptions {
    pub tolerance: f64,
    pub max_iterations: u32,
    pub lambda_scale: f64,
}

impl Default for SolverOptions {
    fn default() -> Self {
        Self {
            tolerance: 1.0e-6,
            max_iterations: 100,
            lambda_scale: 0.5,
        }
    }
}
