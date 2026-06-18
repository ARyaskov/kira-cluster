use std::time::Instant;

#[derive(Debug, Clone)]
pub struct StageTimer {
    pub name: &'static str,
    pub start: Instant,
    pub elapsed_ns: u128,
}

impl StageTimer {
    pub fn new(name: &'static str) -> Self {
        Self {
            name,
            start: Instant::now(),
            elapsed_ns: 0,
        }
    }

    pub fn stop(&mut self) {
        self.elapsed_ns += self.start.elapsed().as_nanos();
    }
}
