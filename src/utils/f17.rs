use ark_ff::{Fp64, MontBackend, MontConfig};
use crate::my_group::MyGroup;

///! Define the field F_17 used in some of our tests

#[derive(MontConfig)]
#[modulus = "17"]
#[generator = "3"]
pub struct F17Config;
pub type F17 = Fp64<MontBackend<F17Config, 1>>;

impl MyGroup for F17 {
    type ScalarField = Self;

    fn generator() -> Self {
        Self::from(1)
    }
}