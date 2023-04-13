//! ECFFT for the scalar field / curve points of ED25519
//! Note this is different from BN254 which is for the base field of BN254

use std::convert::TryInto;

use crate::{ecfft::read_ecfft, ecfft::EcFftParameters, utils::isogeny::Isogeny};
use ark_ff::BigInteger256;
use ark_ec::Group;
use crate::my_group::MyGroup;

pub type F = ark_ed25519::Fr;
pub type G = ark_ed25519::EdwardsProjective;
/// Number of 64-bit limbs needed to represent field elements.
const NUM_LIMBS: usize = 4;

/// ECFFT parameters for the ED22519 scalar field `F`.
/// Computed with the curve `E = EllipticCurve(F, [a, b])` with
/// `a = 7237005577332262213973186563042994240857116359379907606001950938285454250989`,
/// `b = 5612291247948481584627780310922020304781354847659642188369727566000581075360`.
pub struct Ed25519EcFftParameters;

impl EcFftParameters<F> for Ed25519EcFftParameters {
    /// The curve `E` has order multiple of 2^10
    const LOG_N: usize = 10;

    const N: usize = 1 << Self::LOG_N;

    fn coset() -> Vec<F> {
        read_ecfft::read_coset::<F, NUM_LIMBS, _>("ed25519sc_coset", |v| {
            BigInteger256::new(v.try_into().unwrap()).into()
        })
    }

    fn isogenies() -> Vec<Isogeny<F>> {
        read_ecfft::read_isogenies::<F, NUM_LIMBS, _>("ed25519sc_isogenies", |v| {
            BigInteger256::new(v.try_into().unwrap()).into()
        })
    }
}

impl MyGroup for ark_ed25519::Fr {
    type ScalarField = Self;

    fn generator() -> Self {
        Self::from(1)
    }
}

impl MyGroup for ark_ed25519::EdwardsProjective {
    type ScalarField = ark_ed25519::Fr;

    fn generator() -> Self {
        <Self as Group>::generator()
    }
}


#[cfg(test)]
mod tests {
    use ark_std::{rand::Rng, test_rng};
    use super::*;
    use crate::ecfft::ecfft_tests::ecfft_tests;

    ecfft_tests! {
        ed25519sc,
        super::F,
        super::Ed25519EcFftParameters,
        0,
        10
    }

    ecfft_tests! {
        ed25519pt,
        ark_ed25519::EdwardsProjective,
        super::Ed25519EcFftParameters,
        5, // select 0 for full tests (slower)
        10
    }

    #[test]
    pub fn test_scalar_mult_2() {
        // Multiply a random point by 2
        let mut rng = test_rng();
        let pt: G = rng.gen();
        let sc: F = F::from(2);
        let _ = pt * sc;
    }

}
