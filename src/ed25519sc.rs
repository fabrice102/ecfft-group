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
/// `a = 4392976802491101277119858233748628886670447097743165573553979461609125550624`,
/// `b = 2641390116504058046593982364855935054624196913202961355189967522292222358134`.
pub struct Ed25519EcFftParameters;

impl EcFftParameters<F> for Ed25519EcFftParameters {
    /// The curve `E` has order multiple of 2^10
    const LOG_N: usize = 15;

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
        15 // must be equal to LOG_N
    }

    ecfft_tests! {
        ed25519pt,
        ark_ed25519::EdwardsProjective,
        super::Ed25519EcFftParameters,
        5, // select 0 for full tests (slower)
        15 // must be equal to LOG_N
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
