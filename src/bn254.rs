use std::convert::TryInto;

use crate::{ecfft::read_ecfft, ecfft::EcFftParameters, utils::isogeny::Isogeny};
use ark_ff::BigInteger256;
use crate::my_group::MyGroup;

pub type F = ark_bn254::Fq;
/// Number of 64-bit limbs needed to represent field elements.
const NUM_LIMBS: usize = 4;

/// ECFFT parameters for the BN254 base field `F`.
/// Computed with the curve `E = EllipticCurve(F, [a, b])` with
/// `a, b = 1, 5612291247948481584627780310922020304781354847659642188369727566000581075360`.
pub struct Bn254EcFftParameters;

impl EcFftParameters<F> for Bn254EcFftParameters {
    /// The curve `E` has order `21888242871839275222246405745257275088712935808829559400805562964428910444544`
    /// with factorization `2^14 * 3^2 * 229 * 503 * 205460939795467 * 55374745393148401254803 * 113267149255983544517087125127`.
    const LOG_N: usize = 14;

    const N: usize = 1 << Self::LOG_N;

    fn coset() -> Vec<F> {
        read_ecfft::read_coset::<F, NUM_LIMBS, _>("bn254_coset", |v| {
            BigInteger256::new(v.try_into().unwrap()).into()
        })
    }

    fn isogenies() -> Vec<Isogeny<F>> {
        read_ecfft::read_isogenies::<F, NUM_LIMBS, _>("bn254_isogenies", |v| {
            BigInteger256::new(v.try_into().unwrap()).into()
        })
    }
}

impl MyGroup for ark_bn254::Fq {
    type ScalarField = Self;

    fn generator() -> Self {
        Self::from(1)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ecfft::ecfft_tests::ecfft_tests;

    ecfft_tests! {
        bn254,
        super::F,
        super::Bn254EcFftParameters,
        0, // 0 = full tests, higher value allows tests to run faster
        14
    }
}
