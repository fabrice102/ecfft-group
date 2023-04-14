use std::convert::TryInto;

use crate::{ecfft::read_ecfft, ecfft::EcFftParameters, utils::isogeny::Isogeny};
use ark_ff::BigInteger256;
use crate::my_group::MyGroup;

type F = ark_curve25519::Fq;
/// Number of 64-bit limbs needed to represent field elements.
const NUM_LIMBS: usize = 4;

/// ECFFT parameters for the Curve25519 base field `F`.
/// Computed with the curve `E = EllipticCurve(F, [a, b])` with
/// `a, b = 0x1, 0xd63`
pub struct Curve25519Parameters;

impl EcFftParameters<F> for Curve25519Parameters {
    /// The curve `E` has order `57896044618658097711785492504343953926261577544886303154527763127846362480640`
    /// with factorization `2^16 * 883423532389192164791648750371459257908044090955906725380367479367772865`
    const LOG_N: usize = 16;

    const N: usize = 1 << Self::LOG_N;

    fn coset() -> Vec<F> {
        read_ecfft::read_coset::<F, NUM_LIMBS, _>("curve25519_coset", |v| {
            BigInteger256::new(v.try_into().unwrap()).into()
        })
    }

    fn isogenies() -> Vec<Isogeny<F>> {
        read_ecfft::read_isogenies::<F, NUM_LIMBS, _>("curve25519_isogenies", |v| {
            BigInteger256::new(v.try_into().unwrap()).into()
        })
    }
}

impl MyGroup for ark_curve25519::Fq {
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
        curve25519,
        super::F,
        super::Curve25519Parameters,
        0, // 0 = full tests, higher value allows tests to run faster
        16
    }
}
