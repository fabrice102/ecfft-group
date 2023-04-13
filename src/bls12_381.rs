use std::convert::TryInto;

use crate::{ecfft::read_ecfft, ecfft::EcFftParameters, utils::isogeny::Isogeny};
use ark_ff::BigInteger384;
use crate::my_group::MyGroup;

type F = ark_bls12_381::Fq;
/// Number of 64-bit limbs needed to represent field elements.
const NUM_LIMBS: usize = 6;

/// ECFFT parameters for the BLS12-381 base field `F`.
/// Computed with the curve `E = EllipticCurve(F, [a, b])` with
/// `a, b = 0x287cc81c41f14f729fcbc12f57b2dd49bdcfc64938f9ad946c9fe5288aa3e9653670d336b09c058baad66ae717c1df7, 0x33f44f9b6fd7ba0080f0ad4843e076da70b11e6846d41e19792a15a4920e2294f9c971db67257eefea71c70514c6e54`
pub struct Bls12381Parameters;

impl EcFftParameters<F> for Bls12381Parameters {
    /// The curve `E` has order `4002409555221667393417789825735904156556882819939007885330472032889288775404654397856791416969022033299503997812736`
    /// with factorization `2^15 * 122143846289723736371392511771725590715236902463958980875563721706826439679097119075219464629181580606063964777`
    const LOG_N: usize = 15;

    const N: usize = 1 << Self::LOG_N;

    fn coset() -> Vec<F> {
        read_ecfft::read_coset::<F, NUM_LIMBS, _>("bls12-381_coset", |v| {
            BigInteger384::new(v.try_into().unwrap()).into()
        })
    }

    fn isogenies() -> Vec<Isogeny<F>> {
        read_ecfft::read_isogenies::<F, NUM_LIMBS, _>("bls12-381_isogenies", |v| {
            BigInteger384::new(v.try_into().unwrap()).into()
        })
    }
}

impl MyGroup for ark_bls12_381::Fq {
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
        bls12_381,
        super::F,
        super::Bls12381Parameters,
        0, // 0 = full tests, higher value allows tests to run faster
        15
    }
}
