//! Tools to test easily a specific ECFFT implementation

#![cfg(test)]

use ark_std::test_rng;

use crate::ecfft::{
    EcFftCosetPrecomputation, EcFftParameters, EcFftPrecomputation, EcFftPrecomputationStep,
};
use crate::group_polynomial::DenseGroupPolynomial;
use crate::my_group::MyGroup;

/// Tests that precomputations don't panic.
pub fn test_precompute<G: MyGroup, P: EcFftParameters<G::ScalarField>>() {
    P::precompute_on_coset(&P::coset());
    P::precompute_on_coset(&P::coset().into_iter().step_by(2).collect::<Vec<_>>());
}

/// Tests the extend function with a polynomial of degree `2^(log(n)-i) - 1`.
/// 1 <= i <= log(n) - 1
pub fn test_extend_i<G: MyGroup, P: EcFftParameters<G::ScalarField>>(
    i: usize,
    precomputation: &EcFftCosetPrecomputation<G::ScalarField, P>,
) {
    let n = P::N >> i;
    let mut rng = test_rng();
    let coeffs: Vec<G> = (0..n).map(|_| G::rand(&mut rng)).collect();
    let poly = DenseGroupPolynomial { coeffs };
    let EcFftPrecomputationStep { s, s_prime, .. } = &precomputation.steps[i - 1];
    let evals_s = s.iter().map(|x| poly.evaluate(x)).collect::<Vec<_>>();
    let evals_s_prime = s_prime.iter().map(|x| poly.evaluate(x)).collect::<Vec<_>>();
    assert_eq!(evals_s_prime, precomputation.extend(&evals_s));
}

/// Tests the extend function for various degrees.
pub fn test_extend<G: MyGroup, P: EcFftParameters<G::ScalarField>>() {
    let precomputation = P::precompute_on_coset(&P::coset());
    for i in 1..P::LOG_N {
        test_extend_i::<G, P>(i, &precomputation);
    }
}

/// Tests the `evaluate_over_domain` function for various degrees.
pub fn test_eval_i<G: MyGroup, P: EcFftParameters<G::ScalarField>>(
    i: usize,
    precomputation: &EcFftPrecomputation<G::ScalarField, P>,
) {
    let mut rng = test_rng();
    let coeffs: Vec<G> = (0..P::N >> i).map(|_| G::rand(&mut rng)).collect();
    let poly = DenseGroupPolynomial { coeffs };
    let now = std::time::Instant::now();
    let evals = P::sub_coset(i)
        .iter()
        .map(|x| poly.evaluate(x))
        .collect::<Vec<_>>();
    dbg!(now.elapsed().as_secs_f32());
    assert_eq!(evals, precomputation.evaluate_over_domain(&poly));
    dbg!(now.elapsed().as_secs_f32());
}

/// generate tests for an instantiation of EcFftParameters
/// eval and extend are done for i=start_i to i=end_i-1
/// where i=0 means that a polynomial of max degree log_n
macro_rules! ecfft_tests {
    ($name:ident, $my_group_type:ty, $ecfft_parameters:ty, $start_i: expr, $end_i:expr) => {
        mod $name {
            use crate::ecfft::ecfft_tests;
            use crate::ecfft::{EcFftParameters, EcFftPrecomputation};
            use crate::my_group::MyGroup;

            use lazy_static::lazy_static;
            use seq_macro::seq;
            use test_case::test_case;

            type G = $my_group_type;
            type F = <G as MyGroup>::ScalarField;
            type P = $ecfft_parameters;

            // We need to lazily do the precomputation as it takes time to do them

            lazy_static! {
                static ref PRECOMPUTATION: EcFftPrecomputation<F, P> = {
                    P::precompute()
                };
            }

            #[test]
            fn test_precompute() {
                ecfft_tests::test_precompute::<G, P>();
            }

            seq!(N in $start_i..$end_i {
                #(#[test_case(N)])*
                fn test_eval_i(i: usize) {
                    ecfft_tests::test_eval_i::<G, P>(i, &*PRECOMPUTATION);
                }

                #(#[test_case(N)])*
                fn test_extend(i: usize) {
                    if i == 0 {
                        // skip as not supported for i=0
                        return
                    }
                    ecfft_tests::test_extend_i::<G, P>(i, &PRECOMPUTATION.coset_precomputations[0]);
                }
            });
        }
    }
}

pub(crate) use ecfft_tests;
