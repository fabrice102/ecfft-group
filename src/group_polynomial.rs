// Inspired from dense.rs from ark-poly ...

use super::my_group::MyGroup;
use ark_ff::Zero;
use std::ops::{Deref, DerefMut};

/// Stores a polynomial in coefficient form
/// with coefficients in a group G
#[derive(Clone, PartialEq, Eq, Hash, Default)]
pub struct DenseGroupPolynomial<G: MyGroup> {
    /// The coefficient of `x^i` is stored at location `i` in `self.coeffs`.
    pub coeffs: Vec<G>,
}

impl<G: MyGroup> DenseGroupPolynomial<G> {
    /// Returns the total degree of the polynomial
    pub fn degree(&self) -> usize {
        if self.is_zero() {
            0
        } else {
            assert!(self.coeffs.last().map_or(false, |coeff| !coeff.is_zero()));
            self.coeffs.len() - 1
        }
    }

    /// Evaluates `self` at the given `point` in `Self::G::ScalarField`.
    /// More efficient than evaluate when |point| is small and multiplication
    /// of a group element by a small integer (and negating it) is faster
    /// than multiplying by a large integer
    pub fn evaluate_small(&self, point: i32) -> G {
        if self.is_zero() {
            return G::zero();
        } else if point == 0{
            return self.coeffs[0];
        } else if point < 0 {
            return Self::horner_evaluate_neg(&self.coeffs, &G::ScalarField::from((-point) as u128))
        } else {
            return Self::horner_evaluate(&self.coeffs, &G::ScalarField::from(point as u128))
        }
    }

    /// Evaluates `self` at the given `point` in `Self::G::ScalarField`.
    pub fn evaluate(&self, point: &G::ScalarField) -> G {
        if self.is_zero() {
            return G::zero();
        } else if point.is_zero() {
            return self.coeffs[0];
        }
        self.internal_evaluate(point)
    }

    #[inline]
    // Horner's method for polynomial evaluation
    fn horner_evaluate(poly_coeffs: &[G], point: &G::ScalarField) -> G {
        poly_coeffs
            .iter()
            .rfold(G::zero(), move |result, coeff| result * point + coeff)
    }

    #[inline]
    // Horner's method for polynomial evaluation to evaluate
    // the polynomial on -neg_point
    // faster when neg_point is small and multiplying a group element by neg_point
    // and negative the result
    // is faster than multiplying a group element by -neg_pooint
    fn horner_evaluate_neg(poly_coeffs: &[G], neg_point: &G::ScalarField) -> G {
        poly_coeffs
            .iter()
            .rfold(G::zero(), move |result, coeff| - result * neg_point + coeff)
    }

    fn internal_evaluate(&self, point: &G::ScalarField) -> G {
        Self::horner_evaluate(&self.coeffs, point)
    }
}

impl<G: MyGroup> DenseGroupPolynomial<G> {
    /// Returns the zero polynomial.
    pub fn zero() -> Self {
        Self { coeffs: Vec::new() }
    }

    /// Checks if the given polynomial is zero.
    pub fn is_zero(&self) -> bool {
        self.coeffs.is_empty() || self.coeffs.iter().all(|coeff| coeff.is_zero())
    }
}

impl<G: MyGroup> Deref for DenseGroupPolynomial<G> {
    type Target = [G];

    fn deref(&self) -> &[G] {
        &self.coeffs
    }
}

impl<G: MyGroup> DerefMut for DenseGroupPolynomial<G> {
    fn deref_mut(&mut self) -> &mut [G] {
        &mut self.coeffs
    }
}

/// smallest_range returns the vector
/// \[ -floor((n-1)/2),...,floor(n/2) \]
/// that contains the n distinct integers
pub fn smallest_range(n: i32) -> Vec<i32> {
    assert!(n >= 0);
    let mut v: Vec<i32> = Vec::with_capacity(n as usize);
    if n == 0 {
        return v
    }
    let min = -((n-1)/2);
    let max = n/2;
    for i in min..=max {
        v.push(i)
    }
    v
}

#[cfg(test)]
mod tests {
    use ark_std::{test_rng, UniformRand};
    use crate::group_polynomial::{DenseGroupPolynomial, smallest_range};

    /// evaluate a random polynomial with coefficients in ed25519
    /// using both normal evaluation and small one and compare they are equal
    #[test]
    fn evaluate_small_rng_ed25519() {
        type G = ark_ed25519::EdwardsProjective;
        type F = ark_ed25519::Fr;
        let mut rng = test_rng();
        let deg = 10;
        for _i in 0..10 {
            let coeffs: Vec<G> = (0..=deg).map(|_| G::rand(&mut rng)).collect();
            let poly = DenseGroupPolynomial { coeffs };
            let point: i32 = i32::rand(&mut rng);
            assert_eq!(poly.evaluate(&F::from(point)), poly.evaluate_small(point));
            assert_ne!(poly.evaluate(&F::from(point)), poly.evaluate_small(point+1));
        }
    }

    #[test]
    #[should_panic]
    fn smallest_range_neg() {
        smallest_range(-1);
    }

    #[test]
    fn smallest_range_0() {
        assert_eq!(Vec::<i32>::new(), smallest_range(0))
    }

    #[test]
    fn smallest_range_1() {
        assert_eq!(vec![0], smallest_range(1))
    }

    #[test]
    fn smallest_range_2() {
        assert_eq!(vec![0,1], smallest_range(2))
    }

    #[test]
    fn smallest_range_3() {
        assert_eq!(vec![-1,0,1], smallest_range(3))
    }

    #[test]
    fn smallest_range_8() {
        assert_eq!(vec![-3,-2,-1,0,1,2,3,4], smallest_range(8))
    }
}
