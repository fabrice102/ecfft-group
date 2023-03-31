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

    /// Evaluates `self` at the given `point` in `Self::Point`.
    pub fn evaluate(&self, point: &G::ScalarField) -> G {
        if self.is_zero() {
            return G::zero();
        } else if point.is_zero() {
            return self.coeffs[0].clone(); // TODO: while dense.rs does not need clone here? strange
        }
        self.internal_evaluate(point)
    }

    #[inline]
    // Horner's method for polynomial evaluation
    fn horner_evaluate(poly_coeffs: &[G], point: &G::ScalarField) -> G {
        // TODO: not optimal in case of elliptic curve
        // need to use multiscalar mult
        poly_coeffs
            .iter()
            .rfold(G::zero(), move |result, coeff| result * point + coeff)
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
