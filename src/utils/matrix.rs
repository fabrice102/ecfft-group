use ark_ff::PrimeField;

use crate::my_group::MyGroup;

#[derive(Clone, Copy)]
/// 2x2 matrix.
pub struct Matrix<T>(pub [[T; 2]; 2]);

impl<F: PrimeField> Matrix<F> {
    /// Inverse of the matrix. Panics if the matrix is not invertible.
    pub fn inverse(&self) -> Self {
        let [[a, b], [c, d]] = self.0;
        let det = a * d - b * c;
        Self([[d / det, -b / det], [-c / det, a / det]])
    }

    #[allow(clippy::many_single_char_names)]
    /// Multiply a vector of 2 field elements by the matrix.
    pub fn multiply<G: MyGroup<ScalarField = F>>(&self, v: [G; 2]) -> [G; 2] {
        let [[a, b], [c, d]] = self.0;
        let [x, y] = v;
        [x * a + y * b, x * c + y * d]
    }

    #[allow(clippy::many_single_char_names)]
    /// Multiply a vector of 2 field elements by the matrix.
    pub fn multiply_in_place<G: MyGroup<ScalarField = F>>(&self, x: &mut G, y: &mut G) {
        // TODO: use multiscalar
        let [[a, b], [c, d]] = self.0;
        let (a, b) = (*x * a + *y * b, *x * c + *y * d);
        *x = a;
        *y = b;
    }
}

#[cfg(test)]
mod tests {
    use ark_std::{rand::Rng, test_rng};

    use crate::bn254::F;

    use super::Matrix;

    #[test]
    fn test_inverse() {
        let mut rng = test_rng();
        for _ in 0..100 {
            let a: F = rng.gen();
            let b: F = rng.gen();
            let c: F = rng.gen();
            let d: F = rng.gen();
            let mat = Matrix([[a, b], [c, d]]);
            let mat_inv = mat.inverse();
            let x: F = rng.gen();
            let y: F = rng.gen();
            let v = [x, y];

            assert_eq!(v, mat_inv.multiply(mat.multiply(v)));
            assert_eq!(v, mat.multiply(mat_inv.multiply(v)));
        }
    }
}
