///! This module provides functions to multiply Vandermonde matrices by vectors (on the left)

use ark_ff::PrimeField;
use crate::group_polynomial::smallest_range;
use crate::my_group::MyGroup;

/// Represent an "extended" Vandermonde matrix
/// M_{i,j} = points[i] ^ j, where 0 <= j < nb_cols
pub struct VandermondeMatrix<F> {
    pub points: Vec<F>,
    pub nb_cols: usize,
}

impl<F: PrimeField> VandermondeMatrix<F> {
    /// Multiply a Vandermonde matrix by a vector on the left
    /// Result is vector * self, where vector is seen as a row vector
    pub fn left_multiply<G: MyGroup<ScalarField=F>>(&self, vector: &[G]) -> Vec<G> {
        assert_eq!(self.points.len(), vector.len());
        let nb_rows = self.points.len();

        // at step j,
        //   col[i] = vector[i] * vandermonde_{i, j}
        //          = vector[i] * points[i]^j
        let mut col = vector.to_vec();
        // res is the multiplication result vector
        //   so res[j] is the sum of all col[i] at step j
        let mut res: Vec<G> = Vec::with_capacity(self.nb_cols);

        for j in 0..self.nb_cols {
            if j > 0 {
                // Update col to be the value as defined above
                // when j = 0, we start with the value "vector"
                // so no update needed
                for i in 0..nb_rows {
                    col[i] *= self.points[i];
                }
            }

            // Compute res
            res.push(col[0]);
            for i in 1..nb_rows {
                res[j] += col[i]
            }
        }

        return res;
    }
}

/// Same as VandermondeMatrix.left_multiply, where the points are
///   -floor((n-1)/2), ..., floor(n/2)
/// where n = len(vector)
/// For example, for n = 5, the points are:
///   points = -2, -1, 0, 1, 2
/// and for n = 4
///   points = -1, 0, 1, 2
/// Assumes negative of elements in G is very cheap
/// Does not make use of potentially optimized double
pub fn small_vandermonde_left_multiply<F, G>(vector: &[G], nb_cols: usize)
                                             -> Vec<G>
    where
        F: PrimeField,
        G: MyGroup<ScalarField=F>,
{
    let nb_rows = vector.len();
    assert!(nb_rows > 0);
    let zero_i = (nb_rows - 1) / 2; // the index of the point equal to 0
    // points[i] = i - zero_i
    // abs_points[i] = |points[i]| = |i - zero_i|, the absolute value
    let abs_points: Vec<F> = smallest_range(nb_rows as i32)
        .iter().map(|p| F::from(p.abs() as u128)).collect();

    // at step j,
    //   abs_col[i] = vector[i] * abs_points[i]^j
    //   except for abs_col[zero_i] which will always be set to vector[i]
    //   since we then don't add it in res[j] below
    let mut abs_col = vector.to_vec();
    // res is the multiplication result vector
    //   so res[j] is the sum of all +-abs_col[i] at step j
    //      (except the sum skips i = zero_i when j > 0),
    //   where +- depends on i and j
    let mut res: Vec<G> = Vec::with_capacity(nb_cols);

    for j in 0..nb_cols {
        if j > 0 {
            // Update col to be the value as defined above
            // when j = 0, we start with the value "vector"
            // so no update needed
            for i in 0..nb_rows {
                if i != zero_i {
                    abs_col[i] *= abs_points[i];
                }
            }
        }

        // Compute res
        if nb_rows == 1 && j > 0 {
            // special case when there is a single point = 0
            // for any column after j=0, the result is zero
            res.push(G::zero());
            continue
        }
        res.push(abs_col[nb_rows-1]);
        for i in 0..nb_rows-1 {
            if i < zero_i && j % 2 == 1 {
                res[j] -= abs_col[i];
            } else if i == zero_i {
                if j == 0 {
                    res[j] += abs_col[i]
                }
            } else {
                res[j] += abs_col[i];
            }
        }
    }

    return res;
}

#[cfg(test)]
mod tests {
    use crate::utils::f17::F17;
    use test_case::test_case;
    use crate::vandermonde::{small_vandermonde_left_multiply, VandermondeMatrix};

    fn to_vec_f17(v: &[u64]) -> Vec<F17> {
        v.iter().map(|x| F17::from(*x)).collect()
    }

    // First test case manually created
    // Following ones created from tests.ipynb
    // using sage
    #[test_case(vec ! [1, 2], vec ! [3, 4], vec ! [7, 11, 2, 1])]
    #[test_case(vec ! [3], vec ! [10], vec ! [10, 13])]
    #[test_case(vec ! [16], vec ! [16], vec ! [16])]
    #[test_case(vec ! [1, 0], vec ! [2, 9], vec ! [11, 2, 2])]
    #[test_case(vec ! [10, 16], vec ! [4, 10], vec ! [14, 13, 2, 12])]
    #[test_case(vec ! [5, 16, 2, 5, 14], vec ! [14, 1, 14, 3, 14], vec ! [12, 2])]
    #[test_case(vec ! [2, 8, 4, 4, 6], vec ! [12, 14, 3, 13, 10], vec ! [1, 5, 13])]
    pub fn test_vandermonde_left_multiply(
        points: Vec<u64>,
        vector: Vec<u64>,
        expected_res: Vec<u64>,
    ) {
        let points = to_vec_f17(&points[..]);
        let vector = to_vec_f17(&vector[..]);
        let expected_res = to_vec_f17(&expected_res[..]);

        let mat = VandermondeMatrix {
            points,
            nb_cols: expected_res.len(),
        };

        let actual_res = mat.left_multiply(&vector[..]);

        assert_eq!(expected_res, actual_res);
    }

    // Tests created from tests.ipynb
    // using sage
    #[test_case(vec ! [2], vec ! [2])]
    #[test_case(vec ! [14], vec ! [14, 0])]
    #[test_case(vec ! [0], vec ! [0, 0, 0, 0, 0])]
    #[test_case(vec ! [15, 11], vec ! [9])]
    #[test_case(vec ! [10, 16], vec ! [9, 16, 16, 16, 16])]
    #[test_case(vec ! [2, 9, 4], vec ! [15, 2, 6, 2, 6])]
    #[test_case(vec ! [10, 14, 1, 14], vec ! [5, 2, 16, 1, 14])]
    #[test_case(vec ! [3, 14, 12, 14, 3], vec ! [12, 0])]
    pub fn test_small_vandermonde_left_multiply(
        vector: Vec<u64>,
        expected_res: Vec<u64>,
    ) {
        let vector = to_vec_f17(&vector[..]);
        let expected_res = to_vec_f17(&expected_res[..]);

        let actual_res = small_vandermonde_left_multiply(
            &vector[..], expected_res.len());

        assert_eq!(expected_res, actual_res);
    }
}