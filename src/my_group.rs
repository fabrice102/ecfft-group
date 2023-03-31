use ark_ff::{fields::PrimeField, UniformRand};
use ark_std::{
    fmt::{Debug, Display},
    hash::Hash,
    ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign},
};
use num_traits::Zero;
use zeroize::Zeroize;

pub trait MyGroup:
    Eq
    + 'static
    + Sized
    + Copy
    + Clone
    + Default
    + Send
    + Sync
    + Hash
    + Debug
    + Display
    + UniformRand
    + Zeroize
    + Zero
    + Neg<Output = Self>
    + Add<Self, Output = Self>
    + Sub<Self, Output = Self>
    + Mul<<Self as MyGroup>::ScalarField, Output = Self>
    + AddAssign<Self>
    + SubAssign<Self>
    + MulAssign<<Self as MyGroup>::ScalarField>
    + for<'a> Add<&'a Self, Output = Self>
    + for<'a> Sub<&'a Self, Output = Self>
    + for<'a> Mul<&'a <Self as MyGroup>::ScalarField, Output = Self>
    + for<'a> AddAssign<&'a Self>
    + for<'a> SubAssign<&'a Self>
    + for<'a> MulAssign<&'a <Self as MyGroup>::ScalarField>
    + core::iter::Sum<Self>
    + for<'a> core::iter::Sum<&'a Self>
{
    type ScalarField: PrimeField;

    /// Returns a fixed generator of this group.
    #[must_use]
    fn generator() -> Self;
}
