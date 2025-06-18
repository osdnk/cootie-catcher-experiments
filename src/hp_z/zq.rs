use num_traits::{Zero, One};
use rand::Rng;
use std::iter::Sum;
use crate::prime_ring::r#static::MOD_Q;

static Q: u64 = MOD_Q;

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct ZqElement {
    value: u64,
}

impl std::ops::Add for ZqElement {
    type Output = Self;
    fn add(self, other: Self) -> Self {
        ZqElement {
            value: (self.value + other.value) % Q,
        }
    }
}

impl std::ops::Sub for ZqElement {
    type Output = Self;
    fn sub(self, other: Self) -> Self {
        ZqElement {
            value: (self.value + Q - other.value) % Q,
        }
    }
}

impl std::ops::Mul for ZqElement {
    type Output = Self;
    fn mul(self, other: Self) -> Self {
        ZqElement {
            value: ((self.value as u128 * other.value as u128) % (Q as u128)) as u64
        }
    }
}

impl Zero for ZqElement {
    fn zero() -> Self {
        ZqElement { value: 0 }
    }
    fn is_zero(&self) -> bool {
        self.value == 0
    }
}

impl One for ZqElement {
    fn one() -> Self {
        ZqElement { value: 1 }
    }
    fn is_one(&self) -> bool {
        self.value == 1
    }
}

impl Sum for ZqElement {
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.fold(ZqElement { value: 0 }, |acc, x| acc + x)
    }
}

pub struct Zq;

impl Zq {
    pub fn zero() -> ZqElement {
        ZqElement { value: 0 }
    }

    pub fn one() -> ZqElement {
        ZqElement { value: 1 }
    }

    pub fn from_u64(value: u64) -> ZqElement {
        ZqElement { value: value % Q }
    }

    pub fn random() -> ZqElement {
        let mut rng = rand::thread_rng();
        ZqElement { value: rng.gen_range(0..Q) }
    }
}
