#![feature(generic_const_exprs)]
#![allow(incomplete_features, unused)]

use std::cell::UnsafeCell;
use std::iter::Sum;
use std::ops::{Add, Mul, Sub};
use std::sync::Mutex;
use num_traits::{One, Zero};
use once_cell::sync::Lazy;
use crate::prime_ring::r#static::{MOD_Q, PHI, TWO_PHI_MINUS_ONE, TEST_PHI, TEST_TWO_PHI_MINUS_ONE, BASIS};
use crate::cpp::bindings::{eltwise_add_mod, eltwise_sub_mod, eltwise_mult_mod, eltwise_reduce_mod, multiply_poly};
use crate::helpers::println_with_timestamp;
use rand::Rng;
use crate::arithmetic::{add_mod, call_sage_inverse_polynomial, first_n_columns, multiply_mod, random, reduce_mod, reduce_mod_vec, sub_mod};

pub type DPrimeRingElement = PrimeRingElement<PHI, TWO_PHI_MINUS_ONE, MOD_Q>;

#[derive(Clone, Copy, Debug)]
pub struct PrimeRingElement<
    const phi: usize,
    const two_phi_minus_one: usize,
    const mod_q: u64
> {
    pub coeffs: [u64; phi]
}


impl<
    const phi: usize,
    const two_phi_minus_one: usize,
    const mod_q: u64
> Add for PrimeRingElement<phi, two_phi_minus_one, mod_q> {
    type Output = Self;
    #[cfg(any(target_arch = "aarch64", not(feature = "use-hardware")))]
    // Implements the addition of two PrimeRingElements
    fn add(self, other: Self) -> Self {
        // let nof_adds_curr = *nof_adds.lock().unwrap();
        // *nof_adds.lock().unwrap() = nof_adds_curr + 1;
        let mut result = [0; phi];
        for i in 0..phi {
            result[i] = (self.coeffs[i] + other.coeffs[i]) % mod_q;
            if result[i] < 0 { result[i] += mod_q; }
        }
        PrimeRingElement { coeffs: result }
    }
    #[cfg(all(target_arch = "x86_64", feature = "use-hardware"))]
    fn add(self, other: Self) -> Self {
        let mut result = [0; phi];
        unsafe {
            eltwise_add_mod(result.as_mut_ptr(), self.coeffs.as_ptr(), other.coeffs.as_ptr(), phi as u64, mod_q);
        }
        PrimeRingElement { coeffs: result }
    }
}

impl<
    const phi: usize,
    const two_phi_minus_one: usize,
    const mod_q: u64
> Sum for PrimeRingElement<phi, two_phi_minus_one, mod_q> {
    fn sum<I: Iterator<Item=Self>>(iter: I) -> Self {
        iter.fold(PrimeRingElement::zero(), |acc, x| acc + x)
    }
}



impl<
    const phi: usize,
    const two_phi_minus_one: usize,
    const mod_q: u64
> Sub for PrimeRingElement<phi, two_phi_minus_one, mod_q> {
    type Output = Self;
    #[cfg(any(target_arch = "aarch64", not(feature = "use-hardware")))]
    // Implements the addition of two PrimeRingElements
    fn sub(self, other: Self) -> Self {
        let mut result = [0; phi];
        for i in 0..phi {
            result[i] = (mod_q + self.coeffs[i] - other.coeffs[i]) % mod_q;
        }
        PrimeRingElement { coeffs: result }
    }
    #[cfg(all(target_arch = "x86_64", feature = "use-hardware"))]
    fn sub(self, other: Self) -> Self {
        let mut result = [0; phi];
        unsafe {
            eltwise_sub_mod(result.as_mut_ptr(), self.coeffs.as_ptr(), other.coeffs.as_ptr(), phi as u64, mod_q);
        }
        PrimeRingElement { coeffs: result }
    }
}
#[cfg(any(target_arch = "aarch64", not(feature = "use-hardware")))]
pub fn poly_mul_mod<const phi:usize, const two_phi_minus_one:usize>
(a: &[u64; phi], b: &[u64; phi], mod_q: u64) -> [u64; two_phi_minus_one] {
    let mut result = [0u64; two_phi_minus_one];
    for i in 0..phi {
        for j in 0..phi {
            result[i + j] = (((result[i + j] as u128) + (a[i] as u128) * (b[j] as u128)) % (mod_q as u128)) as u64;
        }
    }
    result
}
#[cfg(all(target_arch = "x86_64", feature = "use-hardware"))]
pub fn poly_mul_mod<const phi:usize, const two_phi_minus_one:usize>
(a: &[u64; phi], b: &[u64; phi], mod_q: u64) -> [u64; two_phi_minus_one] {
    let n = a.len().next_power_of_two() * 2;
    let mut padded_a = vec![0u64; n];
    padded_a[..a.len()].copy_from_slice(a);
    let mut padded_b = vec![0u64; n];
    padded_b[..b.len()].copy_from_slice(b);
    let mut result = vec![0u64; n];

    unsafe {
        multiply_poly(result.as_mut_ptr(), padded_a.as_ptr(), padded_b.as_ptr(), n as u64, mod_q);
    }

    let mut shrunk_result = [0u64; two_phi_minus_one];
    shrunk_result.copy_from_slice(&result[0..two_phi_minus_one]);

    shrunk_result
}


impl<
    const phi: usize,
    const two_phi_minus_one: usize,
    const mod_q: u64
> Mul for PrimeRingElement<phi, two_phi_minus_one, mod_q> {
    type Output = Self;

    fn mul(self, other: Self) -> Self {
        // let nof_mul_curr = *nof_mul.lock().unwrap();
        // *nof_mul.lock().unwrap() = nof_mul_curr + 1;

        let mut result = poly_mul_mod::<phi, two_phi_minus_one>(&self.coeffs, &other.coeffs, mod_q);

        // Apply the cyclotomic polynomial reduction
        for i in phi + 1..result.len() {
            result[i - phi - 1] = result[i - phi - 1] + result[i]
        }

        let mut reduced_result = [0; phi];
        for i in 0..phi {
            reduced_result[i] = ((mod_q) + result[i] - result[phi]);
        }


        #[cfg(all(target_arch = "x86_64", feature = "use-hardware"))]
        unsafe {
            eltwise_reduce_mod(reduced_result.as_mut_ptr(), reduced_result.as_mut_ptr(), phi as u64, mod_q);
        }
        #[cfg(any(target_arch = "aarch64", not(feature = "use-hardware")))]
        {
            for i in 0..phi {
                reduced_result[i] = reduced_result[i] % mod_q;
            }
        }


        PrimeRingElement { coeffs: reduced_result.try_into().unwrap() }
    }
}

impl<
    const phi: usize,
    const two_phi_minus_one: usize,
    const mod_q: u64
> Zero for PrimeRingElement<phi, two_phi_minus_one, mod_q> {
    fn zero() -> Self {
        PrimeRingElement {
            coeffs: [0; phi],
        }
    }

    fn is_zero(&self) -> bool {
        self.coeffs.iter().all(|&coeff| coeff == 0)
    }
}


impl<
    const phi: usize,
    const two_phi_minus_one: usize,
    const mod_q: u64
> One for PrimeRingElement<phi, two_phi_minus_one, mod_q> {
    fn one() -> Self {
        let mut coeffs = [0; phi];
        coeffs[0] = 1;
        PrimeRingElement { coeffs }
    }
}

impl<
    const PHI: usize,
    const TWO_PHI_MINUS_ONE: usize,
    const mod_q: u64
> PartialEq for PrimeRingElement<PHI, TWO_PHI_MINUS_ONE, mod_q> {
    fn eq(&self, other: &Self) -> bool {
        self.coeffs == other.coeffs
    }
}


impl<
    const phi: usize,
    const two_phi_minus_one: usize,
    const mod_q: u64
> PrimeRingElement<phi, two_phi_minus_one, mod_q> {
    pub fn conjugate(&self) -> PrimeRingElement<phi, two_phi_minus_one, mod_q> {
        let mut coeffs = [0; phi];

        coeffs[0] += self.coeffs[0];


        #[cfg(all(target_arch = "x86_64", feature = "use-hardware"))]
        unsafe {
            let coeff_1 = [self.coeffs[1]; phi];
            eltwise_sub_mod(
                coeffs.as_mut_ptr(),
                coeffs.as_mut_ptr(),
                coeff_1.as_ptr(),
                phi as u64,
                mod_q
            );
            let mut rev_coeffs = self.coeffs.clone();
            rev_coeffs.reverse();
            eltwise_add_mod(
                coeffs.as_mut_ptr().add(2),             // result[2..phi]
                coeffs.as_ptr().add(2),                 // operand1[2..phi]
                rev_coeffs.as_ptr(),                    // operand2[0..phi-2]
                (phi - 2) as u64,
                mod_q,
            );
        }
        #[cfg(any(target_arch = "aarch64", not(feature = "use-hardware")))]
        {
            for i in 0..coeffs.len() {
                coeffs[i] += mod_q - self.coeffs[1];
            }


            for i in 2..phi {
                coeffs[phi - i + 1] += self.coeffs[i];
            }

            reduce_mod_vec(&mut coeffs, mod_q);
        }



        PrimeRingElement {
            coeffs
        }
    }

    pub fn trace(&self) -> u64 {
        let mut result = 0u64;
        result = multiply_mod(self.coeffs[0], phi as u64 - 1, mod_q);
        for i in 0..PHI {
            result = sub_mod(result, self.coeffs[i], mod_q);
        }
        result
    }

    pub fn inverse(&self) -> PrimeRingElement<phi, two_phi_minus_one, mod_q> {
        call_sage_inverse_polynomial(self)
    }

    pub fn negative(&self) -> PrimeRingElement<phi, two_phi_minus_one, mod_q> {
        PrimeRingElement::<phi, two_phi_minus_one, mod_q>::zero() - *self
    }
}

pub struct PrimeRing <
    const phi: usize = PHI,
    const two_phi_minus_one: usize = TWO_PHI_MINUS_ONE,
    const mod_q: u64 = MOD_Q,
>;

impl PrimeRing {
    // Constructor for PrimeRingElement
    pub fn new(coeffs: [u64; PHI]) -> DPrimeRingElement {
        DPrimeRingElement { coeffs }
    }

    pub fn new_rns(coeffs: [u64; PHI]) -> DPrimeRingElement {
        DPrimeRingElement { coeffs }
    }
    pub fn new_from_larger_vec(coeffs: &Vec<u64>) -> DPrimeRingElement {
        let mut reduced_result = [0; PHI];
            // Apply the cyclotomic polynomial reduction
        let mut result = coeffs[0..PHI + 1].to_vec();

        for i in PHI + 1..coeffs.len() {
            result[i % (PHI + 1)] += coeffs[i]
        }

        for i in 0..PHI {
            reduced_result[i] = (MOD_Q + result[i] - result[PHI]);
        }

        #[cfg(all(target_arch = "x86_64", feature = "use-hardware"))]
        unsafe {
            eltwise_reduce_mod(reduced_result.as_mut_ptr(), reduced_result.as_mut_ptr(), PHI as u64, MOD_Q);
            // eltwise_reduce_mod(reduced_result.as_mut_ptr(), reduced_result.as_mut_ptr(), PHI as u64, MOD_Q);
        }
        #[cfg(any(target_arch = "aarch64", not(feature = "use-hardware")))]
        {
            for i in 0..PHI {
                reduced_result[i] = reduced_result[i] % MOD_Q;
            }
        }


        PrimeRingElement { coeffs: reduced_result }
    }
    pub fn new_small_tests_only(coeffs: [u64; TEST_PHI]) -> PrimeRingElement<TEST_PHI, TEST_TWO_PHI_MINUS_ONE, MOD_Q> {
        PrimeRingElement::<TEST_PHI, TEST_TWO_PHI_MINUS_ONE, MOD_Q> { coeffs }
    }

    pub fn constant(v: u64) -> DPrimeRingElement {
        let mut coeffs = [0; PHI];
        coeffs[0] = v;
        PrimeRingElement { coeffs }
    }
    pub fn all(v: u64) -> DPrimeRingElement<> {
        PrimeRing::new( [v; PHI])
    }
    pub fn random_all() -> DPrimeRingElement<> {
        // TODO
        PrimeRing::all( random(1, MOD_Q)[0])
    }

    pub fn sample_subtractive() -> DPrimeRingElement {
        let mut rng = rand::thread_rng();

        // Randomly select an index i from 0 to PHI
        let i: usize = rng.gen_range(0..=PHI);

        // Initialize coefficients array
        let mut coeffs= [0u64; PHI];

        // Set coefficients to 1 from 0 to index i
        for index in 0..i {
            coeffs[index] = 1;
        }

        PrimeRingElement { coeffs }
    }

    pub fn random() -> DPrimeRingElement {
        PrimeRing::new(random(PHI, MOD_Q).try_into().unwrap())
    }
    pub fn random_real() -> DPrimeRingElement {
        let t = PrimeRing::new(random(PHI, MOD_Q).try_into().unwrap());
        t + t.conjugate()
    }

    pub fn random_constant() -> DPrimeRingElement {
        // TODO
        PrimeRing::constant(
            random(1, MOD_Q)[0]
        )
    }

    pub fn random_short() -> DPrimeRingElement {
        PrimeRing::new(random(PHI, 10).try_into().unwrap())
    }

    pub fn random_bin() -> DPrimeRingElement {
        PrimeRing::new(random(PHI, 2).try_into().unwrap())
    }

    pub fn random_from(l: u64) -> DPrimeRingElement {
        PrimeRing::new(random(PHI, l).try_into().unwrap())
    }

    pub fn random_constant_from(l: u64) -> DPrimeRingElement {
        PrimeRing::constant(
            random(1, l)[0]
        )
    }

    // pub fn embed_in_dual(coeffs: [[u64; PHI]; NOF_PRIMES]) -> DPrimeRingElement {
    //     (0..PHI).into_iter().map(|i| {
    //         let mut c_coeffs = 0u64;
    //         for j in 0..NOF_PRIMES {
    //             c_coeffs[j] = coeffs[i];
    //         }
    //         PrimeRing::constant_rns(c_coeffs) * DUAL_BASIS[i].clone()
    //     }).sum()
    // }

    // pub fn embed_constant_in_dual(coeffs: [u64; NOF_PRIMES]) -> DPrimeRingElement {
    //     PrimeRing::constant_rns(coeffs) * DUAL_BASIS[0].clone()
    // }
    //
    //
    // pub fn random_all_in_dual() -> DPrimeRingElement {
    //     let coeffs = random(NOF_PRIMES, MOD_Q);
    //     (0..PHI).into_iter().map(|i| {
    //         PrimeRing::constant_rns(coeffs.clone().try_into().unwrap()) * DUAL_BASIS[i].clone()
    //     }).sum()
    // }
}





#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_addition() {
        let a_coeffs: [u64; TEST_PHI] = [1, 2, 0, 0, 0, 0];
        let b_coeffs: [u64; TEST_PHI] = [3, 4, 0, 0, 0, 0];

        let a = PrimeRing::new_small_tests_only(a_coeffs);
        let b = PrimeRing::new_small_tests_only(b_coeffs);

        let c = a + b;

        assert_eq!(c.coeffs[0], 4);
        assert_eq!(c.coeffs[1], 6);
        assert_eq!(c.coeffs[2], 0);
        assert_eq!(c.coeffs[3], 0);
        assert_eq!(c.coeffs[4], 0);
        assert_eq!(c.coeffs[5], 0);
    }

    #[test]
    fn test_multiplication() {
        let a_coeffs: [u64; TEST_PHI] = [1, 2, 3, 4, 5, 6];
        let b_coeffs: [u64; TEST_PHI] = [1, 2, 3, 4, 5, 6];

        let a = PrimeRing::new_small_tests_only(a_coeffs);
        let b = PrimeRing::new_small_tests_only(b_coeffs);

        let c = a * b;
        println!("{:?}", c.coeffs);

        assert_eq!(c.coeffs, [7, 7, 0, MOD_Q - 14, MOD_Q - 35, MOD_Q - 14]);
    }



    // #[test]
    // fn test_multiplication_big() {
    //     let a_coeffs: [u64; TEST_PHI] = [4611686019232684273, 4611684019232694273, 4610686019232694273, 4600686019232694273, 4611686019232694240, 2611686019232694273];
    //     let b_coeffs: [u64; TEST_PHI] = [4611686019232684273, 4611684019232694273, 4610686019232694273, 4600686019232694273, 4611686019232694240, 2611686019232694273];
    //     const modulus:u64 = 4611686019232694273;
    //     let a = PrimeRingElement::<TEST_PHI, TEST_TWO_PHI_MINUS_ONE, modulus>{ coeffs: a_coeffs };
    //     let b = PrimeRingElement::<TEST_PHI, TEST_TWO_PHI_MINUS_ONE, modulus>{ coeffs: a_coeffs };
    //
    //     let c = a * b;
    //     println!("{:?}", c.coeffs);
    //
    //     assert_eq!(c.coeffs, [40073398178371627, 1809422541500652510, 1363234956273370231, 1760542567796230419, 184234902343308818, 1834626940863330945]);
    //
    // }
}

// #[test]
// fn test_trace() {
//     let a = PrimeRing::random_from(FHE_Q);
//     let (short_trace, rem) = a.trace_with_remainder();
//     let mut composed = 0u64;
//     composed[j] = add_mod(short_trace[j], multiply_mod(rem[j], FHE_Q, MOD_Q), MOD_Q)
//     println!("{:?}",rem);
//     assert_eq!(composed, a.trace());
// }

//
// #[test]
// fn test_trace_dual() {
//     for i in 0..PHI {
//         for j in 0..PHI {
//             if i == j { continue }
//             assert_eq!(0, (BASIS[j] * DUAL_BASIS[i]).trace()[0]);
//         }
//     }
// }

