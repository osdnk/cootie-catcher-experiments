use std::sync::Mutex;
use std::time::Instant;
use fast_modulo::mulmod_u64;
use rayon::iter::{IntoParallelIterator, IntoParallelRefIterator};
use rayon::prelude::*;
use crate::arithmetic::{last_n_columns, reduce_mod_vec, sample_random_mat, transpose};
use crate::prime_ring::r#static::{MOD_Q, PHI};
use crate::prime_ring::ring::{DPrimeRingElement, PrimeRing, PrimeRingElement};
use crate::ntt::ntt::{inverse_ntt_pow_of_2, ntt_pow_of_2};


pub fn hadamard_64(a: &[u64], b: &[u64], mod_q: u64) -> Vec<u64> {
    let (short, long) = if a.len() < b.len() { (a, b) } else { (b, a) };

    let mut res:Vec<u64> = short
        .iter()
        .zip(long.iter())
        .map(|(aa, ba)| mulmod_u64(*aa, *ba, mod_q))
        .collect();

    res.extend(long[short.len()..].iter().cloned());
    res
}
/// Computes the convolution of a vector `a` such that it maps `a` to `b`.
/// This ensures that poly(a) * negative_poly(conjugate(a)) = poly(b),
/// where poly^(-1) are the coefficients of the polynomial.
/// Uses double-CRT representation for computing the convolution.
///
/// # Arguments
///
/// * `a_ring_el` - The input vector of `RingElement`.
///
/// # Returns
///
/// A vector of vectors containing elements of type `BASE_INT`.
// pub fn convolution(a_ring_el: &Vec<DPrimeRingElement>) -> Vec<DPrimeRingElement> {
//     let b_ring_el: Vec<DPrimeRingElement> = a_ring_el.iter().rev().map(|w| w.clone().conjugate()).collect();
//     let mut a_ring_el_copy = a_ring_el.clone();
//     let mut b_ring_el_copy = b_ring_el.clone();
//
//     let a: Vec<[u64; PHI]> = a_ring_el_copy
//         .iter_mut()
//         .map(|w| {
//             w.coeffs[0]
//         })
//         .collect();
//     let b: Vec<[u64; PHI]> = b_ring_el_copy
//         .iter_mut()
//         .map(|w| {
//             w.coeffs[0]
//         })
//         .collect();
//     let now = Instant::now();
//
//     let n = a[0].len().next_power_of_two() * 2;
//
//     let mut extended_a: Vec<Vec<u64>> = a
//         .par_iter()
//         .map(|t| {
//             let mut v = t.to_vec();
//             let mut v_64: Vec<u64> = v.iter().map(|t| *t as u64).collect();
//             v_64.resize(n, 0);
//             ntt_pow_of_2(&mut v_64, MOD_Q);
//             v_64
//         })
//         .collect();
//
//     let mut extended_b: Vec<Vec<u64>> = b
//         .par_iter()
//         .map(|t| {
//             let mut v = t.to_vec();
//             let mut v_64: Vec<u64> = v.iter().map(|t| *t as u64).collect();
//             v_64.resize(n, 0);
//             ntt_pow_of_2(&mut v_64, MOD_Q);
//             v_64
//         })
//         .collect();
//
//     let n2 = extended_a.len() * 2;
//     let n2_pow = n2.next_power_of_two();
//
//     let extended_c_transposed: Vec<Vec<u64>> = (0..extended_a[0].len()).into_par_iter().map(|t| {
//         let mut c_row = vec![0; extended_a.len()];
//         let mut at: Vec<u64> = extended_a.iter().map(|ai| ai[t] as u64).collect();
//         let mut bt: Vec<u64> = extended_b.iter().map(|bi| bi[t] as u64).collect();
//         at.resize(n2_pow, 0);
//         bt.resize(n2_pow, 0);
//
//         rayon::scope(|s| {
//             s.spawn(|_| {
//                 ntt_pow_of_2(&mut at, MOD_Q as u64);
//             });
//             s.spawn(|_| {
//                 ntt_pow_of_2(&mut bt, MOD_Q as u64);
//             });
//         });
//
//
//         let mut ct = hadamard_64(&at, &bt, MOD_Q as u64);
//
//         inverse_ntt_pow_of_2(&mut at, MOD_Q as u64);
//         inverse_ntt_pow_of_2(&mut bt, MOD_Q as u64);
//
//         reduce_mod_vec(&mut ct, MOD_Q as u64);
//         inverse_ntt_pow_of_2(&mut ct, MOD_Q as u64);
//
//             for (j, value) in ct.iter().enumerate().skip(extended_a.len() - 1).take(extended_a.len()) {
//                 c_row[j + 1 - extended_a.len()] = *value;
//             }
//             c_row.iter().map(|u| u.clone() as u64).collect()
//         }).collect();
//
//     transpose(&extended_c_transposed)
//         .into_par_iter()
//         .map(|mut v| {
//             inverse_ntt_pow_of_2(&mut v, MOD_Q as u64);
//             PrimeRing::new_from_larger_vec(&v)
//         })
//         .collect()
//
// }

pub fn convolution(a_ring_el: &Vec<DPrimeRingElement>) -> Vec<DPrimeRingElement> {
    // TODO. Typically it's balanced RNS so it could be faster then I guess.
    let b_ring_el: Vec<DPrimeRingElement> = a_ring_el.iter().rev().map(|w| w.clone().conjugate()).collect();

    // let result = (0..NOF_PRIMES).into_iter().map(|j| {
    let mut a_ring_el_copy = a_ring_el.clone();
    let mut b_ring_el_copy = b_ring_el.clone();

    let a: Vec<[u64; PHI]> = a_ring_el_copy
        .iter_mut()
        .map(|w| {
            w.coeffs
        })
        .collect();
    let b: Vec<[u64; PHI]> = b_ring_el_copy
        .iter_mut()
        .map(|w| {
            w.coeffs
        })
        .collect();
    let now = Instant::now();

    let n = a[0].len().next_power_of_two() * 2;

    let mut extended_a: Vec<Vec<u64>> = a
        .par_iter()
        .map(|t| {
            let mut v = t.to_vec();
            let mut v_64: Vec<u64> = v.iter().map(|t| *t as u64).collect();
            v_64.resize(n, 0);
            ntt_pow_of_2(&mut v_64, MOD_Q);
            v_64
        })
        .collect();

    let mut extended_b: Vec<Vec<u64>> = b
        .par_iter()
        .map(|t| {
            let mut v = t.to_vec();
            let mut v_64: Vec<u64> = v.iter().map(|t| *t as u64).collect();
            v_64.resize(n, 0);
            ntt_pow_of_2(&mut v_64, MOD_Q);
            v_64
        })
        .collect();

    let n2 = extended_a.len() * 2;
    let n2_pow = n2.next_power_of_two();

    let extended_c_transposed: Vec<Vec<u64>> = (0..extended_a[0].len()).into_par_iter().map(|t| {
        let mut c_row = vec![0; extended_a.len()];
        let mut at: Vec<u64> = extended_a.iter().map(|ai| ai[t] as u64).collect();
        let mut bt: Vec<u64> = extended_b.iter().map(|bi| bi[t] as u64).collect();
        at.resize(n2_pow, 0);
        bt.resize(n2_pow, 0);

        rayon::scope(|s| {
            s.spawn(|_| {
                ntt_pow_of_2(&mut at, MOD_Q as u64);
            });
            s.spawn(|_| {
                ntt_pow_of_2(&mut bt, MOD_Q as u64);
            });
        });


        let mut ct = hadamard_64(&at, &bt, MOD_Q as u64);

        inverse_ntt_pow_of_2(&mut at, MOD_Q as u64);
        inverse_ntt_pow_of_2(&mut bt, MOD_Q as u64);

        reduce_mod_vec(&mut ct, MOD_Q as u64);
        inverse_ntt_pow_of_2(&mut ct, MOD_Q as u64);

        for (j, value) in ct.iter().enumerate().skip(extended_a.len() - 1).take(extended_a.len()) {
            c_row[j + 1 - extended_a.len()] = *value;
        }
        c_row.iter().map(|u| u.clone() as u64).collect()
    }).collect();

    let result = transpose(&extended_c_transposed)
        .into_par_iter()
        .map(|mut v| {
            inverse_ntt_pow_of_2(&mut v, MOD_Q as u64);
            v
            // PrimeRing::new_from_larger_vec(&v)
        })
        .collect::<Vec<Vec<u64>>>();

        // result

    // }).collect::<Vec<Vec<Vec<u64>>>>();

    let mut results_elements = Vec::with_capacity(result[0].len());
    for i in 0 .. result.len() {
        results_elements.push(
            PrimeRing::new_from_larger_vec(&result[i])
        );
    }

    results_elements

}

#[test]
fn test_convolution() {
    let a = PrimeRing::random();
    let b = PrimeRing::random();
    let c = PrimeRing::random();
    let witness = vec![a, b, c];

    let convoluted = convolution(&witness);
    assert_eq!(convoluted.len(), 3);
    assert_eq!(convoluted, vec![
        a * a.conjugate() + b * b.conjugate() + c * c.conjugate(),
        c * b.conjugate() + b * a.conjugate(),
        c * a.conjugate()
    ]);
}
