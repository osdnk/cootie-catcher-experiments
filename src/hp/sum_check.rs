use num_traits::Zero;
use rayon::prelude::*;
use crate::arithmetic::{compute_hp_power_series, compute_one_prefixed_power_series, parallel_dot_matrix_matrix, parallel_dot_vector_matrix, ring_inner_product, row_wise_tensor, transpose};
use crate::helpers::println_with_timestamp;
use crate::hp::protocol_wrappers::prover::HPNorm2Output;
use crate::prime_ring::r#static::{NOF_THREADS_UPPED_BOUND, ONE, TWO};
use crate::prime_ring::ring::{DPrimeRingElement, PrimeRing};


// We assume the MLE in a form s.t. p(bin(k)) = w_k
fn partially_evaluate_hypercube(witness: &Vec<DPrimeRingElement>) -> Vec<DPrimeRingElement> {
    // partial evaluation as last arg = 0
    let even_witness: Vec<DPrimeRingElement> = (0..witness.len())
        .filter(|&i| i % 2 == 0)
        .map(|i| witness[i].clone())
        .collect();

    // partial evaluation as last arg = 1
    let odd_witness: Vec<DPrimeRingElement> = (0..witness.len())
        .filter(|&i| i % 2 != 0)
        .map(|i| witness[i].clone())
        .collect();

    even_witness.iter()
        .zip(odd_witness.iter())
        .map(|(a, b)| *a+ *b)
        .collect()
}

fn partially_evaluate_left(witness: &Vec<DPrimeRingElement>, c: &DPrimeRingElement) -> Vec<DPrimeRingElement> {
    // partial evaluation as last arg = 0
    let left_witness: Vec<DPrimeRingElement> = (0..witness.len())
        .filter(|&i| i < witness.len() / 2)
        .map(|i| witness[i].clone())
        .collect();

    // partial evaluation as last arg = 1
    let right_witness: Vec<DPrimeRingElement> = (0..witness.len())
        .filter(|&i| i >= witness.len() / 2)
        .map(|i| witness[i].clone())
        .collect();

    left_witness.iter()
        .zip(right_witness.iter())
        .map(|(l, r)| *l * (*ONE - *c) + *r * *c)
        .collect()
}

fn gen_row(ts: &Vec<DPrimeRingElement>) -> Vec<DPrimeRingElement> {
    let mut result = vec![*ONE];

    for t in ts.iter().rev() {
        let l_factor = *ONE - t.clone();
        let l_result = result
            .par_iter()
            .map(|m| *m * l_factor)
            .collect::<Vec<DPrimeRingElement>>();
        let r_factor = t.clone();
        let r_result = result
            .par_iter()
            .map(|m| *m * r_factor)
            .collect::<Vec<DPrimeRingElement>>();
        result = [l_result, r_result].concat();
    }
    // create a row in a form (ts_0 ts_0) \tensor (ts_1 ts_1) ...
    result
}

// 1 - x, x
fn sum_check(witness: &Vec<DPrimeRingElement>)-> (DPrimeRingElement, DPrimeRingElement, Vec<DPrimeRingElement>) {
    let mut precomputed = Vec::new();
    precomputed.push(witness.clone()); // TODO
    let mut temp_wit = partially_evaluate_hypercube(&witness);
    while temp_wit.len() != 1 {
        let new_temp_wit = partially_evaluate_hypercube(&temp_wit);
        precomputed.push(temp_wit);
        temp_wit = new_temp_wit;
    }

    let c = *temp_wit.last().unwrap();
    let mut cc = c.clone();

    let mut ts = Vec::with_capacity(witness.len().ilog2() as usize);
    for i in 0..witness.len().ilog2() {
        let p_i = precomputed.get(precomputed.len() - 1 - i as usize).unwrap();
        let c_i_0 = p_i.get(0).unwrap();
        let c_i_1 = p_i.get(1).unwrap();
        assert_eq!(*c_i_0 + *c_i_1, cc);
        let t_i = PrimeRing::random();
        cc = *c_i_0 * (*ONE - t_i) + *c_i_1 * t_i;
        precomputed = precomputed.iter().map(|row|
            if row.len() == 1 {
                row.clone()
            } else {
                partially_evaluate_left(&row, &t_i)
            }
        ).collect();
        ts.push(t_i);
    }

    let expanded_row = gen_row(&ts);
    assert_eq!(ring_inner_product(
        &expanded_row,
        &witness
    ), cc);
    (c, cc, expanded_row)
}


#[test]
fn test_partially_evaluate_hypercube() {
    let a = PrimeRing::random();
    let b = PrimeRing::random();
    let c = PrimeRing::random();
    let d = PrimeRing::random();
    let witness = [a,b,c,d].to_vec();
    let new_wit = partially_evaluate_hypercube(&witness);
    assert_eq!([a + b, c + d].to_vec(), new_wit)
}

#[test]
fn test_sum_check() {
    let a = PrimeRing::random();
    let b = PrimeRing::random();
    let c = PrimeRing::random();
    let d = PrimeRing::random();
    let e = PrimeRing::random();
    let f = PrimeRing::random();
    let g = PrimeRing::random();
    let h = PrimeRing::random();
    let witness = [a,b,c,d,e,f,g,h].to_vec();
    let sum = sum_check(&witness).0;
    assert_eq!(a + b + c + d + e + f + g + h, sum)
}


// fn evaluate_poly(coefs: &Vec<DPrimeRingElement>, point: &DPrimeRingElement) -> DPrimeRingElement {
//     coefs.iter().enumerate().map(|(i, c)| *c * point.clone().pow(i as u32)).sum()
// }
// fn evaluate_poly(coefs: &Vec<DPrimeRingElement>, point: &DPrimeRingElement) -> DPrimeRingElement {
//     coefs.iter().enumerate().map(|(i, c)| c.clone() * point.clone().pow(i as u32)).sum()
// }

fn evaluate_poly(coefs: &Vec<DPrimeRingElement>, point: &DPrimeRingElement) -> DPrimeRingElement {
    let mut result = PrimeRing::constant(0);
    let mut temp = *ONE; // temp will hold point^i during each iteration
    for coef in coefs {
        result = result + coef.clone() * temp.clone();
        temp = temp * point.clone(); // update temp to point^(i+1)
    }
    result
}

// fn poly_from_coefs(coefs: &Vec<DPrimeRingElement>) -> Vec<DPrimeRingElement> {
//     coefs.iter().enumerate().map(|(i, c)| *c * PrimeRing::constant(i as u64)).collect()
// }

fn sq_sum_check(witness: &Vec<DPrimeRingElement>) -> DPrimeRingElement {
    let mut n = witness.len();
    let mut curr_c = witness.iter().map(|w| *w * *w).sum();
    let mut curr_witness = witness.clone();
    let mut ts = Vec::new();

    while curr_witness.len() > 1 {
        let half = n / 2;
        let left = curr_witness.iter().take(half).cloned().collect::<Vec<_>>();
        let right = curr_witness.iter().skip(half).take(half).cloned().collect::<Vec<_>>();

        // Calculate the polynomial coefficients for p(t) where p(t) = sum(((1 - t) * left[i] + t * right[i])**2 for i in range(half))
        let mut p_coefs = vec![DPrimeRingElement::zero(); 3];
        for i in 0..half {
            let a = left[i].clone();
            let b = right[i].clone();
            p_coefs[0] = p_coefs[0] + a * a;
            p_coefs[1] = p_coefs[1] + *TWO * a * b - *TWO * a * a;
            p_coefs[2] = p_coefs[2] + a * a - *TWO * a * b + b * b;
        }

        // Evaluate p(t) at t=0 and t=1
        let p0 = evaluate_poly(&p_coefs, &PrimeRing::constant(0));
        let p1 = evaluate_poly(&p_coefs, &*ONE);
        assert_eq!(p0 + p1, curr_c, "Verification failed");

        // Generate random challenge t_i
        let t_i = PrimeRing::random();
        ts.push(t_i.clone());

        // Compute the next claim value
        curr_c = evaluate_poly(&p_coefs, &t_i);

        // Update the current witness
        curr_witness = (0..half).map(|i| {
            let left_i = left[i].clone();
            let right_i = right[i].clone();
            (*ONE - t_i.clone()) * left_i + t_i.clone() * right_i
        }).collect();

        n = half;
    }

    // Generate the expanded row and verify the final condition
    let expanded_row = gen_row(&ts);
    let ip = ring_inner_product(&expanded_row, &witness);
    assert_eq!(ip * ip, curr_c);
    witness.iter().map(|w| *w * *w).sum()
}



#[test]
fn test_s_sum_check_varying_exponents() {
    fn test_s_sum_check(exponent: u32) {
        // Generate a witness of length 2^exponent
        let n = 2_u32.pow(exponent) as usize;
        let witness = (0..n).map(|_| PrimeRing::random()).collect::<Vec<_>>();
        let expected_sum = witness.iter().map(|w| *w * *w).sum::<DPrimeRingElement>();
        let s_sum_result = sq_sum_check(&witness);
        assert_eq!(s_sum_result, expected_sum, "Test failed for witness length {}", n);
    }
    test_s_sum_check(3);
    test_s_sum_check(5);
}

fn verifier_sum_check() {

}

fn sq_sum_check_batched_old(witness: &Vec<Vec<DPrimeRingElement>>) -> Vec<DPrimeRingElement> {
    let mut n = witness.len();
    let columns = transpose(&witness);
    let mut res = vec![DPrimeRingElement::zero(); columns.len()];
    let mut claims :Vec<DPrimeRingElement> = columns
        .iter()
        .map(|column| column.iter().map(|w| *w * *w).sum())
        .collect();

    let mut current_witnesses = columns.clone();

    let mut ts = Vec::new();
    while current_witnesses.first().map_or(0, |w| w.len()) > 1 {
        let half = n / 2;
        let t_i = PrimeRing::random(); // Generate random challenge
        ts.push(t_i.clone());

        claims = current_witnesses.par_iter_mut().enumerate().map(|(i, curr_witness)| {
            let mut p_coefs = vec![DPrimeRingElement::zero(); 3];
            for j in 0..half {
                let a = curr_witness[j];
                let b = curr_witness[j + half];
                p_coefs[0] = p_coefs[0] + a * a;
                p_coefs[1] = p_coefs[1] + *TWO * a * b - *TWO * a * a;
                p_coefs[2] = p_coefs[2] + a * a - *TWO * a * b + b * b;
            }

            // Evaluate p(t) at t=0 and t=1
            let p0 = evaluate_poly(&p_coefs, &PrimeRing::constant(0));
            let p1 = evaluate_poly(&p_coefs, &*ONE);
            assert_eq!(p0 + p1, claims[i], "Verification failed");

            let new_claim = evaluate_poly(&p_coefs, &t_i);

            *curr_witness = (0..half).map(|j| {
                let left_j = curr_witness[j];
                let right_j = curr_witness[j + half];
                (*ONE - t_i) * left_j + t_i * right_j
            }).collect();
            new_claim
        }).collect();

        n = half;
    }

    let expanded_row = gen_row(&ts);
    for (i, column) in columns.iter().enumerate() {
        let ip = ring_inner_product(&expanded_row, column);
        assert_eq!(ip * ip, claims[i]);
        res[i] = column.iter().map(|w| *w * *w).sum();
    }
    columns
        .iter()
        .map(|column| column.iter().map(|w| *w * *w).sum())
        .collect()
}

static mut g_claim :Vec<DPrimeRingElement> = Vec::new();
static mut g_claim_curr :Vec<DPrimeRingElement> = Vec::new();
pub static mut g_ts :Vec<DPrimeRingElement> = Vec::new();
static mut g_last_claims:Vec<DPrimeRingElement> = Vec::new();
pub fn receive_initial_claim(claim: &Vec<DPrimeRingElement>) {
    unsafe {
        g_claim = claim.clone();
        g_claim_curr = claim.clone();
        g_ts.clear();
    }
}

// Function to verify polynomial claims and return the challenge
pub fn verify_polynomials_and_generate_challenge(
    p_coefs_all: &Vec<Vec<DPrimeRingElement>>,
) -> DPrimeRingElement {
    // Generate random challenge t_i
    // sample from real so only oce check needs to be passed
    let t_i = PrimeRing::random_real();

    p_coefs_all.par_iter().enumerate().for_each(|(index, p_coefs)| {
        let p0 = evaluate_poly(p_coefs, &PrimeRing::constant(0));
        let p1 = evaluate_poly(p_coefs, &*ONE);
        unsafe {
            assert_eq!(p0 + p1, *g_claim_curr.get(index).unwrap(), "Verification failed");
            g_claim_curr[index] = evaluate_poly(p_coefs, &t_i);
        }
    });


    unsafe {
        g_ts.push(t_i.clone());
    }
    t_i
}

pub fn receive_unsquared_claims(unsquared_claims: &Vec<DPrimeRingElement>, unsquared_claims_2: &Vec<DPrimeRingElement>) {
    unsafe {
        unsquared_claims.par_iter().zip(unsquared_claims_2.par_iter()).zip(g_last_claims.par_iter()).for_each(|((a, b), c)| { assert_eq!(*a * *b, *c) });
    }
}

pub fn sq_sum_check_batched_init(witness: &Vec<Vec<DPrimeRingElement>>, witness_2: &Vec<Vec<DPrimeRingElement>>) -> Vec<DPrimeRingElement> {
    let columns = transpose(&witness);
    let columns_2 = transpose(&witness_2);
    columns
        .par_iter()
        .zip(columns_2.par_iter())
        .map(|(column, columm_2)| column
            .par_iter()
            .zip(columm_2.par_iter())
            .with_min_len(columns[0].len() / NOF_THREADS_UPPED_BOUND * columns.len())
            .map(|(w, w_2)| *w * *w_2).sum()
        )
        .collect()
}

pub fn sq_sum_check_batched_par(
    witness: &Vec<Vec<DPrimeRingElement>>,
    witness_2: &Vec<Vec<DPrimeRingElement>>,
    mut verify: impl FnMut(&Vec<Vec<DPrimeRingElement>>) -> DPrimeRingElement
) -> (Vec<DPrimeRingElement>, Vec<DPrimeRingElement>) {
    let mut n = witness.len();
    let columns = transpose(&witness);
    let columns_2 = transpose(&witness_2);
    let mut claims = unsafe { g_claim .clone() };



    let mut current_witnesses = columns;
    let mut current_witnesses_2 = columns_2;
    let mut ts = Vec::new();

    while current_witnesses.first().map_or(0, |w| w.len()) > 1 {
        let half = n / 2;

        let p_coefs_all = current_witnesses.par_iter().zip(current_witnesses_2.par_iter()).map(|(curr_witness, curr_witness_2)|{
            let p_coefs = (0..half)
                .into_par_iter()
                .with_min_len(current_witnesses[0].len() * current_witnesses.len() / NOF_THREADS_UPPED_BOUND)
                .map(|j| {
                    let a = curr_witness[j];
                    let b = curr_witness[j + half];
                    let a_2 = curr_witness_2[j];
                    let b_2 = curr_witness_2[j + half];
                    let p0 = a * a_2;
                    let p1 = a * b_2 + a_2 * b -  a * a_2 - a * a_2;
                    let p2 = a * a_2 - a_2 * b - a * b_2 + b * b_2;
                (p0, p1, p2)
            }).reduce(
                || (DPrimeRingElement::zero(), DPrimeRingElement::zero(), DPrimeRingElement::zero()),
                |(p0_acc, p1_acc, p2_acc), (p0, p1, p2)| {
                    (
                        p0_acc + p0,
                        p1_acc + p1,
                        p2_acc + p2,
                    )
                }
            );
            vec![p_coefs.0,p_coefs.1,p_coefs.2]
        }).collect();

        // Verify all polynomials and generate a single challenge
        let t_i = verify(&p_coefs_all);
        ts.push(t_i);

        let current_witnesses_len = current_witnesses.len();
        current_witnesses
            .par_iter_mut()
            .enumerate()
            .for_each(|(i, curr_witness)| {
                let half = curr_witness.len() / 2;
                let new_witness: Vec<_> = (0..half)
                    .into_par_iter()
                    .with_min_len(curr_witness.len() * current_witnesses_len / NOF_THREADS_UPPED_BOUND)
                    .map(|j| {
                        let left_j = &curr_witness[j];
                        let right_j = &curr_witness[j + half];
                        (*ONE - t_i) * *left_j + t_i * *right_j
                    })
                    .collect();
                *curr_witness = new_witness;
            });

        current_witnesses_2
            .par_iter_mut()
            .enumerate()
            .for_each(|(i, curr_witness)| {
                let half = curr_witness.len() / 2;
                let new_witness: Vec<_> = (0..half)
                    .into_par_iter()
                    .with_min_len(curr_witness.len() * current_witnesses_len / NOF_THREADS_UPPED_BOUND)
                    .map(|j| {
                        let left_j = &curr_witness[j];
                        let right_j = &curr_witness[j + half];
                        (*ONE - t_i) * *left_j + t_i * *right_j
                    })
                    .collect();
                *curr_witness = new_witness;
            });

        // Step 2: Compute claims from updated witnesses
        claims = (0..current_witnesses_len)
            .into_par_iter()
            .map(|i| {
                let p_coefs = &p_coefs_all[i];
                evaluate_poly(p_coefs, &t_i)
            })
            .collect();
        n = half;

    }

    let expanded_row = gen_row(&ts);
    let unsquared_claims  = parallel_dot_vector_matrix(&expanded_row, &witness);
    let unsquared_claims_2  = parallel_dot_vector_matrix(&expanded_row, &witness_2);
    (unsquared_claims, unsquared_claims_2)
}

#[test]
fn test_sq_sum_check_batched() {
    let witness = vec![
        vec![PrimeRing::random(), PrimeRing::random()],
        vec![PrimeRing::random(), PrimeRing::random()],
        vec![PrimeRing::random(), PrimeRing::random()],
        vec![PrimeRing::random(), PrimeRing::random()],
    ];

    let witness_2 = vec![
        vec![PrimeRing::random(), PrimeRing::random()],
        vec![PrimeRing::random(), PrimeRing::random()],
        vec![PrimeRing::random(), PrimeRing::random()],
        vec![PrimeRing::random(), PrimeRing::random()],
    ];

    // let witness_2 = witness.clone();

    let transposed_witness = transpose(&witness);
    let transposed_witness_2 = transpose(&witness_2);

    let expected_sums: Vec<DPrimeRingElement> = transposed_witness
        .iter().zip(transposed_witness_2.iter())
        .map(|(column, column_2)| column.iter().zip(column_2.iter()).map(|(w, w_2)| *w * *w_2).sum())
        .collect();

    let claims = sq_sum_check_batched_init(&witness, &witness_2);
    receive_initial_claim(&claims);
    let (unsquared_claims, unsquared_claims_2) = sq_sum_check_batched_par(
        &witness,
        &witness_2,
        verify_polynomials_and_generate_challenge
    );
    receive_unsquared_claims(&unsquared_claims, &unsquared_claims_2);

    for (expected_sum, s_sum_result) in expected_sums.iter().zip(claims.iter()) {
        assert_eq!(expected_sum, s_sum_result, "Batched Test failed");
    }
}

