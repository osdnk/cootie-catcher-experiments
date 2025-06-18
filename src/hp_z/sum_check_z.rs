use num_traits::Zero;
use rayon::prelude::*;
use crate::hp_z::zq::{Zq, ZqElement};
use crate::prime_ring::r#static::TWO;

fn gen_row(ts: &Vec<ZqElement>) -> Vec<ZqElement> {
    let mut result = vec![Zq::one()];
    for t in ts.iter().rev() {
        let l_factor = Zq::one() - *t;
        let r_factor = *t;
        let l_result: Vec<_> = result.par_iter().map(|m| *m * l_factor).collect();
        let r_result: Vec<_> = result.par_iter().map(|m| *m * r_factor).collect();
        result = [l_result, r_result].concat();
    }
    result
}

fn evaluate_poly(coefs: &Vec<ZqElement>, point: &ZqElement) -> ZqElement {
    let mut result = Zq::zero();
    let mut temp = Zq::one();
    for coef in coefs {
        result = result + *coef * temp;
        temp = temp * *point;
    }
    result
}

static mut G_CLAIM: Option<ZqElement> = None;
static mut G_CLAIM_CURR: Option<ZqElement> = None;
pub static mut G_TS: Vec<ZqElement> = Vec::new();
static mut G_LAST_CLAIM: Option<ZqElement> = None;

pub fn receive_initial_claim(claim: ZqElement) {
    unsafe {
        G_CLAIM = Some(claim);
        G_CLAIM_CURR = Some(claim);
        G_TS.clear();
    }
}

pub fn receive_unsquared_claim(unsquared_claim: ZqElement, unsquared_claim_2: ZqElement) {
    unsafe {
        match G_LAST_CLAIM {
            Some(val) => assert_eq!(unsquared_claim * unsquared_claim_2, val),
            None => panic!("G_LAST_CLAIM not set before calling receive_unsquared_claim"),
        }
    }
}

pub fn verify_polynomial_and_generate_challenge(p_coefs: &Vec<ZqElement>) -> ZqElement {
    let t_i = Zq::random();
    let p0 = evaluate_poly(p_coefs, &Zq::zero());
    let p1 = evaluate_poly(p_coefs, &Zq::one());
    unsafe {
        match G_CLAIM_CURR {
            Some(curr) => {
                assert_eq!(p0 + p1, curr, "Verification failed");
                G_CLAIM_CURR = Some(evaluate_poly(p_coefs, &t_i));
                G_TS.push(t_i);
            }
            None => panic!("G_CLAIM_CURR not set before verification"),
        }
    }
    t_i
}

pub fn sq_sum_check_single_init(witness: &Vec<ZqElement>, witness_2: &Vec<ZqElement>) -> ZqElement {
    witness.iter().zip(witness_2.iter()).map(|(a, b)| *a * *b).sum()
}

pub fn sq_sum_check_single_par(
    witness: &Vec<ZqElement>,
    witness_2: &Vec<ZqElement>,
    mut verify: impl FnMut(&Vec<ZqElement>) -> ZqElement,
) -> (ZqElement, ZqElement) {
    let mut n = witness.len();
    let mut current = witness.clone();
    let mut current_2 = witness_2.clone();
    let mut ts = Vec::new();
    let mut claim = unsafe { G_CLAIM.expect("G_CLAIM not initialized") };

    while current.len() > 1 {
        let half = n / 2;

        // Calculate the polynomial coefficients for p(t) where p(t) = sum(((1 - t) * left[i] + t * right[i]) * ((1 - t) * left_2[i] + t * right_2[i])  for i in range(half))
        let (p0, p1, p2) = (0..half)
            .into_par_iter()
            .map(|j| {
                let a = current[j];
                let b = current[j + half];
                let a_2 = current_2[j];
                let b_2 = current_2[j + half];
                let p0 = a * a_2;
                let p1 = a * b_2 + a_2 * b -  a * a_2 - a * a_2;
                let p2 = a * a_2 - a_2 * b - a * b_2 + b * b_2;
                (p0, p1, p2)
            })
            .reduce(
                || (ZqElement::zero(), ZqElement::zero(), ZqElement::zero()),
                |(acc0, acc1, acc2), (x0, x1, x2)| (acc0 + x0, acc1 + x1, acc2 + x2),
            );

        let poly = vec![p0, p1, p2];
        let t_i = verify(&poly);
        ts.push(t_i);

        current = (0..half)
            .into_par_iter()
            .map(|j| {
                let a = current[j];
                let b = current[j + half];
                (Zq::one() - t_i) * a + t_i * b
            })
            .collect();


        current_2 = (0..half)
            .into_par_iter()
            .map(|j| {
                let a = current_2[j];
                let b = current_2[j + half];
                (Zq::one() - t_i) * a + t_i * b
            })
            .collect();

        claim = evaluate_poly(&poly, &t_i);
        n = half;
    }

    let expanded_row = gen_row(&ts);
    let unsquared_claim: ZqElement = expanded_row
        .par_iter()
        .zip(witness.par_iter())
        .map(|(r, w)| *r * *w)
        .sum();

    let unsquared_claim_2: ZqElement = expanded_row
        .par_iter()
        .zip(witness_2.par_iter())
        .map(|(r, w)| *r * *w)
        .sum();

    unsafe {
        G_LAST_CLAIM = Some(claim);
    }

    (unsquared_claim, unsquared_claim_2)
}

#[test]
fn test_sq_sum_check_single() {
    let witness: Vec<ZqElement> = (0..8).map(|_| Zq::random()).collect();
    let witness_2: Vec<ZqElement> = (0..8).map(|_| Zq::random()).collect();
    // let witness_2: Vec<ZqElement> = witness.clone();
    let expected_sum: ZqElement = witness.iter().zip(witness_2.iter()).map(|(a, b)| *a * *b).sum();

    let claim = sq_sum_check_single_init(&witness, &witness_2);
    receive_initial_claim(claim);

    let (unsquared_claim, unsquared_claim_2) = sq_sum_check_single_par(&witness, &witness_2, verify_polynomial_and_generate_challenge);
    receive_unsquared_claim(unsquared_claim, unsquared_claim_2);

    assert_eq!(expected_sum, claim, "Single-vector sum-check failed");
}
