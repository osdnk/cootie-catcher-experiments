use crate::arithmetic::{compute_hp_power_series, conjugate_vector, parallel_dot_matrix_vector, parallel_dot_vector_matrix, transpose, PowerSeries};
use crate::hp::sum_check::{g_ts, sq_sum_check_batched_init, sq_sum_check_batched_par, verify_polynomials_and_generate_challenge};
use crate::prime_ring::ring::DPrimeRingElement;
use rayon::prelude::*;

pub struct HPNorm1Output {
    pub claims: Vec<DPrimeRingElement>
}

pub struct HPNorm2Output {
    pub p_coefs_all: Vec<Vec<DPrimeRingElement>>
}

pub struct HPNorm3Output {
    pub unsquared_claims: Vec<DPrimeRingElement>,
    pub unsquared_claims_conj: Vec<DPrimeRingElement>
}
pub fn norm_1(
    witness: &Vec<Vec<DPrimeRingElement>>,
    // challenges: Vec<DPrimeRingElement>,
) -> (HPNorm1Output, Vec<Vec<DPrimeRingElement>>) {
    // let squashed_wit = parallel_dot_vector_matrix(&challenges, &witness);
    let conj_witness = witness.iter().map(|r| {
        r.iter().map(|e| e.conjugate()).collect()
    }).collect();

    (HPNorm1Output {
        claims: sq_sum_check_batched_init(&witness, &conj_witness)
    }, conj_witness)
}

pub fn norm_2(
    power_series: &Vec<PowerSeries>,
    witness: &Vec<Vec<DPrimeRingElement>>,
    conj_witness: &Vec<Vec<DPrimeRingElement>>,
    mut verify_norm_2_wrapped: impl FnMut(&HPNorm2Output) -> DPrimeRingElement
) -> (Vec<PowerSeries>, HPNorm3Output) {


    let verify_norm_2 = |p_coefs_all: &Vec<Vec<DPrimeRingElement>>| {
        let output = HPNorm2Output {
            p_coefs_all: p_coefs_all.clone()
        };
        verify_norm_2_wrapped(&output)
    };

    let (unsquared_claims, unsquared_claims_conj) = sq_sum_check_batched_par(
        &witness,
        &conj_witness,
        verify_norm_2
    );


    unsafe {
        (vec![power_series.clone(),
              vec![
                  compute_hp_power_series(&g_ts),
              ]].concat(),
         HPNorm3Output {
             unsquared_claims,
             unsquared_claims_conj
         })
    }
}

