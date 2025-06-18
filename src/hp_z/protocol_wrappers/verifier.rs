use rayon::iter::IntoParallelRefIterator;
use crate::arithmetic::{random, ring_inner_product};
use crate::hp_z::protocol_wrappers::prover::{HPNorm1Output, HPNorm3Output};
use crate::hp_z::protocol_wrappers::prover::HPNorm2Output;
use crate::hp_z::sum_check_z::{receive_initial_claim, receive_unsquared_claim, verify_polynomials_and_generate_challenge};
use crate::hp_z::zq::{Zq, ZqElement};
use crate::prime_ring::r#static::MOD_Q;
use crate::prime_ring::ring::{DPrimeRingElement, PrimeRing};
use crate::subroutines::norm_first_round::Norm1Output;
use crate::subroutines::norm_second_round::Norm2Output;
use crate::subroutines::verifier::VerifierState;
use rayon::prelude::*;


pub fn verify_norm_1(norm_1_output: &HPNorm1Output) {
    let challenges = random(norm_1_output.claims.len(), MOD_Q).iter().map(
        |p| Zq::from_u64(*p)
    ).collect::<Vec<ZqElement>>();

    let claim = challenges.par_iter()
        .zip(norm_1_output.claims.par_iter())
        .map(|(a_i, b_i)| *a_i * *b_i)
        .reduce(Zq::zero, |acc, prod| acc + prod);


    receive_initial_claim(&claim);

    challenges
}


pub fn verify_norm_2(norm_2_output: &HPNorm2Output) -> DPrimeRingElement {
    verify_polynomials_and_generate_challenge(&norm_2_output.p_coefs_all)
}

pub fn verify_norm_3(norm_3_output: &HPNorm3Output, state: &VerifierState) -> VerifierState {
    receive_unsquared_claims(&norm_3_output.unsquared_claims);
    VerifierState {
        wit_cols: state.wit_cols,
        wit_rows: state.wit_rows,
        rhs: vec![state.rhs.clone(), vec![norm_3_output.unsquared_claims.clone()]].concat()
    }
}

