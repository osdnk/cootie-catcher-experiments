use crate::arithmetic::{random, ring_inner_product, sample_random_vector};
use crate::hp::protocol_wrappers::prover::{HPNorm1Output, HPNorm3Output};
use crate::hp::protocol_wrappers::prover::HPNorm2Output;
use crate::hp::sum_check::{g_ts, receive_initial_claim, receive_unsquared_claims, verify_polynomials_and_generate_challenge};
use crate::prime_ring::r#static::MOD_Q;
use crate::prime_ring::ring::DPrimeRingElement;
use crate::subroutines::norm_first_round::Norm1Output;
use crate::subroutines::norm_second_round::Norm2Output;
use crate::subroutines::verifier::VerifierState;
use rayon::prelude::*;

pub fn verify_norm_1(norm_1_output: &HPNorm1Output) {
    receive_initial_claim(&norm_1_output.claims);
    // TODO check norms
}


pub fn verify_norm_2(norm_2_output: &HPNorm2Output) -> DPrimeRingElement {
    verify_polynomials_and_generate_challenge(&norm_2_output.p_coefs_all)
}

pub fn verify_norm_3(norm_3_output: &HPNorm3Output, state: &VerifierState) -> VerifierState {
    // technically those should be batched with u
    receive_unsquared_claims(&norm_3_output.unsquared_claims, &norm_3_output.unsquared_claims_conj);
    norm_3_output
        .unsquared_claims
        .par_iter()
        .zip(norm_3_output.unsquared_claims_conj.par_iter())
        .for_each(|(a, b)| { assert_eq!(*a, b.conjugate()) });


    VerifierState {
        wit_cols: state.wit_cols,
        wit_rows: state.wit_rows,
        // rhs: state.rhs.clone(),
        rhs: vec![
            state.rhs.clone(),
            vec![
                norm_3_output.unsquared_claims.clone(),
            ]
        ].concat()
    }
}

