use std::ops::Index;
use std::time::Instant;
use crate::arithmetic::{add_matrices, compose_with_radix_mod, fast_power, first_n_columns, join_matrices_horizontally, last_n_columns, parallel_dot_matrix_matrix, PowerSeries, ring_inner_product, row_wise_tensor, sample_random_ss_mat, transpose};
use crate::helpers::println_with_timestamp;
use crate::prime_ring::r#static::{MOD_Q, REP};
use crate::prime_ring::ring::{DPrimeRingElement, PrimeRing, PrimeRingElement};
use crate::subroutines::decomp::DecompOutput;
use crate::subroutines::norm_first_round::Norm1Output;
use crate::subroutines::norm_second_round::Norm2Output;
use crate::subroutines::split::SplitOutput;

pub(crate) struct VerifierState {
    pub(crate) wit_cols: usize,
    pub(crate) wit_rows: usize,
    pub(crate) rhs: Vec<Vec<DPrimeRingElement>>
}

pub fn verifier_split(
    power_series: &Vec<PowerSeries>,
    output_split: SplitOutput,
    verifier_state: &VerifierState,
) -> VerifierState {

    let now = Instant::now();
    let ck_l = first_n_columns(&output_split.rhs, verifier_state.wit_cols);
    let ck_r = last_n_columns(&output_split.rhs, verifier_state.wit_cols);
    let elapsed = now.elapsed();
    println_with_timestamp!("  Time to extract left and right columns from RHS: {:.2?}", elapsed);
    // Compute the multiplier and row-wise tensor product
    let now = Instant::now();

    let mut multiplier_l = Vec::with_capacity(power_series.len());
    let mut multiplier_r = Vec::with_capacity(power_series.len());

    power_series.iter().for_each(|ps| {
        let index = ps.expanded_layers.iter().position(|el| { el.len() == verifier_state.wit_rows }).unwrap();
        multiplier_l.push(ps.tensors[index][0]);
        multiplier_r.push(ps.tensors[index][1]);
    });

    let elapsed = now.elapsed();
    println_with_timestamp!("  Time to compute the multiplier: {:.2?}", elapsed);

    let now = Instant::now();
    let ck_l_multiplied = row_wise_tensor(&ck_l, &transpose(&vec![multiplier_l]));
    let ck_r_multiplied = row_wise_tensor(&ck_r, &transpose(&vec![multiplier_r]));
    let elapsed = now.elapsed();
    println_with_timestamp!("  Time for row-wise tensor product: {:.2?}", elapsed);

    let now = Instant::now();
    let ck_lr = add_matrices(&ck_r_multiplied, &ck_l_multiplied);
    assert_eq!(ck_lr, verifier_state.rhs);

    let elapsed = now.elapsed();
    println_with_timestamp!("  Time verification: {:.2?}", elapsed);

    VerifierState {
        wit_cols: verifier_state.wit_cols * 2,
        wit_rows: verifier_state.wit_rows / 2,
        rhs: output_split.rhs
    }
}
pub fn verifier_fold(verifier_state: &VerifierState, challenge: &Vec<Vec<DPrimeRingElement>>) -> VerifierState {
    VerifierState {
        wit_cols: REP,
        wit_rows: verifier_state.wit_rows,
        rhs: parallel_dot_matrix_matrix(&verifier_state.rhs, &challenge)
    }
}

pub fn challenge_for_fold(verifier_state: &VerifierState) -> Vec<Vec<DPrimeRingElement>> {
    sample_random_ss_mat(verifier_state.wit_cols, REP)
}

pub fn norm_challenge(norm_1_output: &Norm1Output, verifier_state: &VerifierState) -> (VerifierState, DPrimeRingElement, DPrimeRingElement) {
    let challenge = PrimeRing::random();
    let inverse_challenge = challenge.inverse().conjugate();
    let state = VerifierState {
        wit_cols: verifier_state.wit_cols * 3,
        wit_rows: verifier_state.wit_rows,
        rhs: join_matrices_horizontally(&verifier_state.rhs, &norm_1_output.new_rhs.clone())
    };
    (state, challenge, inverse_challenge)
}

pub fn verify_norm_2(norm_1_output: &Norm1Output, norm_2_output: &Norm2Output, state: &VerifierState) -> VerifierState {
    let new_evaluations = compose_with_radix_mod(
        &last_n_columns(&norm_2_output.new_rhs, state.wit_cols / 3 * 2),
        norm_1_output.radix,
        2,
    );

    for i in 0..state.wit_cols / 3 {
        let ip = new_evaluations[2][i];
        assert_eq!(ip, new_evaluations[2][i]);
        // FIXME check something about the inner product
        assert_eq!(
            new_evaluations[0][i] + new_evaluations[1][i].conjugate(),
            norm_2_output.new_rhs[0][i] * norm_2_output.new_rhs[1][i].conjugate() + ip
        );
    }

    println_with_timestamp!("{:?} {:?} {:?} {:?}", state.rhs.len(),state.rhs[0].len(), norm_2_output.new_rhs.len(), norm_2_output.new_rhs[0].len());
    VerifierState {
        wit_cols: state.wit_cols,
        wit_rows: state.wit_rows,
        rhs: vec![state.rhs.clone(), norm_2_output.new_rhs.clone()].concat()
    }
}


pub fn verify_decomp(decomp_output: DecompOutput, verifier_state: &VerifierState) -> VerifierState {
    let composed_rhs = compose_with_radix_mod(&decomp_output.rhs, decomp_output.radix, decomp_output.parts);

    assert_eq!(composed_rhs, verifier_state.rhs);

    VerifierState {
        wit_rows: verifier_state.wit_rows,
        wit_cols: verifier_state.wit_cols * decomp_output.parts,
        rhs: decomp_output.rhs,
    }

}
