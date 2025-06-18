use std::time::Instant;
use crate::arithmetic::{parallel_dot_series_matrix, sample_bin_random_mat, sample_short_random_mat};
use crate::helpers::println_with_timestamp;
use crate::hp;
use crate::hp::protocol_wrappers::prover::HPNorm2Output;
use crate::prime_ring::r#static::{DECOMP_PARTS, MODULE_SIZE, NOF_ROUNDS, INIT_REP, SUM_CHECK, WIT_DIM};
use crate::prime_ring::ring::{DPrimeRingElement};
use crate::subroutines::crs::CRS;
use crate::subroutines::decomp::{decomp, decomp_ell};
use crate::subroutines::fold::fold;
use crate::subroutines::norm_first_round::norm_1;
use crate::subroutines::norm_second_round::norm_2;
use crate::subroutines::split::split;
use crate::subroutines::verifier::{VerifierState, challenge_for_fold, norm_challenge, verifier_fold, verifier_split, verify_decomp, verify_norm_2};

pub fn protocol() {
    println_with_timestamp!("Start CRS");

    let crs = CRS::gen_crs(WIT_DIM, MODULE_SIZE);
    println_with_timestamp!("end CRS");
    let mut witness = sample_bin_random_mat(WIT_DIM, INIT_REP);
    println_with_timestamp!("end sampling witness");


    let now = Instant::now();
    let commitment = parallel_dot_series_matrix(&crs.ck, &witness);
    let elapsed = now.elapsed();

    let mut verifier_state = VerifierState {
        wit_cols: INIT_REP,
        wit_rows: WIT_DIM,
        rhs: commitment
    };

    let mut verifier_runtime = Instant::now().elapsed();
    let mut prover_runtime = Instant::now().elapsed();

    println_with_timestamp!("Time for parallel_dot_matrix_matrix (commitment): {:.2?}", elapsed);
    prover_runtime = prover_runtime + elapsed;


    let mut statement = crs.ck; // TODO -- will be more things
    for i in 0..NOF_ROUNDS {

        let parts = DECOMP_PARTS[i as usize];

        if parts > 1 {
            let now = Instant::now();
            let (new_witness, bdecomp_output) = decomp_ell(&statement, &witness, parts);
            witness = new_witness;
            let elapsed = now.elapsed();
            println_with_timestamp!("Time for decomp: {:.2?}", elapsed);
            prover_runtime = prover_runtime + elapsed;


            let now = Instant::now();
            let new_verifier_state = verify_decomp(bdecomp_output, &verifier_state);
            verifier_state = new_verifier_state;
            let elapsed = now.elapsed();
            println_with_timestamp!("Time for verify_decomp: {:.2?}", elapsed);
            verifier_runtime = verifier_runtime + elapsed;
        }



        if SUM_CHECK {
            let now = Instant::now();
            let (claim, conj_witness) = hp::protocol_wrappers::prover::norm_1(&witness);
            let elapsed = now.elapsed();
            println_with_timestamp!("Time for prover::norm_1: {:.2?}", elapsed);
            prover_runtime = prover_runtime + elapsed;

            let now = Instant::now();
            hp::protocol_wrappers::verifier::verify_norm_1(&claim);
            let elapsed = now.elapsed();
            println_with_timestamp!("Time for verifier::verify_norm_1: {:.2?}", elapsed);
            verifier_runtime = verifier_runtime + elapsed;

            let mut inner_verifier_runtime = Instant::now().elapsed();

            let now = Instant::now();
            let (new_power_series, norm_3_output) = hp::protocol_wrappers::prover::norm_2(&statement, &witness, &conj_witness, |norm_2_output: &HPNorm2Output| {
                let now_inner = Instant::now();
                let challenge = hp::protocol_wrappers::verifier::verify_norm_2(&norm_2_output);
                let elapsed_inner = now_inner.elapsed();
                inner_verifier_runtime = inner_verifier_runtime + elapsed_inner;
                challenge
            });
            let elapsed = now.elapsed() - inner_verifier_runtime;
            println_with_timestamp!("Time for prover::norm_2: {:.2?}", elapsed);
            println_with_timestamp!("Time for verifier::verify_norm_2: {:.2?}", inner_verifier_runtime);
            prover_runtime = prover_runtime + elapsed;
            verifier_runtime = verifier_runtime + inner_verifier_runtime;

            let now = Instant::now();
            verifier_state =
                hp::protocol_wrappers::verifier::verify_norm_3(&norm_3_output, &verifier_state);
            let elapsed = now.elapsed();
            println_with_timestamp!("Time for verifier::verify_norm_3: {:.2?}", elapsed);
            verifier_runtime = verifier_runtime + elapsed;

            statement = new_power_series;
        } else {

            let now = Instant::now();
            let (new_witness, norm_1_output) = norm_1(&statement, &witness);
            witness = new_witness;
            let elapsed = now.elapsed();
            println_with_timestamp!("Time for norm_1: {:.2?}", elapsed);
            prover_runtime = prover_runtime + elapsed;

            let now = Instant::now();
            let (new_verifier_state, challenge, inverse_challenge) = norm_challenge(&norm_1_output, &verifier_state);
            verifier_state = new_verifier_state;
            let elapsed = now.elapsed();
            println_with_timestamp!("Time for norm_challenge: {:.2?}", elapsed);
            verifier_runtime = verifier_runtime + elapsed;
            //
            let now = Instant::now();
            let (new_power_series, norm_2_output) = norm_2(&statement, &witness, &challenge, &inverse_challenge);
            statement = new_power_series;
            let elapsed = now.elapsed();
            println_with_timestamp!("Time for norm_2: {:.2?}", elapsed);
            prover_runtime = prover_runtime + elapsed;

            let now = Instant::now();
            verifier_state = verify_norm_2(&norm_1_output, &norm_2_output, &verifier_state);
            let elapsed = now.elapsed();
            println_with_timestamp!("Time for verify_norm_2: {:.2?}", elapsed);
            verifier_runtime = verifier_runtime + elapsed;
        }


        let now = Instant::now();
        let (new_witness, split_output) = split(&mut statement, &witness);
        witness = new_witness;
        let elapsed = now.elapsed();
        println_with_timestamp!("Time for split: {:.2?}", elapsed);
        prover_runtime = prover_runtime + elapsed;


        let now = Instant::now();
        verifier_state = verifier_split(&statement, split_output, &verifier_state);
        let elapsed = now.elapsed();
        println_with_timestamp!("Time for split verifier: {:.2?}", elapsed);

        let now = Instant::now();
        let challenge = challenge_for_fold(&verifier_state);
        let elapsed = now.elapsed();
        println_with_timestamp!("Time for challenge fold: {:.2?}", elapsed);
        verifier_runtime = verifier_runtime + elapsed;

        let now = Instant::now();
        let new_witness = fold(&witness, &challenge);
        witness = new_witness;
        let elapsed = now.elapsed();
        println_with_timestamp!("Time for fold: {:.2?}", elapsed);
        prover_runtime = prover_runtime + elapsed;

        let now = Instant::now();
        verifier_state = verifier_fold(&verifier_state, &challenge);
        let elapsed = now.elapsed();
        println_with_timestamp!("Time for fold verifier: {:.2?}", elapsed);
        verifier_runtime = verifier_runtime + elapsed;
    }

    let now = Instant::now();
    assert_eq!(parallel_dot_series_matrix(
        &statement,
        &witness
    ), verifier_state.rhs);
    let elapsed = now.elapsed();
    println_with_timestamp!("Time for final assert_eq: {:.2?}", elapsed);
    verifier_runtime = verifier_runtime + elapsed;


    println_with_timestamp!("PRV: {:.2?}", prover_runtime);
    println_with_timestamp!("VER: {:.2?}", verifier_runtime);
}
