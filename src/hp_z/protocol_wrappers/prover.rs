use log::trace;
use crate::arithmetic::{compute_hp_power_series, ring_inner_product, transpose, PowerSeries};
use crate::hp_z::sum_check_z::{g_ts, sq_sum_check_batched_init, sq_sum_check_batched_par, sq_sum_check_single_init, verify_polynomials_and_generate_challenge};
use crate::hp_z::zq::{Zq, ZqElement};
use crate::prime_ring::ring::DPrimeRingElement;

pub struct HPNorm1Output {
    pub claims: Vec<ZqElement>
}

pub struct HPNorm2Output {
    pub p_coefs_all: Vec<Vec<DPrimeRingElement>>
}

pub struct HPNorm3Output {
    pub unsquared_claims: Vec<DPrimeRingElement>
}


pub fn norm_1(
    witness: &Vec<Vec<DPrimeRingElement>>
) -> HPNorm1Output {
    // let traced_witness = transpose(&witness.iter().map(|r| {
    //     r.iter().map(|e| Zq::from_u64(e.trace())).collect::<Vec<ZqElement>>()
    // }).collect());

    let claims = witness.iter().map(|c| {
        let c_conj = c.iter().map(DPrimeRingElement::conjugate).collect();
        let t = ring_inner_product(&c, &c_conj);
        Zq::from_u64(t.trace())
    }).collect();



    HPNorm1Output {
        claims
    }
}

pub fn norm_2(
    power_series: &Vec<PowerSeries>,
    witness: &Vec<Vec<DPrimeRingElement>>,
    mut verify_norm_2_wrapped: impl FnMut(&HPNorm2Output) -> DPrimeRingElement
) -> (Vec<PowerSeries>, HPNorm3Output) {

    let verify_norm_2 = |p_coefs_all: &Vec<Vec<DPrimeRingElement>>| {
        let output = HPNorm2Output {
            p_coefs_all: p_coefs_all.clone()
        };
        verify_norm_2_wrapped(&output)
    };
    let unsquared_claims = sq_sum_check_batched_par(
        &witness,
        verify_norm_2
    );


    unsafe {
        (vec![power_series.clone(), vec![compute_hp_power_series(&g_ts)]].concat(),
         HPNorm3Output {
             unsquared_claims
         })
    }
}

