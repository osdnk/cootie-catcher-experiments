use lazy_static::lazy_static;
use once_cell::sync::Lazy;
use rug::Integer;
use crate::prime_ring::ring::{DPrimeRingElement, PrimeRing, PrimeRingElement};
use crate::ntt::root_of_unity::mod_inverse;

// pub const PHI: usize = 6;

pub const n_fhe: usize = 723;
// pub const PHI: usize = 6;
// pub const n_fhe: usize = 5;
pub const LOG_CONDUCTOR: usize = 10;
pub const TWO_PHI_MINUS_ONE: usize = 2 * PHI - 1;
pub const TEST_PHI: usize = 6;
pub const TEST_TWO_PHI_MINUS_ONE: usize = 11;
pub static NOF_THREADS_UPPED_BOUND: usize = 120;

//pub const MOD_Q: u64 = 4611686078556930049;
//pub const MOD_Q: u64 = 1152921493869428737;
//pub const MOD_Q: u64 = 2305842979148922881;
//pub const MOD_Q: u64 = 2305843095113039873;
pub const MOD_Q: u64 = 4546383823830515713;
pub static ONE: Lazy<DPrimeRingElement> = Lazy::new(|| {
    PrimeRing::constant(1)
});

pub static TWO: Lazy<DPrimeRingElement> = Lazy::new(|| {
    PrimeRing::constant(2)
});



pub const BASIS: Lazy<Vec<DPrimeRingElement>> = Lazy::new(|| {
    (0..PHI).into_iter().map(|i| {
        let mut coeffs_t = [0u64; CONDUCTOR];
        coeffs_t[i] = 1;
        let t = PrimeRing::new_from_larger_vec(&coeffs_t.try_into().unwrap());
        t
    }).collect()
});

pub static  LOG_Q: usize = 148;

pub const MAX_THREADS: usize = 120;

// pub static WIT_DIM: usize = 16384;


cfg_if::cfg_if! {
    if #[cfg(feature = "A_OLD")] {
            pub const PHI: usize = 16;
            pub static MODULE_SIZE: usize = 65;
            pub static WIT_DIM: usize = 8388608;
            pub static REP: usize = 22;
            pub static INIT_REP: usize = 22;
            pub static NOF_ROUNDS: i32 = 13;
            pub static DECOMP_PARTS: [usize; 13] = [1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2];
            pub static SUM_CHECK: bool = false;
    } else if #[cfg(feature = "A_NEW")] {
            pub const PHI: usize = 16;
            pub static MODULE_SIZE: usize = 59;
            pub static WIT_DIM: usize = 8388608;
            pub static REP: usize = 20;
            pub static INIT_REP: usize = 20;
            pub static NOF_ROUNDS: i32 = 13;
            pub static DECOMP_PARTS: [usize; 13] = [1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2];
            pub static SUM_CHECK: bool = true;

    } else if #[cfg(feature = "B_OLD")] {
            pub const PHI: usize = 30;
            pub static MODULE_SIZE: usize = 39;
            pub static WIT_DIM: usize = 4194304;
            pub static REP: usize = 18;
            pub static INIT_REP: usize = 18;
            pub static NOF_ROUNDS: i32 = 13;
            pub static DECOMP_PARTS: [usize; 13] = [1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2];
            pub static SUM_CHECK: bool = false;

    } else if #[cfg(feature = "B_NEW")] {
            pub const PHI: usize = 30;
            pub static MODULE_SIZE: usize = 34;
            pub static WIT_DIM: usize = 4194304;
            pub static REP: usize = 18;
            pub static INIT_REP: usize = 18;
            pub static NOF_ROUNDS: i32 = 13;
            pub static DECOMP_PARTS: [usize; 13] = [1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2];
            pub static SUM_CHECK: bool = true;

    } else if #[cfg(feature = "C_OLD")] {
            pub const PHI: usize = 112;
            pub static MODULE_SIZE: usize = 14;
            pub static WIT_DIM: usize = 2097152;
            pub static REP: usize = 14;
            pub static INIT_REP: usize = 14;
            pub static NOF_ROUNDS: i32 = 13;
            pub static DECOMP_PARTS: [usize; 13] = [1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2];
            pub static SUM_CHECK: bool = false;

    } else if #[cfg(feature = "C_NEW")] {
            pub const PHI: usize = 112;
            pub static MODULE_SIZE: usize = 12;
            pub static WIT_DIM: usize = 2097152;
            pub static REP: usize = 14;
            pub static INIT_REP: usize = 14;
            pub static NOF_ROUNDS: i32 = 14;
            pub static DECOMP_PARTS: [usize; 14] = [1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2];
            pub static SUM_CHECK: bool = true;

    } else if #[cfg(feature = "D_OLD")] { // DD
            pub const PHI: usize = 250;
            pub static MODULE_SIZE: usize = 7;
            pub static WIT_DIM: usize = 1048576;
            pub static REP: usize = 12;
            pub static INIT_REP: usize = 12;
            pub static NOF_ROUNDS: i32 = 13;
            pub static DECOMP_PARTS: [usize; 13] = [1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2];
            pub static SUM_CHECK: bool = false;

    } else if #[cfg(feature = "D_NEW")] {
            pub const PHI: usize = 250;
            pub static MODULE_SIZE: usize = 7;
            pub static WIT_DIM: usize = 1048576;
            pub static REP: usize = 12;
            pub static INIT_REP: usize = 12;
            pub static NOF_ROUNDS: i32 = 13;
            pub static DECOMP_PARTS: [usize; 13] = [1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2];
            pub static SUM_CHECK: bool = true;

    } else if #[cfg(feature = "E_OLD")] {
            pub const PHI: usize = 508;
            pub static MODULE_SIZE: usize = 4;
            pub static WIT_DIM: usize = 524288;
            pub static REP: usize = 10;
            pub static INIT_REP: usize = 10;
            pub static NOF_ROUNDS: i32 = 12;
            pub static DECOMP_PARTS: [usize; 12] = [1, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2];
            pub static SUM_CHECK: bool = false;

    } else if #[cfg(feature = "E_NEW")] {
            pub const PHI: usize = 508;
            pub static MODULE_SIZE: usize = 4;
            pub static WIT_DIM: usize = 524288;
            pub static REP: usize = 10;
            pub static INIT_REP: usize = 10;
            pub static NOF_ROUNDS: i32 = 13;
            pub static DECOMP_PARTS: [usize; 13] = [1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2];
            pub static SUM_CHECK: bool = true;

    } else if #[cfg(feature = "F_OLD")] {
            pub const PHI: usize = 1020;
            pub static MODULE_SIZE: usize = 3;
            pub static WIT_DIM: usize = 262144;
            pub static REP: usize = 9;
            pub static INIT_REP: usize = 9;
            pub static NOF_ROUNDS: i32 = 12;
            pub static DECOMP_PARTS: [usize; 12] = [1, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2];
            pub static SUM_CHECK: bool = false;

    } else if #[cfg(feature = "F_NEW")] {
            pub const PHI: usize = 1020;
            pub static MODULE_SIZE: usize = 2;
            pub static WIT_DIM: usize = 262144;
            pub static REP: usize = 9;
            pub static INIT_REP: usize = 9;
            pub static NOF_ROUNDS: i32 = 12;
            pub static DECOMP_PARTS: [usize; 12] = [1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2];
            pub static SUM_CHECK: bool = true;

    } else if #[cfg(feature = "F_OLD_2")] {
            pub const PHI: usize = 1020;
            pub static MODULE_SIZE: usize = 2;
            pub static WIT_DIM: usize = 131072;
            pub static REP: usize = 9;
            pub static INIT_REP: usize = 9;
            pub static NOF_ROUNDS: i32 = 11;
            pub static DECOMP_PARTS: [usize; 11] = [1, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2];
            pub static SUM_CHECK: bool = false;

    } else if #[cfg(feature = "F_NEW_2")] {
            pub const PHI: usize = 1020;
            pub static MODULE_SIZE: usize = 2;
            pub static WIT_DIM: usize = 131072;
            pub static REP: usize = 9;
            pub static INIT_REP: usize = 9;
            pub static NOF_ROUNDS: i32 = 11;
            pub static DECOMP_PARTS: [usize; 11] = [1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2];
            pub static SUM_CHECK: bool = true;

    } else if #[cfg(feature = "F_OLD_3")] {
            pub const PHI: usize = 1020;
            pub static MODULE_SIZE: usize = 2;
            pub static WIT_DIM: usize = 65536;
            pub static REP: usize = 9;
            pub static INIT_REP: usize = 9;
            pub static NOF_ROUNDS: i32 = 10;
            pub static DECOMP_PARTS: [usize; 10] = [1, 3, 3, 3, 3, 2, 2, 2, 2, 2];
            pub static SUM_CHECK: bool = false;

    } else if #[cfg(feature = "F_NEW_3")] {
            pub const PHI: usize = 1020;
            pub static MODULE_SIZE: usize = 2;
            pub static WIT_DIM: usize = 65536;
            pub static REP: usize = 9;
            pub static INIT_REP: usize = 9;
            pub static NOF_ROUNDS: i32 = 10;
            pub static DECOMP_PARTS: [usize; 10] = [1, 2, 2, 2, 2, 2, 2, 2, 2, 2];
            pub static SUM_CHECK: bool = true;

    } else if #[cfg(feature = "F_OLD_4")] {
            pub const PHI: usize = 1020;
            pub static MODULE_SIZE: usize = 2;
            pub static WIT_DIM: usize = 32768;
            pub static REP: usize = 9;
            pub static INIT_REP: usize = 9;
            pub static NOF_ROUNDS: i32 = 9;
            pub static DECOMP_PARTS: [usize; 9] = [1, 3, 3, 3, 2, 2, 2, 2, 2];
            pub static SUM_CHECK: bool = false;

    } else if #[cfg(feature = "F_NEW_4")] {
            pub const PHI: usize = 1020;
            pub static MODULE_SIZE: usize = 2;
            pub static WIT_DIM: usize = 32768;
            pub static REP: usize = 9;
            pub static INIT_REP: usize = 9;
            pub static NOF_ROUNDS: i32 = 9;
            pub static DECOMP_PARTS: [usize; 9] = [1, 2, 2, 2, 2, 2, 2, 2, 2];
            pub static SUM_CHECK: bool = true;

    }else if #[cfg(feature = "F_OLD_0")] {
            pub const PHI: usize = 1020;
            pub static MODULE_SIZE: usize = 3;
            pub static WIT_DIM: usize = 524288;
            pub static REP: usize = 9;
            pub static INIT_REP: usize = 9;
            pub static NOF_ROUNDS: i32 = 13;
            pub static DECOMP_PARTS: [usize; 13] = [1, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2];

            pub static SUM_CHECK: bool = false;

    } else if #[cfg(feature = "F_NEW_0")] {
            pub const PHI: usize = 1020;
            pub static MODULE_SIZE: usize = 2;
            pub static WIT_DIM: usize = 524288;
            pub static REP: usize = 9;
            pub static INIT_REP: usize = 9;
            pub static NOF_ROUNDS: i32 = 13;
            pub static DECOMP_PARTS: [usize; 13] = [1, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2];

            pub static SUM_CHECK: bool = true;


    }else if #[cfg(feature = "F_OLD_00")] {
            pub const PHI: usize = 1020;
            pub static MODULE_SIZE: usize = 3;
            pub static WIT_DIM: usize = 1048576;
            pub static REP: usize = 9;
            pub static INIT_REP: usize = 9;
            pub static NOF_ROUNDS: i32 = 13;
            pub static DECOMP_PARTS: [usize; 13] = [1, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2];
            pub static SUM_CHECK: bool = false;

    } else if #[cfg(feature = "F_NEW_00")] {
            pub const PHI: usize = 1020;
            pub static MODULE_SIZE: usize = 2;
            pub static WIT_DIM: usize = 1048576;
            pub static REP: usize = 9;
            pub static INIT_REP: usize = 9;
            pub static NOF_ROUNDS: i32 = 13;
            pub static DECOMP_PARTS: [usize; 13] = [1, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2];
            pub static SUM_CHECK: bool = true;

    } else {
            pub const PHI: usize = 112;
            pub static MODULE_SIZE: usize = 1;
            pub static WIT_DIM: usize = 2097152 / 8;
            pub static REP: usize = 1;
            pub static INIT_REP: usize = 14;
            pub static NOF_ROUNDS: i32 = 13;
            pub static DECOMP_PARTS: [usize; 13] = [1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2];
            pub static SUM_CHECK: bool = false;
    }
}

pub const CONDUCTOR: usize = PHI + 1;
