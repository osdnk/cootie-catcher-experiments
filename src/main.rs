#![feature(adt_const_params)]
#![feature(const_trait_impl)]


use crate::helpers::println_with_timestamp;
use crate::prime_ring::r#static::{DECOMP_PARTS, MODULE_SIZE, MOD_Q, REP, WIT_DIM, PHI, SUM_CHECK};
use crate::prime_ring::ring::PrimeRing;
use crate::protocol::protocol;

mod prime_ring;

mod ntt;
mod cpp;
mod arithmetic;
mod subroutines;
//
mod helpers;
mod hp;
// mod hp_z;
mod protocol;

fn main() {
    cfg_if::cfg_if! {
        if #[cfg(feature = "A_OLD")] {
                println_with_timestamp!("A_OLD");
        } else if #[cfg(feature = "A_NEW")] {
                println_with_timestamp!("A_NEW");
        } else if #[cfg(feature = "B_OLD")] {
                println_with_timestamp!("B_OLD");
        } else if #[cfg(feature = "B_NEW")] {
                println_with_timestamp!("B_NEW");
        } else if #[cfg(feature = "C_OLD")] {
                println_with_timestamp!("C_OLD");
        } else if #[cfg(feature = "C_NEW")] {
                println_with_timestamp!("C_NEW");
        } else if #[cfg(feature = "D_OLD")] {
                println_with_timestamp!("D_OLD");
        } else if #[cfg(feature = "D_NEW")] {
                println_with_timestamp!("D_NEW");
        } else if #[cfg(feature = "E_OLD")] {
                println_with_timestamp!("E_OLD");
        } else if #[cfg(feature = "E_NEW")] {
                println_with_timestamp!("E_NEW");
        } else if #[cfg(feature = "F_OLD")] {
                println_with_timestamp!("F_OLD");
        } else if #[cfg(feature = "F_NEW")] {
                println_with_timestamp!("F_NEW");
        } else if #[cfg(feature = "F_OLD_2")] {
                println_with_timestamp!("F_OLD_2");
        } else if #[cfg(feature = "F_NEW_2")] {
                println_with_timestamp!("F_NEW_2");
        } else if #[cfg(feature = "F_OLD_3")] {
                println_with_timestamp!("F_OLD_3");
        } else if #[cfg(feature = "F_NEW_3")] {
                println_with_timestamp!("F_NEW_3");
        } else if #[cfg(feature = "F_OLD_4")] {
                println_with_timestamp!("F_OLD_4");
        } else if #[cfg(feature = "F_NEW_4")] {
                println_with_timestamp!("F_NEW_4");
        } else {
                println_with_timestamp!("default");
        }
    }

    cfg_if::cfg_if! {
        if #[cfg(feature = "use-hardware")] {
                println_with_timestamp!("hardware acc full");
        } else {
                cfg_if::cfg_if! {
                    if #[cfg(feature = "partial-hardware")] {
                            println_with_timestamp!("partial hardware acc");
                    } else {
                            println_with_timestamp!("no hardware acc");
                    }
                }
        }
    }



    println_with_timestamp!(
        "PARAMS: NEW: {:?}, MODULE: {:?}, WIT_DIM: {:?}, REP: {:?}, Q: {:?}, PHI: {:?}, DECOMP_PARTS: {:?}",
        SUM_CHECK, MODULE_SIZE, WIT_DIM, REP, MOD_Q, PHI, DECOMP_PARTS);
    let a = PrimeRing::random();
    let b = a.inverse();
    assert_eq!(a * b, PrimeRing::constant(1));
    println_with_timestamp!("OK sage");

    // println!("hi");
    protocol();
    // print_stats();
}
