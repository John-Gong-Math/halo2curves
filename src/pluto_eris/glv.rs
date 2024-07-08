use super::G1;
use crate::arithmetic::mul_896;
use crate::arithmetic::CurveEndo;
use crate::arithmetic::EndoParameters;
use crate::ff::WithSmallOrderMulGroup;
use crate::group::{Curve, Group};
use crate::pluto_eris::Fq;
use crate::CurveExt;
use ethnum::U256;
use ff::{Field, PrimeField};
use rand_core::OsRng;
use subtle::{Choice, ConditionallySelectable, ConstantTimeEq, CtOption};

const GAMMA1: [u64; 7] = [
    0xdf52222a6a19d56d,
    0x2aae415b8e5c9500,
    0xaaa955554a0aaaab,
    0x00000002aaaaaaaa,
    0x0000000000000000,
    0x0000000000000000,
    0x0000000000000000,
];

const GAMMA2: [u64; 7] = [
    0xe38e224e38e4e395,
    0x00038e38e38e38e0,
    0x0000000000000000,
    0x0000000000000000,
    0x0000000000000000,
    0x0000000000000000,
    0x0000000000000000,
];

const B1: [u64; 7] = [
    0x3ffffde200000001,
    0x05ff0065a001ae51,
    0x0000300001968000,
    0x0000000060000000,
    0x0000000000000000,
    0x0000000000000000,
    0x0000000000000000,
];

const B2: [u64; 7] = [
    0x2000010effffffff,
    0x0000800000000000,
    0x0000000000000000,
    0x0000000000000000,
    0x0000000000000000,
    0x0000000000000000,
    0x0000000000000000,
];

fn decompose_scalar(k: &Fq) -> (U256, bool, U256, bool) {
    let to_limbs = |e: &Fq| {
        let repr = e.to_repr();
        let repr = repr.as_ref();
        let tmp0 = u64::from_le_bytes(repr[0..8].try_into().unwrap());
        let tmp1 = u64::from_le_bytes(repr[8..16].try_into().unwrap());
        let tmp2 = u64::from_le_bytes(repr[16..24].try_into().unwrap());
        let tmp3 = u64::from_le_bytes(repr[24..32].try_into().unwrap());
        let tmp4 = u64::from_le_bytes(repr[32..40].try_into().unwrap());
        let tmp5 = u64::from_le_bytes(repr[40..48].try_into().unwrap());
        let tmp6 = u64::from_le_bytes(repr[48..56].try_into().unwrap());
        [tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6]
    };

    let get_low_high_224 = |e: &Fq| {
        let repr = e.to_repr();
        let repr = repr.as_ref();
        let mut low_224: [u8; 32] = [0; 32];
        let (one, two) = low_224.split_at_mut(28);
        let padding: [u8; 4] = [0; 4];
        one.copy_from_slice(&repr[0..28]);
        two.copy_from_slice(&padding);
        let low_224 = U256::from_le_bytes(low_224);

        let mut high_224: [u8; 32] = [0; 32];
        let (one, two) = high_224.split_at_mut(28);

        one.copy_from_slice(&repr[28..56]);
        two.copy_from_slice(&padding);
        let high_224 = U256::from_le_bytes(high_224);

        (low_224, high_224)
    };

    let is_neg = |e: &Fq| {
        let (_, high_224) = get_low_high_224(e);
        if high_224 != 0 {
            return true;
        }

        false
    };

    let input = to_limbs(k);
    let c1 = mul_896(GAMMA2, input);
    let c2 = mul_896(GAMMA1, input);
    let c1 = [c1[7], c1[8], c1[9], c1[10], c1[11], c1[12], c1[13]];
    let c2 = [c2[7], c2[8], c2[9], c2[10], c2[11], c2[12], c2[13]];
    let q1 = mul_896(c1, B1);
    let q2 = mul_896(c2, B2);
    let q1 = Fq::from_raw([q1[0], q1[1], q1[2], q1[3], q1[4], q1[5], q1[6]]);
    let q2 = Fq::from_raw([q2[0], q2[1], q2[2], q2[3], q2[4], q2[5], q2[6]]);
    let k2 = q2 - q1;
    let k1 = k + k2 * Fq::ZETA;
    let k1_neg = is_neg(&k1);
    let k2_neg = is_neg(&k2);
    let k1 = if k1_neg { -k1 } else { k1 };
    let k2 = if k2_neg { -k2 } else { k2 };

    (
        get_low_high_224(&k1).0,
        k1_neg,
        get_low_high_224(&k2).0,
        k2_neg,
    )
}

#[cfg(test)]
impl G1 {
    fn mul_U256(&self, scalar: U256) -> Self {
        let mut acc = Self::identity();
        for bit in scalar
            .to_be_bytes()
            // .to_repr()
            // .as_ref()
            .iter()
            // .rev()
            .flat_map(|byte| (0..8).rev().map(move |i| Choice::from((byte >> i) & 1u8)))
        {
            acc = acc.double();
            acc = Self::conditional_select(&acc, &(acc + self), bit);
        }

        acc
    }

    fn glv(&self, s: &Fq) -> G1 {
        let (k1, k1_neg, k2, k2_neg) = decompose_scalar(s);

        // should msm this
        let p1 = self.mul_U256(k1);
        let p2 = self.mul_U256(k2).endo();

        if k1_neg & k2_neg {
            p2 - p1
        } else if k1_neg {
            -(p1 + p2)
        } else if k2_neg {
            p1 + p2
        } else {
            p1 - p2
        }
    }
}

#[test]
fn test_glv() {
    for _ in 1..1000 {
        let s1 = <G1 as CurveExt>::ScalarExt::random(OsRng);

        let (k1, k1_neg, k2, k2_neg) = decompose_scalar(&s1);
        eprintln!(
            "Decomposition:
            K1: {k1},
            K1_neg: {k1_neg},
            K2: {k2},
            K2_neg: {k2_neg}"
        );

        let p = G1::random(OsRng);

        let r1 = p * s1;
        let r2 = p.glv(&s1);

        assert_eq!(r1.to_affine(), r2.to_affine())
    }
}
