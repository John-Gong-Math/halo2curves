use crate::arithmetic::mul_896;
use crate::pluto_eris::Fp;
use crate::pluto_eris::Fq;

const GAMMA1: [u32; 14] = [0x6a19d56d,
0xdf52222a,
0x8e5c9500,
0x2aae415b,
0x4a0aaaab,
0xaaa95555,
0xaaaaaaaa,
0x00000002,
0x00000000,
0x00000000,
0x00000000,
0x00000000,
0x00000000,
0x00000000];

const GAMMA2: [u32; 14] = [0x38e4e395,
0xe38e224e,
0xe38e38e0,
0x00038e38,
0x00000000,
0x00000000,
0x00000000,
0x00000000,
0x00000000,
0x00000000,
0x00000000,
0x00000000,
0x00000000,
0x00000000];

const B1: [u32; 14] = [0x00000001,
0x3ffffde2,
0xa001ae51,
0x05ff0065,
0x01968000,
0x00003000,
0x60000000,
0x00000000,
0x00000000,
0x00000000,
0x00000000,
0x00000000,
0x00000000,
0x00000000];

const B2: [u32; 14] = [0xffffffff,
0x2000010e,
0x00000000,
0x00008000,
0x00000000,
0x00000000,
0x00000000,
0x00000000,
0x00000000,
0x00000000,
0x00000000,
0x00000000,
0x00000000,
0x00000000];