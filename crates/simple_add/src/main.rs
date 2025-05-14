// BaseAir: inform the framework “I have 3 columns”
// Air: write the constraints using AirBuilder
// AirBuilder: toolbox for grabbing variables, declaring equations, and adding conditions
// AirBuilderWithPublicValues (optional): expose public columns directly to the verifier
use p3_air::{Air, AirBuilder, AirBuilderWithPublicValues, BaseAir};
use p3_matrix::Matrix;
use p3_matrix::dense::RowMajorMatrix;

use p3_baby_bear::{BabyBear, Poseidon2BabyBear};
use p3_field::PrimeCharacteristicRing;

use p3_challenger::DuplexChallenger;
use p3_commit::ExtensionMmcs;
use p3_dft::Radix2DitParallel;
use p3_field::extension::BinomialExtensionField;
use p3_field::{Field, PrimeField64};
use p3_fri::{create_test_fri_config, TwoAdicFriPcs};

use p3_merkle_tree::MerkleTreeMmcs;
use p3_symmetric::{PaddingFreeSponge, TruncatedPermutation};
use p3_uni_stark::{prove, verify, StarkConfig};
use rand::rng;

use bincode;

// -----------------------------------------------------------------------------
// 0. type alise
// -----------------------------------------------------------------------------
// Base field
type Val            = BabyBear;   
// Oracle
type Perm           = Poseidon2BabyBear<16>;
// For reduce the size of Merkle Tree
type MyHash         = PaddingFreeSponge<Perm, 16, 8, 8>;
// TODO: compression function how it works and why we need it
type MyCompress     = TruncatedPermutation<Perm, 2, 8, 16>;
// Merkle Tree on Base field with compression function
type ValMmcs        = MerkleTreeMmcs<
                          <Val as Field>::Packing,
                          <Val as Field>::Packing,
                          MyHash,
                          MyCompress,
                          8>;
// Extension field log_blowup is the extension parameter
type Challenge      = BinomialExtensionField<Val, 4>;
// Merkle Tree on Extension field with compression function
type ChallengeMmcs  = ExtensionMmcs<Val, Challenge, ValMmcs>;
// Fiat-shamir
type Challenger     = DuplexChallenger<Val, Perm, 16, 8>;
// Choose one from vairous FFT implemention
type Dft            = Radix2DitParallel<Val>;
// Pcs
type Pcs            = TwoAdicFriPcs<Val, Dft, ValMmcs, ChallengeMmcs>;
// Stark Config 
type MyConfig       = StarkConfig<Pcs, Challenge, Challenger>;

// -----------------------------------------------------------------------------
// 1. AIR
// -----------------------------------------------------------------------------
const COLS: usize = 3;

pub struct AddAir;

impl<F: Field> BaseAir<F> for AddAir {
    fn width(&self) -> usize { COLS }
}

impl<AB: AirBuilderWithPublicValues> Air<AB> for AddAir {
    fn eval(&self, builder: &mut AB) {
        let main = builder.main();
        let now= main.row_slice(0).unwrap();
        let pis      = builder.public_values();
        let c_pub    = pis[0];

        builder.when_first_row().assert_eq(now[2], c_pub);
        builder.when_first_row().assert_eq(now[0] + now[1], now[2]);

    }
}

// -----------------------------------------------------------------------------
// 2. generate trace by the table and
// -----------------------------------------------------------------------------
fn gen_trace<F: PrimeField64>(a: u64, b: u64, c: u64) -> RowMajorMatrix<F> {
    let mut v = Vec::with_capacity(COLS);
    v.push(F::from_u64(a));
    v.push(F::from_u64(b));
    v.push(F::from_u64(c));
    RowMajorMatrix::new(v, COLS)
}

// -----------------------------------------------------------------------------
// 3. config → prove → verify, print witness / proof
// -----------------------------------------------------------------------------
fn main() {
    let (a, b) = (13_u64, 55_u64);
    // let c = 34_u64;
    let c = a + b;

    // --- witness (= trace + pis) -------------------------------------------------
    let trace = gen_trace::<Val>(a, b, c);
    println!("Trace (first row shown): {:?}", &trace.values[..COLS]);

    let pis = vec![Val::from_u64(c)];
    println!("Public inputs          : {:?}", pis.iter().map(|x| x.as_canonical_u64()).collect::<Vec<_>>());

    // --- STARK config --------------------------------------------------------------
    let perm           = Perm::new_from_rng_128(&mut rng());
    let hash           = MyHash::new(perm.clone());
    let compress       = MyCompress::new(perm.clone());
    let val_mmcs       = ValMmcs::new(hash, compress);
    let challenge_mmcs = ChallengeMmcs::new(val_mmcs.clone());
    let dft            = Dft::default();
    // log_final_poly_len = 0 the degree of final fold polynomial
    let fri_cfg        = create_test_fri_config(challenge_mmcs, 0);
    let pcs            = Pcs::new(dft, val_mmcs, fri_cfg);
    let challenger     = Challenger::new(perm);
    let cfg            = MyConfig::new(pcs, challenger);

    // --- Proof & Verify -------------------------------------------------------------
    let proof = prove(&cfg, &AddAir, trace, &pis);
    let proof_bytes = bincode::serialize(&proof).expect("Failed to serialize proof");
    println!("Proof length (bytes)   : {}", proof_bytes.len());

    verify(&cfg, &AddAir, &proof, &pis).expect("verify failed");
    println!("Verification passed!");
}