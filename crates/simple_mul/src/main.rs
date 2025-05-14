
use p3_air::{Air, AirBuilder, AirBuilderWithPublicValues, BaseAir};
use p3_matrix::Matrix;
use p3_matrix::dense::RowMajorMatrix;

use p3_baby_bear::{BabyBear, Poseidon2BabyBear};
use p3_field::{Field, PrimeField64, PrimeCharacteristicRing};

use p3_challenger::DuplexChallenger;
use p3_commit::ExtensionMmcs;
use p3_dft::Radix2DitParallel;

use p3_field::extension::BinomialExtensionField;
use p3_fri::{create_test_fri_config, TwoAdicFriPcs};

use p3_merkle_tree::MerkleTreeMmcs;
use p3_symmetric::{PaddingFreeSponge, TruncatedPermutation};
use p3_uni_stark::{prove, verify, StarkConfig};
use rand::rng;

type Val = BabyBear;
type Perm = Poseidon2BabyBear<16>;
type MyHash = PaddingFreeSponge<Perm, 16, 8, 8>;
type MyCompress     = TruncatedPermutation<Perm, 2, 8, 16>;  
type ValMmcs        = MerkleTreeMmcs<
                          <Val as Field>::Packing,
                          <Val as Field>::Packing,
                          MyHash,
                          MyCompress,
                          8>;
type Challenge = BinomialExtensionField<Val, 4>;
type ChallengeMmcs = ExtensionMmcs<Val, Challenge, ValMmcs>;
type Challenger = DuplexChallenger<Val, Perm, 16, 8>;
type Dft = Radix2DitParallel<Val>;
type Pcs = TwoAdicFriPcs<Val, Dft, ValMmcs, ChallengeMmcs>;
type MyConfig = StarkConfig<Pcs, Challenge, Challenger>;

const COLS: usize = 3;

pub struct MulAir;

impl<F: Field> BaseAir<F> for MulAir {
    fn width(&self) -> usize { COLS }
}

impl<AB: AirBuilderWithPublicValues> Air<AB> for MulAir {
    fn eval(&self, builder: &mut AB) {
        let main = builder.main();
        let now = main.row_slice(0).unwrap();
        let pis = builder.public_values();
        let c_pub = pis[0];

        builder.when_first_row().assert_eq(now[2], c_pub);
        builder.when_first_row().assert_eq(now[0] * now[1], now[2]);
    }
}

fn gen_trace<F: PrimeField64>(a: u64, b: u64, c: u64) -> RowMajorMatrix<F> {
    let mut v = Vec::with_capacity(COLS);
    v.push(F::from_u64(a));
    v.push(F::from_u64(b));
    v.push(F::from_u64(c));
    RowMajorMatrix::new(v, COLS)
}

fn main() {
    let (a, b) = (13_u64, 55_u64);
    let c = 715_u64;

    let trace = gen_trace::<Val>(a, b, c);
    println!("▶ Trace (first row shown): {:?}", &trace.values[..COLS]);

    let pis = vec![Val::from_u64(c)];
    println!("▶ Public inputs          : {:?}", pis.iter().map(|x| x.as_canonical_u64()).collect::<Vec<_>>());
    
    let perm           = Perm::new_from_rng_128(&mut rng());
    let hash           = MyHash::new(perm.clone());
    let compress       = MyCompress::new(perm.clone());
    let val_mmcs       = ValMmcs::new(hash, compress);
    let challenge_mmcs = ChallengeMmcs::new(val_mmcs.clone());
    let dft            = Dft::default();
    let fri_cfg        = create_test_fri_config(challenge_mmcs, 0);
    let pcs            = Pcs::new(dft, val_mmcs, fri_cfg);
    let challenger     = Challenger::new(perm);
    let cfg            = MyConfig::new(pcs, challenger);

    let proof = prove(&cfg, &MulAir, trace, &pis);
    verify(&cfg, &MulAir, &proof, &pis).expect("verify failed");
    println!("✅ Verification passed!");
}
