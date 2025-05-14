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

type Val            = BabyBear;
type Perm           = Poseidon2BabyBear<16>;
type MyHash         = PaddingFreeSponge<Perm, 16, 8, 8>;
type MyCompress     = TruncatedPermutation<Perm, 2, 8, 16>;
type ValMmcs        = MerkleTreeMmcs<
                          <Val as Field>::Packing,
                          <Val as Field>::Packing,
                          MyHash,
                          MyCompress,
                          8>;
type Challenge      = BinomialExtensionField<Val, 4>;
type ChallengeMmcs  = ExtensionMmcs<Val, Challenge, ValMmcs>;
type Challenger     = DuplexChallenger<Val, Perm, 16, 8>;
type Dft            = Radix2DitParallel<Val>;
type Pcs            = TwoAdicFriPcs<Val, Dft, ValMmcs, ChallengeMmcs>;
type MyConfig       = StarkConfig<Pcs, Challenge, Challenger>;

const COLS: usize = 3;

pub struct Num2BitsAir;

impl<F: Field> BaseAir<F> for Num2BitsAir {
    fn width(&self) -> usize { COLS }
}

impl<AB: AirBuilderWithPublicValues> Air<AB> for Num2BitsAir {
    fn eval(&self, builder: &mut AB) {
        let main = builder.main();
        let (now, nxt) = (main.row_slice(0).unwrap(), main.row_slice(1).unwrap());
        let pis      = builder.public_values();
        let a_pub    = pis[0];

        builder.when_first_row().assert_eq(now[0], a_pub);

        /*
        def n2bits_recursive(n: int) -> None:
            """
            (LSB to MSB)
            """
            if n == 0:
                return
            print(n % 2, end='')
            n2bits_recursive(n // 2)
            19 9 1
            9 4 1
            4 2 0
            2 1 0
            1 0 1
            0 0 0
        */
        builder.when_transition().assert_eq(nxt[0] * AB::Expr::from_u64(2) + now[2], now[0]);
        builder.when_transition().assert_eq(nxt[0], now[1]);

        // Termination state: n == 0, so all three columns are zero.
        builder.when_last_row().assert_zero(now[0]);
        builder.when_last_row().assert_zero(now[1]);
        builder.when_last_row().assert_zero(now[2]);

    }
}

fn gen_trace<F: PrimeField64>(n: u64) -> RowMajorMatrix<F> {
    // Build the core rows following the recursion n → (n/2, n%2).
    let mut n_local = n;
    let mut rows: Vec<[F; COLS]> = Vec::new();
    while n_local > 0 {
        rows.push([
            F::from_u64(n_local),
            F::from_u64(n_local / 2),
            F::from_u64(n_local % 2),
        ]);
        n_local /= 2;
    }

    // Add the terminating state n == 0.
    rows.push([F::from_u64(0); COLS]);

    // Pad with additional zero‑rows so the total height is a power of two.
    let target_len = rows.len().next_power_of_two();
    rows.resize(target_len, [F::from_u64(0); COLS]);

    // Flatten into a single vector for RowMajorMatrix.
    RowMajorMatrix::new(rows.into_iter().flatten().collect(), COLS)
}

fn num2bits(n: u64) {
    // --- witness (= trace + pis) -------------------------------------------------
    let trace = gen_trace::<Val>(n);
    // println!("Trace (first row shown): {:?}", &trace.values[..COLS]);
    println!("▶ Trace (all rows):");
    for row_idx in 0..trace.height() {
        let row = trace.row(row_idx).unwrap();
        let values: Vec<u64> = row.into_iter()
                                  .map(|x| x.as_canonical_u64()) // BabyBear → u64
                                  .collect();
        println!("Row {:>3}: {:?}", row_idx, values);
    }

    let pis = vec![Val::from_u64(n)];
    println!("▶ Public inputs          : {:?}", pis.iter().map(|x| x.as_canonical_u64()).collect::<Vec<_>>());

    // --- STARK config --------------------------------------------------------------
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

    // --- Proof & Verify -------------------------------------------------------------
    let proof = prove(&cfg, &Num2BitsAir, trace, &pis);

    verify(&cfg, &Num2BitsAir, &proof, &pis).expect("verify failed");
    println!("Verification passed!");
}

fn main() {
    for i in 0..255 {
        num2bits(i);
    }
}
