mod common;

use cheq::Atom;
use common::{TestCase, run_group_test};

fn make_diatomic(z1: u8, z2: u8, dist: f64) -> Vec<Atom> {
    vec![
        Atom {
            atomic_number: z1,
            position: [0.0, 0.0, 0.0],
        },
        Atom {
            atomic_number: z2,
            position: [dist, 0.0, 0.0],
        },
    ]
}

#[test]
fn test_alkali_halides_group() {
    let cases = vec![
        TestCase {
            name: "LiF",
            atoms: make_diatomic(3, 9, 1.564),
            expected: vec![(0, 0.791)],
        },
        TestCase {
            name: "LiCl",
            atoms: make_diatomic(3, 17, 2.021),
            expected: vec![(0, 0.939)],
        },
        TestCase {
            name: "LiBr",
            atoms: make_diatomic(3, 35, 2.170),
            expected: vec![(0, 0.902)],
        },
        TestCase {
            name: "LiI",
            atoms: make_diatomic(3, 53, 2.392),
            expected: vec![(0, 0.841)],
        },
        TestCase {
            name: "NaF",
            atoms: make_diatomic(11, 9, 1.926),
            expected: vec![(0, 0.665)],
        },
        TestCase {
            name: "NaCl",
            atoms: make_diatomic(11, 17, 2.361),
            expected: vec![(0, 0.766)],
        },
        TestCase {
            name: "NaBr",
            atoms: make_diatomic(11, 35, 2.502),
            expected: vec![(0, 0.745)],
        },
        TestCase {
            name: "NaI",
            atoms: make_diatomic(11, 53, 2.711),
            expected: vec![(0, 0.709)],
        },
        TestCase {
            name: "KF",
            atoms: make_diatomic(19, 9, 2.171),
            expected: vec![(0, 0.662)],
        },
        TestCase {
            name: "KCl",
            atoms: make_diatomic(19, 17, 2.667),
            expected: vec![(0, 0.775)],
        },
        TestCase {
            name: "KBr",
            atoms: make_diatomic(19, 35, 2.821),
            expected: vec![(0, 0.768)],
        },
        TestCase {
            name: "KI",
            atoms: make_diatomic(19, 53, 3.048),
            expected: vec![(0, 0.754)],
        },
        TestCase {
            name: "RbF",
            atoms: make_diatomic(37, 9, 2.270),
            expected: vec![(0, 0.653)],
        },
        TestCase {
            name: "RbCl",
            atoms: make_diatomic(37, 17, 2.787),
            expected: vec![(0, 0.763)],
        },
        TestCase {
            name: "RbBr",
            atoms: make_diatomic(37, 35, 2.945),
            expected: vec![(0, 0.757)],
        },
        TestCase {
            name: "RbI",
            atoms: make_diatomic(37, 53, 3.177),
            expected: vec![(0, 0.747)],
        },
        TestCase {
            name: "CsF",
            atoms: make_diatomic(55, 9, 2.345),
            expected: vec![(0, 0.655)],
        },
        TestCase {
            name: "CsCl",
            atoms: make_diatomic(55, 17, 2.906),
            expected: vec![(0, 0.769)],
        },
        TestCase {
            name: "CsBr",
            atoms: make_diatomic(55, 35, 3.072),
            expected: vec![(0, 0.767)],
        },
        TestCase {
            name: "CsI",
            atoms: make_diatomic(55, 53, 3.315),
            expected: vec![(0, 0.763)],
        },
        TestCase {
            name: "CO2",
            atoms: vec![
                Atom {
                    atomic_number: 6,
                    position: [0.0, 0.0, 0.0],
                },
                Atom {
                    atomic_number: 8,
                    position: [1.160, 0.0, 0.0],
                },
                Atom {
                    atomic_number: 8,
                    position: [-1.160, 0.0, 0.0],
                },
            ],
            expected: vec![(1, -0.45)],
        },
        TestCase {
            name: "Ketene",
            atoms: {
                let r_co = 1.160;
                let r_cc = 1.314;
                let r_ch = 1.079;
                let angle_hch = 122.0f64.to_radians();
                let half_angle = angle_hch / 2.0;

                vec![
                    Atom {
                        atomic_number: 6,
                        position: [0.0, 0.0, 0.0],
                    },
                    Atom {
                        atomic_number: 8,
                        position: [r_co, 0.0, 0.0],
                    },
                    Atom {
                        atomic_number: 6,
                        position: [-r_cc, 0.0, 0.0],
                    },
                    Atom {
                        atomic_number: 1,
                        position: [
                            -r_cc - r_ch * half_angle.cos(),
                            r_ch * half_angle.sin(),
                            0.0,
                        ],
                    },
                    Atom {
                        atomic_number: 1,
                        position: [
                            -r_cc - r_ch * half_angle.cos(),
                            -r_ch * half_angle.sin(),
                            0.0,
                        ],
                    },
                ]
            },
            expected: vec![(1, -0.45), (0, 0.42), (2, -0.23)],
        },
    ];

    run_group_test("Alkali Halides & Others", cases, 0.01, 0.2);
}
