mod common;

use cheq::Atom;
use common::{TestCase, run_group_test};

fn make_nacl_monomer(dist: f64) -> Vec<Atom> {
    vec![
        Atom {
            atomic_number: 11,
            position: [0.0, 0.0, 0.0],
        },
        Atom {
            atomic_number: 17,
            position: [dist, 0.0, 0.0],
        },
    ]
}

fn make_nacl_dimer_square(dist: f64) -> Vec<Atom> {
    vec![
        Atom {
            atomic_number: 11,
            position: [0.0, 0.0, 0.0],
        },
        Atom {
            atomic_number: 17,
            position: [dist, 0.0, 0.0],
        },
        Atom {
            atomic_number: 17,
            position: [0.0, dist, 0.0],
        },
        Atom {
            atomic_number: 11,
            position: [dist, dist, 0.0],
        },
    ]
}

fn make_nacl_tetramer_cube(dist: f64) -> Vec<Atom> {
    vec![
        Atom {
            atomic_number: 11,
            position: [0.0, 0.0, 0.0],
        },
        Atom {
            atomic_number: 17,
            position: [dist, 0.0, 0.0],
        },
        Atom {
            atomic_number: 17,
            position: [0.0, dist, 0.0],
        },
        Atom {
            atomic_number: 11,
            position: [dist, dist, 0.0],
        },
        Atom {
            atomic_number: 17,
            position: [0.0, 0.0, dist],
        },
        Atom {
            atomic_number: 11,
            position: [dist, 0.0, dist],
        },
        Atom {
            atomic_number: 11,
            position: [0.0, dist, dist],
        },
        Atom {
            atomic_number: 17,
            position: [dist, dist, dist],
        },
    ]
}

#[test]
fn test_clusters_group() {
    let d = 2.361;

    let cases = vec![
        TestCase {
            name: "NaCl Monomer",
            atoms: make_nacl_monomer(d),
            expected: vec![(0, 0.749), (1, -0.749)],
        },
        TestCase {
            name: "(NaCl)2 Square",
            atoms: make_nacl_dimer_square(d),
            expected: vec![(0, 0.826), (1, -0.826)],
        },
        TestCase {
            name: "(NaCl)4 Cube",
            atoms: make_nacl_tetramer_cube(d),
            expected: vec![(0, 0.823), (4, -0.823)],
        },
    ];

    run_group_test("Ionic Clusters", cases, 0.05, 0.10);
}
