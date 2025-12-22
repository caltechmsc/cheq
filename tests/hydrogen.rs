mod common;

use cheq::Atom;
use common::{TestCase, run_group_test};
use std::f64::consts::PI;

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

fn make_water_like(z_center: u8, z_outer: u8, r: f64, angle_deg: f64) -> Vec<Atom> {
    let angle_rad = angle_deg.to_radians();
    let half_angle = angle_rad / 2.0;
    vec![
        Atom {
            atomic_number: z_center,
            position: [0.0, 0.0, 0.0],
        },
        Atom {
            atomic_number: z_outer,
            position: [r * half_angle.cos(), r * half_angle.sin(), 0.0],
        },
        Atom {
            atomic_number: z_outer,
            position: [r * half_angle.cos(), -r * half_angle.sin(), 0.0],
        },
    ]
}

fn make_pyramidal_xy3(z_center: u8, z_outer: u8, r: f64, angle_hxh_deg: f64) -> Vec<Atom> {
    let theta_rad = angle_hxh_deg.to_radians();
    let d_hh = 2.0 * r * (theta_rad / 2.0).sin();
    let big_r = d_hh / 3.0f64.sqrt();
    let z_h = (r.powi(2) - big_r.powi(2)).sqrt();

    let z_coord = if angle_hxh_deg >= 119.99 { 0.0 } else { -z_h };

    vec![
        Atom {
            atomic_number: z_center,
            position: [0.0, 0.0, 0.0],
        },
        Atom {
            atomic_number: z_outer,
            position: [big_r, 0.0, z_coord],
        },
        Atom {
            atomic_number: z_outer,
            position: [-big_r * 0.5, big_r * 3.0f64.sqrt() * 0.5, z_coord],
        },
        Atom {
            atomic_number: z_outer,
            position: [-big_r * 0.5, -big_r * 3.0f64.sqrt() * 0.5, z_coord],
        },
    ]
}

fn make_tetrahedral_xy4(z_center: u8, z_outer: u8, r: f64) -> Vec<Atom> {
    let s = r / 3.0f64.sqrt();
    vec![
        Atom {
            atomic_number: z_center,
            position: [0.0, 0.0, 0.0],
        },
        Atom {
            atomic_number: z_outer,
            position: [s, s, s],
        },
        Atom {
            atomic_number: z_outer,
            position: [s, -s, -s],
        },
        Atom {
            atomic_number: z_outer,
            position: [-s, s, -s],
        },
        Atom {
            atomic_number: z_outer,
            position: [-s, -s, s],
        },
    ]
}

#[test]
fn test_hydrogen_molecules_group() {
    let cases = vec![
        TestCase {
            name: "HF",
            atoms: make_diatomic(1, 9, 0.917),
            expected: vec![(0, 0.46)],
        },
        TestCase {
            name: "H2O",
            atoms: make_water_like(8, 1, 0.958, 104.5),
            expected: vec![(1, 0.35)],
        },
        TestCase {
            name: "NH3",
            atoms: make_pyramidal_xy3(7, 1, 1.012, 106.7),
            expected: vec![(1, 0.24)],
        },
        TestCase {
            name: "CH4",
            atoms: make_tetrahedral_xy4(6, 1, 1.087),
            expected: vec![(1, 0.15)],
        },
        TestCase {
            name: "C2H6",
            atoms: {
                let r_cc = 1.535;
                let r_ch = 1.094;
                let angle_hcc = 111.2f64.to_radians();

                let z_h1 = -r_ch * angle_hcc.cos();
                let r_h1 = r_ch * angle_hcc.sin();

                let z_h2 = r_cc + r_ch * angle_hcc.cos();

                let mut atoms = vec![
                    Atom {
                        atomic_number: 6,
                        position: [0.0, 0.0, 0.0],
                    },
                    Atom {
                        atomic_number: 6,
                        position: [0.0, 0.0, r_cc],
                    },
                ];

                for i in 0..3 {
                    let phi = (i as f64 * 120.0).to_radians();
                    atoms.push(Atom {
                        atomic_number: 1,
                        position: [r_h1 * phi.cos(), r_h1 * phi.sin(), z_h1],
                    });
                }

                for i in 0..3 {
                    let phi = (i as f64 * 120.0 + 60.0).to_radians();
                    atoms.push(Atom {
                        atomic_number: 1,
                        position: [r_h1 * phi.cos(), r_h1 * phi.sin(), z_h2],
                    });
                }
                atoms
            },
            expected: vec![(2, 0.16)],
        },
        TestCase {
            name: "C2H2",
            atoms: vec![
                Atom {
                    atomic_number: 6,
                    position: [0.0, 0.0, 0.0],
                },
                Atom {
                    atomic_number: 6,
                    position: [1.203, 0.0, 0.0],
                },
                Atom {
                    atomic_number: 1,
                    position: [-1.061, 0.0, 0.0],
                },
                Atom {
                    atomic_number: 1,
                    position: [1.203 + 1.061, 0.0, 0.0],
                },
            ],
            expected: vec![(2, 0.13)],
        },
        TestCase {
            name: "C2H4",
            atoms: {
                let r_cc = 1.339;
                let r_ch = 1.086;
                let _angle_hch = 117.6f64.to_radians();
                let angle_cch = 121.2f64.to_radians();

                vec![
                    Atom {
                        atomic_number: 6,
                        position: [0.0, 0.0, 0.0],
                    },
                    Atom {
                        atomic_number: 6,
                        position: [r_cc, 0.0, 0.0],
                    },
                    Atom {
                        atomic_number: 1,
                        position: [-r_ch * angle_cch.cos(), r_ch * angle_cch.sin(), 0.0],
                    },
                    Atom {
                        atomic_number: 1,
                        position: [-r_ch * angle_cch.cos(), -r_ch * angle_cch.sin(), 0.0],
                    },
                    Atom {
                        atomic_number: 1,
                        position: [r_cc + r_ch * angle_cch.cos(), r_ch * angle_cch.sin(), 0.0],
                    },
                    Atom {
                        atomic_number: 1,
                        position: [r_cc + r_ch * angle_cch.cos(), -r_ch * angle_cch.sin(), 0.0],
                    },
                ]
            },
            expected: vec![(2, 0.15)],
        },
        TestCase {
            name: "C6H6",
            atoms: {
                let r_cc = 1.397;
                let r_ch = 1.084;
                let mut atoms = Vec::new();
                for i in 0..6 {
                    let phi = (i as f64 * 60.0).to_radians();
                    atoms.push(Atom {
                        atomic_number: 6,
                        position: [r_cc * phi.cos(), r_cc * phi.sin(), 0.0],
                    });
                    atoms.push(Atom {
                        atomic_number: 1,
                        position: [(r_cc + r_ch) * phi.cos(), (r_cc + r_ch) * phi.sin(), 0.0],
                    });
                }
                atoms
            },
            expected: vec![(1, 0.10)],
        },
        TestCase {
            name: "H2CO",
            atoms: {
                let r_co = 1.205;
                let r_ch = 1.111;
                let _angle_hch = 116.5f64.to_radians();
                let angle_och = 121.75f64.to_radians();

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
                        atomic_number: 1,
                        position: [r_ch * angle_och.cos(), r_ch * angle_och.sin(), 0.0],
                    },
                    Atom {
                        atomic_number: 1,
                        position: [r_ch * angle_och.cos(), -r_ch * angle_och.sin(), 0.0],
                    },
                ]
            },
            expected: vec![(1, -0.43), (0, 0.19), (2, 0.12)],
        },
        TestCase {
            name: "H3COH",
            atoms: {
                let r_co = 1.427;
                let r_oh = 0.956;
                let r_ch_trans = 1.096;
                let r_ch_gauche = 1.094;
                let angle_coh = 108.9f64.to_radians();

                let h_o_x = r_co + r_oh * (PI - angle_coh).cos();
                let h_o_y = r_oh * (PI - angle_coh).sin();

                let angle_hco = 109.5f64.to_radians();

                let h_trans_x = r_ch_trans * angle_hco.cos();
                let h_trans_y = -r_ch_trans * angle_hco.sin();

                let cos120 = -0.5f64;
                let sin120 = 3.0f64.sqrt() / 2.0;

                let h_g1_x = r_ch_gauche * angle_hco.cos();
                let h_g1_y = (-r_ch_gauche * angle_hco.sin()) * cos120;
                let h_g1_z = (-r_ch_gauche * angle_hco.sin()) * sin120;

                let h_g2_x = r_ch_gauche * angle_hco.cos();
                let h_g2_y = (-r_ch_gauche * angle_hco.sin()) * cos120;
                let h_g2_z = -h_g1_z;

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
                        atomic_number: 1,
                        position: [h_o_x, h_o_y, 0.0],
                    },
                    Atom {
                        atomic_number: 1,
                        position: [h_trans_x, h_trans_y, 0.0],
                    },
                    Atom {
                        atomic_number: 1,
                        position: [h_g1_x, h_g1_y, h_g1_z],
                    },
                    Atom {
                        atomic_number: 1,
                        position: [h_g2_x, h_g2_y, h_g2_z],
                    },
                ]
            },
            expected: vec![(2, 0.36), (1, -0.66), (0, -0.15), (4, 0.14), (3, 0.18)],
        },
        TestCase {
            name: "HOC(O)H",
            atoms: {
                let r_co_double = 1.202;
                let r_co_single = 1.343;
                let r_oh = 0.972;
                let r_ch = 1.097;
                let angle_oco = 125.0f64.to_radians();
                let angle_hco_double = 124.1f64.to_radians();

                let o_double_pos = [r_co_double, 0.0, 0.0];

                let o_single_pos = [
                    r_co_single * angle_oco.cos(),
                    r_co_single * angle_oco.sin(),
                    0.0,
                ];

                let h_c_pos = [
                    r_ch * angle_hco_double.cos(),
                    -r_ch * angle_hco_double.sin(),
                    0.0,
                ];

                let angle_oh_global = (125.0f64 + 180.0 - 106.3).to_radians();
                let h_o_pos = [
                    o_single_pos[0] + r_oh * angle_oh_global.cos(),
                    o_single_pos[1] + r_oh * angle_oh_global.sin(),
                    0.0,
                ];

                vec![
                    Atom {
                        atomic_number: 6,
                        position: [0.0, 0.0, 0.0],
                    },
                    Atom {
                        atomic_number: 8,
                        position: o_double_pos,
                    },
                    Atom {
                        atomic_number: 8,
                        position: o_single_pos,
                    },
                    Atom {
                        atomic_number: 1,
                        position: h_c_pos,
                    },
                    Atom {
                        atomic_number: 1,
                        position: h_o_pos,
                    },
                ]
            },
            expected: vec![(1, -0.44), (0, 0.56), (3, 0.16), (2, -0.65), (4, 0.38)],
        },
        TestCase {
            name: "H3CCN",
            atoms: {
                let r_cc = 1.458;
                let r_cn = 1.157;
                let r_ch = 1.104;
                let angle_hcc = 109.5f64.to_radians();

                let h_x = r_ch * angle_hcc.cos();
                let h_r = r_ch * angle_hcc.sin();

                let mut atoms = vec![
                    Atom {
                        atomic_number: 6,
                        position: [0.0, 0.0, 0.0],
                    },
                    Atom {
                        atomic_number: 6,
                        position: [r_cc, 0.0, 0.0],
                    },
                    Atom {
                        atomic_number: 7,
                        position: [r_cc + r_cn, 0.0, 0.0],
                    },
                ];

                for i in 0..3 {
                    let phi = (i as f64 * 120.0).to_radians();
                    atoms.push(Atom {
                        atomic_number: 1,
                        position: [h_x, h_r * phi.cos(), h_r * phi.sin()],
                    });
                }
                atoms
            },
            expected: vec![(2, -0.24), (1, 0.22), (0, -0.37), (3, 0.13)],
        },
        TestCase {
            name: "SiH4",
            atoms: make_tetrahedral_xy4(14, 1, 1.480),
            expected: vec![(1, 0.13)],
        },
        TestCase {
            name: "PH3",
            atoms: make_pyramidal_xy3(15, 1, 1.421, 93.3),
            expected: vec![(1, 0.08)],
        },
        TestCase {
            name: "HCl",
            atoms: make_diatomic(1, 17, 1.275),
            expected: vec![(0, 0.32)],
        },
    ];

    run_group_test("Small Molecules with Hydrogen", cases, 0.035, 0.20);
}
