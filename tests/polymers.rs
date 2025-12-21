mod common;

use cheq::Atom;
use common::{TestCase, run_group_test};

struct Turtle {
    position: [f64; 3],
    heading: [f64; 3],
    up: [f64; 3],
    left: [f64; 3],
}

impl Turtle {
    fn new() -> Self {
        Self {
            position: [0.0, 0.0, 0.0],
            heading: [1.0, 0.0, 0.0],
            up: [0.0, 0.0, 1.0],
            left: [0.0, 1.0, 0.0],
        }
    }

    fn forward(&mut self, dist: f64) {
        self.position[0] += self.heading[0] * dist;
        self.position[1] += self.heading[1] * dist;
        self.position[2] += self.heading[2] * dist;
    }

    fn roll(&mut self, angle_deg: f64) {
        let rad = angle_deg.to_radians();
        let (s, c) = rad.sin_cos();

        let new_up = [
            self.up[0] * c + self.left[0] * s,
            self.up[1] * c + self.left[1] * s,
            self.up[2] * c + self.left[2] * s,
        ];

        let new_left = [
            -self.up[0] * s + self.left[0] * c,
            -self.up[1] * s + self.left[1] * c,
            -self.up[2] * s + self.left[2] * c,
        ];

        self.up = new_up;
        self.left = new_left;
        self.normalize();
    }

    fn pitch(&mut self, angle_deg: f64) {
        let rad = angle_deg.to_radians();
        let (s, c) = rad.sin_cos();

        let new_heading = [
            self.heading[0] * c + self.up[0] * s,
            self.heading[1] * c + self.up[1] * s,
            self.heading[2] * c + self.up[2] * s,
        ];

        let new_up = [
            -self.heading[0] * s + self.up[0] * c,
            -self.heading[1] * s + self.up[1] * c,
            -self.heading[2] * s + self.up[2] * c,
        ];

        self.heading = new_heading;
        self.up = new_up;
        self.normalize();
    }

    fn normalize(&mut self) {
        let mag_h =
            (self.heading[0].powi(2) + self.heading[1].powi(2) + self.heading[2].powi(2)).sqrt();
        self.heading.iter_mut().for_each(|x| *x /= mag_h);

        let mag_u = (self.up[0].powi(2) + self.up[1].powi(2) + self.up[2].powi(2)).sqrt();
        self.up.iter_mut().for_each(|x| *x /= mag_u);

        let mag_l = (self.left[0].powi(2) + self.left[1].powi(2) + self.left[2].powi(2)).sqrt();
        self.left.iter_mut().for_each(|x| *x /= mag_l);
    }

    fn extend(&mut self, r: f64, theta: f64, phi: f64) {
        self.roll(phi);
        self.pitch(180.0 - theta);
        self.forward(r);
    }

    fn save(&self) -> Self {
        Self {
            position: self.position,
            heading: self.heading,
            up: self.up,
            left: self.left,
        }
    }

    fn restore(&mut self, state: &Self) {
        self.position = state.position;
        self.heading = state.heading;
        self.up = state.up;
        self.left = state.left;
    }
}

fn make_polymer_chain(n_units: usize, is_pvdf: bool) -> Vec<Atom> {
    let mut atoms = Vec::new();
    let mut turtle = Turtle::new();

    atoms.push(Atom {
        atomic_number: 6,
        position: turtle.position,
    });

    let mut c_states = Vec::new();
    c_states.push(turtle.save());

    let r_cc = 1.54;
    let angle_ccc = 109.5;
    let torsion_trans = 180.0;

    for _ in 1..n_units {
        turtle.extend(r_cc, angle_ccc, torsion_trans);
        atoms.push(Atom {
            atomic_number: 6,
            position: turtle.position,
        });
        c_states.push(turtle.save());
    }

    for (i, state) in c_states.iter().enumerate() {
        let mut t = Turtle::new();
        t.restore(state);

        let is_cf2 = is_pvdf && (i % 2 != 0);
        let z_sub = if is_cf2 { 9 } else { 1 };
        let r_sub = if is_cf2 { 1.35 } else { 1.09 };

        let mut t1 = t.save();
        t1.extend(r_sub, 109.5, 60.0);
        atoms.push(Atom {
            atomic_number: z_sub,
            position: t1.position,
        });

        let mut t2 = t.save();
        t2.extend(r_sub, 109.5, -60.0);
        atoms.push(Atom {
            atomic_number: z_sub,
            position: t2.position,
        });

        if i == 0 {
            let mut t_cap = t.save();
            t_cap.extend(1.09, 0.0, 0.0);
            atoms.push(Atom {
                atomic_number: 1,
                position: t_cap.position,
            });
        }
        if i == n_units - 1 {
            let mut t_cap = t.save();
            t_cap.extend(1.09, 109.5, 180.0);
            atoms.push(Atom {
                atomic_number: 1,
                position: t_cap.position,
            });
        }
    }

    atoms
}

fn make_ala_his_ala(protonated: bool) -> Vec<Atom> {
    let mut atoms = Vec::new();
    let mut t = Turtle::new();

    let r_n_ca = 1.46;
    let r_ca_c = 1.51;
    let r_c_n = 1.33;
    let r_ca_cb = 1.54;
    let r_c_o = 1.23;
    let r_n_h = 1.01;

    let a_n_ca_c = 111.0;
    let a_ca_c_n = 116.0;
    let a_c_n_ca = 122.0;

    let phi = -135.0;
    let psi = 135.0;
    let omega = 180.0;

    atoms.push(Atom {
        atomic_number: 7,
        position: t.position,
    });

    let mut t_h = t.save();
    t_h.roll(180.0);
    t_h.pitch(180.0 - 109.5);
    t_h.forward(r_n_h);
    atoms.push(Atom {
        atomic_number: 1,
        position: t_h.position,
    });

    let mut t_h2 = t.save();
    t_h2.roll(180.0 + 120.0);
    t_h2.pitch(180.0 - 109.5);
    t_h2.forward(r_n_h);
    atoms.push(Atom {
        atomic_number: 1,
        position: t_h2.position,
    });

    t.forward(r_n_ca);
    atoms.push(Atom {
        atomic_number: 6,
        position: t.position,
    });
    let ca1_state = t.save();

    let mut t_ha = t.save();
    t_ha.extend(1.09, 109.5, psi + 120.0);
    atoms.push(Atom {
        atomic_number: 1,
        position: t_ha.position,
    });

    let mut t_cb = t.save();
    t_cb.extend(r_ca_cb, 109.5, psi - 120.0);
    atoms.push(Atom {
        atomic_number: 6,
        position: t_cb.position,
    });
    let cb_state = t_cb.save();
    for k in 0..3 {
        let mut t_m = cb_state.save();
        t_m.extend(1.09, 109.5, k as f64 * 120.0 + 60.0);
        atoms.push(Atom {
            atomic_number: 1,
            position: t_m.position,
        });
    }

    t.restore(&ca1_state);
    t.extend(r_ca_c, a_n_ca_c, psi);
    atoms.push(Atom {
        atomic_number: 6,
        position: t.position,
    });

    let mut t_o = t.save();
    t_o.extend(r_c_o, 121.0, 0.0);
    atoms.push(Atom {
        atomic_number: 8,
        position: t_o.position,
    });

    t.extend(r_c_n, a_ca_c_n, omega);
    atoms.push(Atom {
        atomic_number: 7,
        position: t.position,
    });
    let mut t_hn = t.save();
    t_hn.extend(r_n_h, 120.0, 0.0);
    atoms.push(Atom {
        atomic_number: 1,
        position: t_hn.position,
    });
    t.extend(r_n_ca, a_c_n_ca, phi);
    atoms.push(Atom {
        atomic_number: 6,
        position: t.position,
    });
    let ca2_state = t.save();

    let mut t_ha2 = t.save();
    t_ha2.extend(1.09, 109.5, psi + 120.0);
    atoms.push(Atom {
        atomic_number: 1,
        position: t_ha2.position,
    });

    let mut t_cb2 = t.save();
    t_cb2.extend(r_ca_cb, 109.5, psi - 120.0);
    atoms.push(Atom {
        atomic_number: 6,
        position: t_cb2.position,
    });

    let mut t_hb1 = t_cb2.save();
    t_hb1.extend(1.09, 109.5, 120.0);
    atoms.push(Atom {
        atomic_number: 1,
        position: t_hb1.position,
    });
    let mut t_hb2 = t_cb2.save();
    t_hb2.extend(1.09, 109.5, 240.0);
    atoms.push(Atom {
        atomic_number: 1,
        position: t_hb2.position,
    });

    t_cb2.extend(1.50, 109.5, -60.0);
    atoms.push(Atom {
        atomic_number: 6,
        position: t_cb2.position,
    });
    let cg_state = t_cb2.save();

    let mut t_ring = cg_state.save();
    t_ring.extend(1.38, 126.0, 90.0);
    atoms.push(Atom {
        atomic_number: 7,
        position: t_ring.position,
    });
    let nd1_state = t_ring.save();

    t_ring.extend(1.32, 108.0, 180.0);
    atoms.push(Atom {
        atomic_number: 6,
        position: t_ring.position,
    });
    let ce1_state = t_ring.save();

    t_ring.extend(1.32, 108.0, 0.0);
    atoms.push(Atom {
        atomic_number: 7,
        position: t_ring.position,
    });
    let ne2_state = t_ring.save();

    t_ring.extend(1.38, 108.0, 0.0);
    atoms.push(Atom {
        atomic_number: 6,
        position: t_ring.position,
    });
    let cd2_state = t_ring.save();

    let mut t_hd2 = cd2_state.save();
    t_hd2.extend(1.08, 126.0, 180.0);
    atoms.push(Atom {
        atomic_number: 1,
        position: t_hd2.position,
    });

    let mut t_he1 = ce1_state.save();
    t_he1.extend(1.08, 126.0, 180.0);
    atoms.push(Atom {
        atomic_number: 1,
        position: t_he1.position,
    });

    let mut t_hd1 = nd1_state.save();
    t_hd1.extend(1.01, 126.0, 0.0);
    atoms.push(Atom {
        atomic_number: 1,
        position: t_hd1.position,
    });

    if protonated {
        let mut t_he2 = ne2_state.save();
        t_he2.extend(1.01, 126.0, 180.0);
        atoms.push(Atom {
            atomic_number: 1,
            position: t_he2.position,
        });
    }

    t.restore(&ca2_state);
    t.extend(r_ca_c, a_n_ca_c, psi);
    atoms.push(Atom {
        atomic_number: 6,
        position: t.position,
    });
    let c2_state = t.save();

    let mut t_o2 = c2_state.save();
    t_o2.extend(r_c_o, 121.0, 0.0);
    atoms.push(Atom {
        atomic_number: 8,
        position: t_o2.position,
    });

    t.extend(r_c_n, a_ca_c_n, omega);
    atoms.push(Atom {
        atomic_number: 7,
        position: t.position,
    });
    let mut t_hn3 = t.save();
    t_hn3.extend(r_n_h, 120.0, 0.0);
    atoms.push(Atom {
        atomic_number: 1,
        position: t_hn3.position,
    });

    t.extend(r_n_ca, a_c_n_ca, phi);
    atoms.push(Atom {
        atomic_number: 6,
        position: t.position,
    });
    let ca3_state = t.save();

    let mut t_ha3 = t.save();
    t_ha3.extend(1.09, 109.5, psi + 120.0);
    atoms.push(Atom {
        atomic_number: 1,
        position: t_ha3.position,
    });

    let mut t_cb3 = t.save();
    t_cb3.extend(r_ca_cb, 109.5, psi - 120.0);
    atoms.push(Atom {
        atomic_number: 6,
        position: t_cb3.position,
    });
    let cb3_state = t_cb3.save();
    for k in 0..3 {
        let mut t_m = cb3_state.save();
        t_m.extend(1.09, 109.5, k as f64 * 120.0 + 60.0);
        atoms.push(Atom {
            atomic_number: 1,
            position: t_m.position,
        });
    }

    t.restore(&ca3_state);
    t.extend(r_ca_c, a_n_ca_c, psi);
    atoms.push(Atom {
        atomic_number: 6,
        position: t.position,
    });
    let c3_state = t.save();

    let mut t_o3 = c3_state.save();
    t_o3.extend(r_c_o, 121.0, 0.0);
    atoms.push(Atom {
        atomic_number: 8,
        position: t_o3.position,
    });

    let mut t_oh = c3_state.save();
    t_oh.extend(1.34, 115.0, 180.0);
    atoms.push(Atom {
        atomic_number: 8,
        position: t_oh.position,
    });

    let mut t_ho = t_oh.save();
    t_ho.extend(0.97, 109.0, 180.0);
    atoms.push(Atom {
        atomic_number: 1,
        position: t_ho.position,
    });

    atoms
}

#[test]
fn test_polymers_group() {
    let cases = vec![
        TestCase {
            name: "PE (5 units)",
            atoms: make_polymer_chain(5, false),
            expected: vec![(2, -0.284), (10, 0.143)],
        },
        TestCase {
            name: "PVDF (5 units)",
            atoms: make_polymer_chain(5, true),
            expected: vec![(2, -0.082), (10, 0.050), (3, 0.776), (12, -0.415)],
        },
        TestCase {
            name: "Ala-His-Ala (Neutral)",
            atoms: make_ala_his_ala(false),
            expected: vec![
                (0, -0.81),
                (38, -0.45),
                (37, -0.46),
                (11, -0.46),
                (12, 0.27),
                (9, 0.47),
                (10, -0.51),
            ],
        },
        TestCase {
            name: "Ala-His-Ala (Protonated)",
            atoms: make_ala_his_ala(true),
            expected: vec![(19, -0.50), (21, -0.50), (25, 0.33), (26, 0.33)],
        },
    ];

    run_group_test("Polymers & Biomolecules", cases, 0.10, 0.35);
}
