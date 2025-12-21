use cheq::{Atom, QEqSolver, get_default_parameters};

pub struct TestCase<'a> {
    pub name: &'a str,
    pub atoms: Vec<Atom>,
    pub expected: Vec<(usize, f64)>,
}

pub fn run_group_test(
    group_name: &str,
    cases: Vec<TestCase>,
    group_avg_limit: f64,
    group_max_limit: f64,
) {
    let params = get_default_parameters();
    let solver = QEqSolver::new(params);

    let mut group_total_error = 0.0;
    let mut group_max_error = 0.0;
    let mut total_data_points = 0;

    println!("\nRunning Group Test: {}", group_name);
    println!("{:-<80}", "");
    println!(
        "{:<20} | {:<10} | {:<10} | {:<10}",
        "Molecule", "Atom Idx", "Expected", "Calculated"
    );

    for case in cases {
        let result = solver.solve(&case.atoms, 0.0).expect("Solver failed");

        for (index, expected_q) in &case.expected {
            let calculated_q = result.charges[*index];
            let error = (calculated_q - expected_q).abs();

            println!(
                "{:<20} | {:<10} | {:<10.4} | {:<10.4} (Err: {:.4})",
                case.name, index, expected_q, calculated_q, error
            );

            group_total_error += error;
            if error > group_max_error {
                group_max_error = error;
            }
            total_data_points += 1;
        }
    }

    let group_avg_error = if total_data_points > 0 {
        group_total_error / total_data_points as f64
    } else {
        0.0
    };

    println!("{:-<80}", "");
    println!("Group Statistics for '{}':", group_name);
    println!("  Total Data Points: {}", total_data_points);
    println!(
        "  Group Avg Error:   {:.4} (Limit: {:.4})",
        group_avg_error, group_avg_limit
    );
    println!(
        "  Group Max Error:   {:.4} (Limit: {:.4})",
        group_max_error, group_max_limit
    );
    println!("{:-<80}\n", "");

    assert!(
        group_avg_error <= group_avg_limit,
        "Group average error {:.4} exceeds limit {:.4}",
        group_avg_error,
        group_avg_limit
    );

    assert!(
        group_max_error <= group_max_limit,
        "Group maximum error {:.4} exceeds limit {:.4}",
        group_max_error,
        group_max_limit
    );
}
