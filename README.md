# Cheq

**Cheq** is a high-performance, foundational software library for calculating dynamic partial atomic charges using the Charge Equilibration (QEq) method. It provides a modern, robust solution for translating molecular geometry and elemental properties into self-consistent charge distributions essential for accurate molecular dynamics simulations. This library is engineered from the ground up in **Rust** for exceptional performance, memory safety, and strict adherence to the principles of modular software design.

The core mission of Cheq is to provide a reliable, predictable, and easy-to-integrate tool for developers and researchers building the next generation of simulation tools for general chemistry, materials science, and biological systems.

## Features

- **Dynamic Charge Calculation**: Implements the canonical Rappé & Goddard (1991) QEq method.
- **Geometry Dependent**: Charges respond naturally to changes in molecular conformation.
- **Decoupled Architecture**: Agnostic to basic data structures via the flexible `AtomView` trait.
- **Memory Safe & Fast**: Built in Rust with optimized linear algebra for high performance.
- **Configurable Parameters**: Includes standard parameters with support for custom TOML-based sets.

## Getting Started

To get started with Cheq, add it as a dependency in your `Cargo.toml`:

```toml
[dependencies]
cheq = "0.1.0"
```

Then, you can use it in your Rust code as follows:

```rust
use cheq::{get_default_parameters, QEqSolver, Atom, SolverOptions};
use approx::assert_relative_eq;

fn main() {
    // 1. Set up the parameters and solver.
    let params = get_default_parameters();
    // Use default solver options (can be customized for tolerance/iterations)
    let solver = QEqSolver::new(params);

    // 2. Define the molecular system (e.g., a Water molecule).
    // Cheq is agnostic to your data structure; here we use the simple `Atom` struct.
    let bond_length = 0.9575; // Angstroms
    let angle_rad = 104.45f64.to_radians();

    let atoms = vec![
        Atom { atomic_number: 8, position: [0.0, 0.0, 0.0] },         // Oxygen
        Atom { atomic_number: 1, position: [bond_length, 0.0, 0.0] }, // Hydrogen 1
        Atom { atomic_number: 1, position: [
            bond_length * angle_rad.cos(),
            bond_length * angle_rad.sin(),
            0.0,
        ]},                                                           // Hydrogen 2
    ];

    // 3. Solve for partial charges (imposing a total charge of 0.0).
    let result = solver.solve(&atoms, 0.0).expect("QEq calculation failed to converge");

    // 4. Inspect the results.
    println!("Equilibrated Chemical Potential: {:.4} eV", result.equilibrated_potential);
    println!("Charges: O={:.4}, H1={:.4}, H2={:.4}",
        result.charges[0], result.charges[1], result.charges[2]);

    // Verify physics:
    assert!(result.charges[0] < 0.0); // Oxygen should be negative
    assert_relative_eq!(result.charges.iter().sum::<f64>(), 0.0, epsilon = 1e-9); // Conservation
}
```

> **Note**: This is a simplified example. For more advanced usage, including custom parameter loading and integration with existing molecular data structures, please refer to the [API Documentation](https://docs.rs/cheq).

## Documentation

- [API Documentation](https://docs.rs/cheq) - Comprehensive reference for all public types and functions.
- [Theory Background](https://pubs.acs.org/doi/10.1021/j100161a070) - The original QEq paper by Rappé and Goddard (J. Phys. Chem. 1991, 95, 3358-3363).

## Tech Stack

- **Core Language**: Rust
- **Linear Algebra**: `faer`
- **Serialization**: `serde` & `toml`

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
