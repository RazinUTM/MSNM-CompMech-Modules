# MSNM Computational Mechanics Modules

## P1 — Two springs in series

P1 assembles two linear springs by the direct stiffness method, prescribes the displacement at node 1, and calculates the nodal displacements and reaction. The spring stiffnesses are 2,000 and 500 N/m. A force of 100 N acts at node 3. The prescribed displacement is 0 m in the first case and 0.01 m in the second.

### Run in VS Code

Open the repository folder and run these commands from its root. For MATLAB, use the MATLAB Command Window or run the commands from the integrated terminal with `matlab -batch "COMMAND"`.

```matlab
run('problems/SpringSystem/MAIN_SpringSystem.m')
addpath('tests','tests/kernels')
TEST_SV_Prescribed
TEST_SpringSystem
```

For Python, use Python 3.12 or later and select the interpreter for this project in VS Code.

```sh
python -m pip install -r requirements.txt
python problems/SpringSystem/main_spring_system.py
python -m pytest -q tests
```

Both drivers print the nodal displacements, reaction, spring extensions and forces, and the largest nodal equilibrium residual. The four plot panels show the displacement for each support case, the spring forces, and the applied forces and reactions.

### Check the results

The end load puts both springs in tension at 100 N. Their extensions are `100/2000 = 0.05` m and `100/500 = 0.20` m. Add these extensions in sequence to the prescribed displacement to obtain the nodal values.

| Prescribed displacement (m) | Node 2 (m) | Node 3 (m) | Reaction at node 1 (N) |
|---:|---:|---:|---:|
| 0 | 0.05 | 0.25 | −100 |
| 0.01 | 0.06 | 0.26 | −100 |

The partitioned solution should agree with these values to floating-point precision. In the second case, dropping the support correction `-K(f,c)*U(c)` leaves the free displacements 0.01 m too low. The penalty method gives a close approximation and raises the condition number from about 6.85 to 6.55 million.

### Export the results

To save a PNG and two CSV files, set `P1_OUTPUT_DIR` to an absolute directory before running a driver. For example, in MATLAB at the repository root, run `setenv('P1_OUTPUT_DIR',fullfile(pwd,'results'))`. In a Bash terminal, run `export P1_OUTPUT_DIR="$PWD/results"`. In PowerShell, run `$env:P1_OUTPUT_DIR = Join-Path (Get-Location).Path 'results'`.

Each nodal CSV has one row per node and support case. The residual is the internal nodal force minus the applied force and reaction. Each spring CSV has one row per spring and support case. Positive extension and axial force indicate tension for the stated node order. Files use `_matlab` or `_python` suffixes and are replaced on rerun.

### Shared solver

The P1 drivers call the partitioned solver in `core/solver/`. In Python, use
`solve_prescribed(K, F, constrained, prescribed)` with zero-based DOF indices.
In MATLAB, use `[U,R] = SV_Prescribed(K,F,c,values)` with one-based indices.
Each call handles one load case. The solver computes the free displacements from
`Kff Uf = Ff - Kfc Uc` and returns reactions in the same order as the constrained
DOFs. Assembly, the closed-form comparison and plots are part of the P1 example.

### Provenance

The P1 drivers and validation checks were adapted from [MSNM Numerical Modules 3 at commit 768eb88](https://github.com/RazinUTM/MSNM-Numerical-Modules-3/tree/768eb886069900ef318701a6a9776beb10d2a6db). This version adds post-processing and runs from this checkout. The original authorship notices and [MIT license](LICENSE) accompany the code.

## P2 — Axial bar by 1D finite elements

The example uses an axially loaded bar of length 2 m with a linearly tapered
area, a distributed axial load of 50 kN/m, and an end force of 50 kN. Start
with eight linear elements. The separate convergence study compares
linear and quadratic elements on five meshes, for uniform and tapered bars.

```matlab
run('problems/AxialBar/MAIN_AxialBar.m')
run('studies/AxialBar/STUDY_AxialBarConvergence.m')
```

```sh
python problems/AxialBar/main_axial_bar.py
python studies/AxialBar/study_axial_bar_convergence.py
python -m pytest -q tests
```

To run both P2 MATLAB checks from the VS Code terminal at the repository root:

```sh
matlab -batch "addpath(fullfile(pwd,'tests')); TEST_AxialBar; TEST_AxialBarStudy"
```

### Follow the calculation

1. `SF_FEMLG1D` and `GQ_Points` in `core/` supply the 1D Lagrange functions and
   Gauss quadrature. `KM_Line1D` integrates element matrices and loads.
2. `AS_Global` accumulates the element matrices by DOF number. P1 also uses
   this routine. The shared prescribed-displacement solver returns the
   displacements and support reaction.
3. `SF_Line1D` evaluates the solution and its derivative in each element.
   `UT_AxialBar` supplies the closed form. `PL_AxialBar` plots the nodal
   displacements and the axial force on each element against that reference.

For the eight-element example, the free-end displacement is about
**1.30734 mm**, the exact value is **1.30685 mm**, and the reaction is
**−150 kN**. The relative displacement L2 error is about `7.504e-4`
and the relative energy error is about `1.344e-2`. For the tapered bar,
the study gives L2 rates near 2 and 3 and energy rates near 1 and 2 for
linear and quadratic elements, respectively.

Set `P2_OUTPUT_DIR` to an absolute directory before running the main example
to save the displacement and axial force CSVs and the profile plot PNG. Each
force row includes its element number so that the two element-side values at
a shared node can be distinguished. As with P1, each language uses its own
filename suffix.

The P2 example, study, numerical kernels and validations were adapted from
[MSNM Numerical Modules 3 at commit 7a17b00](https://github.com/RazinUTM/MSNM-Numerical-Modules-3/tree/7a17b007fc6b9f83790b6f8b7b4520a8c05629a2).
MATLAB results in the source headers refer to validation in MSNM3.
