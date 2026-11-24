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
