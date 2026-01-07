# Basis Set Flattening Logic for MD Matrices

The `GenerateMDMatrices` subroutine transforms the basis set structure into a flat, row-indexed format required for the MacMurchie-Davidson integration scheme.

### 1. Cartesian Exponents Matrix
Defines the angular momentum state $(l, m, n)$ for each basis function $i$.

| Row ($b\_idx$) | Type | $l$ (x) | $m$ (y) | $n$ (z) | Physical Meaning |
| :--- | :--- | :--- | :--- | :--- | :--- |
| 1 | s | 0.0 | 0.0 | 0.0 | Shell 1, AO 1 (1s - Core) |
| 2 | s | 0.0 | 0.0 | 0.0 | Shell 1, AO 2 (2s - Inner) |
| 3 | s | 0.0 | 0.0 | 0.0 | Shell 1, AO 3 (3s - Outer) |
| 4 | p | 1.0 | 0.0 | 0.0 | Shell 2, AO 1 (Inner), px |
| 5 | p | 0.0 | 1.0 | 0.0 | Shell 2, AO 1 (Inner), py |
| 6 | p | 0.0 | 0.0 | 1.0 | Shell 2, AO 1 (Inner), pz |
| 7 | p | 1.0 | 0.0 | 0.0 | Shell 2, AO 2 (Outer), px |

### 2. Primitive Exponents Matrix
Contains Gaussian exponents ($\alpha$). Note that orbitals within the same shell share the same exponent pool.

| Row ($i$) | Physical Meaning | Col 1 ($\alpha_1$) | Col 2 ($\alpha_2$) | ... | Col 10 ($\alpha_{10}$) |
| :--- | :--- | :--- | :--- | :--- | :--- |
| 1 | 1s (Core) | 3047.52 | 457.36 | ... | 0.1687 |
| 2 | 2s (Inner) | 3047.52 | 457.36 | ... | 0.1687 |
| 3 | 3s (Outer) | 3047.52 | 457.36 | ... | 0.1687 |
| 4 | 2px (Inner) | 7.8682 | 1.8812 | ... | 0.0000 |

### 3. Contraction Coefficients Matrix
Contains the raw contraction coefficients ($d$). This differentiates between core and valence orbitals sharing the same exponents.

| Row ($i$) | Physical Meaning | Col 1 ($d_1$) | Col 2 ($d_2$) | ... | Col 10 ($d_{10}$) |
| :--- | :--- | :--- | :--- | :--- | :--- |
| 1 | 1s (Core) | 0.00183 | 0.01403 | ... | 0.00000 |
| 2 | 2s (Inner) | 0.00000 | 0.00000 | ... | 0.00000 |
| 3 | 3s (Outer) | 0.00000 | 0.00000 | ... | 1.00000 |
| 4 | 2px (Inner) | 0.06899 | 0.31642 | ... | 0.00000 |