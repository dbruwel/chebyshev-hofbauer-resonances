# Computing Resonances of Chaotic Systems via Chebyshev Polynomials and Hofbauer Extensions

This project was undertaken during the summer of 2024/2025 under the supervision of Dr. Caroline Wormell at the University of Sydney.
It was a short research project focused primarily on the conceptual development of the technique. This repository is not a polished Python package, nor is the accompanying paper complete or peer-reviewed. It is simply a collection of exploratory ideas, many of which remain unrefined or incomplete. The code is provided as part of the development process and is not intended for reuse by others.

This project was completed by myself (Daniel Bruwel), and Owen Marscgner, and Adrian Zhao.

## Project Overview

We study the resonances (eigenvalues) of the transfer operator associated with chaotic dynamical systems. For a one-dimensional discrete-time system defined by

$$
x_{t+1} = f(x_t),
$$

if $x_t \sim \psi$, then the distribution at time $t+1$ is given by $\mathcal{L} \psi$, where $\mathcal{L}$ is the transfer operator.

Eigenfunctions of $\mathcal{L}$ evolve only by scaling:
- If $|\lambda| < 1$, the corresponding distribution decays over time.
- If $|\lambda| = 1$, it corresponds to a stable mode of the system.

This project develops numerical techniques for computing the spectrum of $\mathcal{L}$, particularly in non-Markovian cases. The method uses Chebyshev polynomial approximations and Hofbauer extensions to achieve stable and accurate results.

---

## Installation

Clone the repository and install in editable mode:

```bash
git clone https://github.com/dbruwel/chebyshev-hofbauer-resonances.git
cd chebyshev-hofbauer-resonances
pip install -e .
```

---

## Repository Structure

```
chebyshev-hofbauer-resonances/
│
├── notebooks/                # Jupyter notebooks for exploration and results
│   └── *.ipynb
│
├── paper/                    # Draft paper, figures, and supplementary material
│   ├── main.pdf
│   └── ...
│
├── src/
│   └── chebyshev_hofbauer_resonances/
│       ├── __init__.py
│       └── ...
│
├── setup.py
├── pyproject.toml
└── README.md
```

---

## Usage

Import core functionality from the package:

```python
from chebyshev_hofbauer_resonances.operator_approx import cheb_op_ap
```

More detailed examples can be found in the `notebooks/` directory.

---

## License

MIT License. See `LICENSE` for details.
