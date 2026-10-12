# Computing Resonances of Chaotic Systems via Chebyshev Polynomials and Hofbauer Extensions

> **Note:** This repository is a clone of the active working project, shared for the purposes of the Denison scholarship. It reflects the state of the work at the time of submission and is not the finished product.

This project began as a summer research project in 2024/2025 at the University of Sydney, completed by Daniel Bruwel, Owen Marschner, and Adrian Zhao under the supervision of Dr. Caroline Wormell. That initial project focused primarily on the conceptual development of the technique for the linear tent map.

**Ongoing work:** Daniel and Caroline are now extending this work to general tent maps. This repository therefore contains both the original summer project code and work-in-progress code from the extension. The preliminary report under `papers/preliminary_report/` was prepared in fulfilment of the Denison scholarship requirements; it is not a complete or peer-reviewed paper. A more complete write-up will be produced as part of the extended work.

This repository is not a polished Python package. Much of the code is exploratory, some ideas remain unrefined or incomplete, and it is not intended for reuse by others.

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

First, clone the repository and navigate into the directory:

```bash
git clone https://github.com/dbruwel/chebyshev-hofbauer-resonances.git
cd chebyshev-hofbauer-resonances
```

You can install the package using either pip or Poetry. Poetry is recommended for more stable dependency management.

**Option 1: Using pip**

You can install the package directly using pip. Use the -e flag for an "editable" install, which is useful for development.

```bash
# Standard install
pip install .

# Or, for an editable install
pip install -e .
```

**Option 2: Using Poetry (Recommended)**
1. **Install Poetry** (if you don't already have it):
```bash
pip install poetry
```
2. **Install the package and dependencies:**
```bash
poetry install
```

Poetry will automatically install the project into a virtual environment.
* **Using with conda:** If you have an active conda environment, Poetry will detect it and install the dependencies into that environment. This is the recommended approach.
* **Without conda:** If you are not using conda, Poetry will create and manage its own separate virtual environment.

---

## Repository Structure

```
chebyshev-hofbauer-resonances/
│
├── experiments/
│   ├── notebooks/            # Jupyter notebooks for exploration and results
│   │   ├── linear_tent_map/  # Original summer project
│   │   └── general_tent_map/ # Ongoing extension (work in progress)
│   └── scripts/
│
├── papers/
│   └── preliminary_report/   # Denison scholarship report, figures
│       ├── main.pdf
│       └── ...
│
├── src/
│   └── chebyshev_hofbauer_resonances/
│       ├── __init__.py
│       ├── linear_tent_map/
│       └── general_tent_map/
│
├── pyproject.toml
└── README.md
```

---

## Usage

Import core functionality from the package:

```python
from chebyshev_hofbauer_resonances.linear_tent_map.operator_approx import cheb_op_ap
```

More detailed examples can be found in the `experiments/notebooks/` directory.
