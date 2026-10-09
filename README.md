# Hybrid Quantum–Classical Grid Optimisation Using QAOA

### Renewable-Aware Unit Commitment and Power Dispatch

This project investigates a hybrid quantum–classical approach to **unit commitment**, the problem of determining which electricity generators should operate, and at what output, to satisfy demand at minimum operating cost.

Using the **PGLib-UC RTS-GMLC** benchmark, we implement classical optimisation models alongside Qiskit-based **Quantum Approximate Optimization Algorithm (QAOA)** and **Conditional Value at Risk QAOA (CVaR-QAOA)**.

The project progresses from a six-qubit, single-period demonstration to a 24-period renewable-aware scheduling model, with a 12-qubit quantum experiment on a smaller time-dependent subproblem.

**The objective is to study quantum optimisation methods and evaluate their performance against classical baselines—not to claim quantum computational advantage.**

---

## 1. Problem Statement

Electricity grids must continuously balance generation and consumption. This becomes more challenging with renewable energy because wind and solar generation vary over time.

Thermal generators introduce additional constraints, including:

- Minimum and maximum generation capacity
- Startup costs and production costs
- Ramp-up and ramp-down limits
- Minimum ON/OFF durations
- Reserve capacity requirements

The objective is to minimise total operating cost:

\[
\min \sum_{t}\sum_i \left[C_i(P_{i,t})+S_i u_{i,t}\right]
\]

subject to power balance, reserve requirements, generator operating limits, and applicable scheduling constraints.

Here, \(P_{i,t}\) is generator output, \(C_i\) is its production cost, and \(u_{i,t}\) indicates a startup event.

The problem combines **binary commitment variables** with **continuous generation decisions**, making it a mixed-integer optimisation problem.

## 2. Dataset

We use the [PGLib-UC benchmark library](https://github.com/power-grid-lib/pglib-uc), specifically the [RTS-GMLC instance for 2020-01-27](https://github.com/power-grid-lib/pglib-uc/blob/master/rts_gmlc/2020-01-27.json).

The original instance contains:

| Component | Description |
|---|---|
| Thermal generators | 73 |
| Renewable generators | 81 |
| Scheduling periods | 48 |
| Demand | Time-dependent electricity demand |
| Reserves | Time-dependent reserve requirements |
| Generator parameters | Capacity, ramp limits, operating constraints and costs |

We construct smaller benchmark-derived scenarios by selecting subsets of generators and scaling demand and reserves according to selected thermal capacity. Renewable availability profiles are retained for the selected units.

These reduced scenarios are **constructed test instances**, not geographically equivalent miniature versions of the original grid.

## 3. Methodology

### Experiment A: Single-Period Unit Commitment

The initial experiment uses six thermal generators, one wind generator, and one solar generator for a single scheduling period.

With six binary commitment decisions, there are:

\[
2^6=64
\]

possible thermal ON/OFF configurations.

**Classical baseline**

Each commitment is evaluated using SciPy's linear programming solver for continuous economic dispatch. The lowest-cost feasible configuration is the exact optimum of the simplified single-period model.

**Quantum optimisation**

A six-variable QUBO is constructed using approximate generator costs and a quadratic capacity-adequacy penalty. It is converted into an Ising Hamiltonian and optimised using Qiskit.

We compare:

- Standard QAOA, minimising expected QUBO energy
- CVaR-QAOA, minimising the lowest-energy 20% of probability mass
- Uniform random sampling

The QUBO is a surrogate objective. Physical feasibility and operating cost are evaluated separately using the classical dispatch model.

### Experiment B: Multi-Period Classical Scheduling

The second experiment extends the model to:

- 20 thermal generators
- Seven renewable generators
- 24 consecutive periods

We construct a **mixed-integer linear program (MILP)** using SciPy's HiGHS optimiser.

It includes binary commitment, startup and shutdown decisions, continuous thermal and renewable dispatch, piecewise-linear generation costs, reserve requirements, ramping constraints, and minimum operating durations.

The resulting schedules and plots are saved in `results/multi_period/`.

![Multi-period electricity demand and dispatch](results/multi_period/dispatch_plot.png)

### Experiment C: Time-Dependent QAOA

Directly simulating a quantum circuit for all 20 generators across 24 periods is computationally impractical.

Instead, we select four thermal generators over periods 16–18, while fixing the other 16 thermal generators to their classical MILP schedules.

The quantum subproblem therefore has:

\[
4\times3=12\text{ binary variables}
\]

and requires **12 qubits**.

The QUBO incorporates approximate operating costs, startup transitions between periods, and capacity-adequacy penalties. The penalty formulation is inspired by slack-free or *unbalanced penalisation* methods.

Every candidate quantum schedule is checked using a separate classical dispatch LP that enforces physical constraints within the selected time window.

All \(2^{12}=4096\) commitment patterns are enumerated to establish the exact classical reference solution for this conditional subproblem.

## 4. Qiskit Implementation

Both quantum experiments use Qiskit to construct parameterised QAOA circuits.

The implementation includes:

1. QUBO-to-Ising conversion using \(x_i=(1-Z_i)/2\).
2. Hadamard gates for initial state preparation.
3. `RZ` and `RZZ` gates for the cost Hamiltonian.
4. `RX` gates for the mixer.
5. Classical COBYLA optimisation of circuit parameters.
6. Statevector simulation for training and Qiskit Aer for finite-shot measurements.

Standard QAOA and CVaR-QAOA use the same circuit ansatz but different classical optimisation objectives.

The experiments use ideal quantum simulators, not physical quantum hardware.

## 5. Experimental Results

### 5.1 Single-Period Results

The exact classical optimum is commitment `011001`, with operating cost **4,278.46 benchmark cost units**.

Sampling results with 4,096 shots per method:

| Method | Optimal solution sampled | Feasible solutions |
|---|---:|---:|
| Standard QAOA | 6.10% | 64.50% |
| CVaR-QAOA | 14.75% | 69.36% |
| Uniform random | 1.71% | 67.09% |

All three methods sampled the classical optimum at least once, giving a best-found optimality gap of 0%.

![Six-qubit QAOA comparison](results/qaoa_comparison.png)

CVaR-QAOA concentrated more probability on the optimal commitment in this experiment. However, the original training runs had different optimiser budgets, and the later matched-settings runs still differed in actual function evaluations. The results should therefore not be interpreted as proof that CVaR-QAOA is universally superior.

### 5.2 Multi-Period Classical Results

| Metric | Result |
|---|---:|
| Thermal generators | 20 |
| Renewable generators | 7 |
| Scheduling periods | 24 |
| MILP variables | 3,528 |
| Binary variables | 1,440 |
| Constraints | 3,888 |
| Total operating cost | 166,249.36 |
| Renewable curtailment | 795.71 MW-periods |
| Reported MILP optimality gap | 0.0606% |

The solver produced a feasible multi-period schedule satisfying the implemented power-balance and reserve constraints. The reported MIP gap indicates termination within the solver tolerance, rather than a certified zero-gap optimum.

### 5.3 Twelve-Qubit Temporal QAOA Results

For the conditional three-period scheduling experiment:

| Metric | Result |
|---|---:|
| Qubits | 12 |
| QAOA layers | 1 |
| Candidate commitment schedules | 4,096 |
| Physically feasible schedules | 10 |
| Exact classical subproblem optimum | 40,103.16 |
| Physical cost of the QUBO optimum | 41,342.09 |
| QUBO approximation gap | 3.09% |

The quantum algorithms were evaluated with 4,096 measurement shots each.

| Metric | Standard QAOA | CVaR-QAOA |
|---|---:|---:|
| Physically feasible samples | 9.69% | 13.01% |
| Physically optimal samples | 1.07% | 2.44% |
| Best sampled physical cost | 40,103.16 | 40,103.16 |
| Function evaluations | 148 | 175 |
| Circuit depth | 15 | 15 |

Only 10 of the 4,096 candidate schedules are physically feasible, corresponding to a uniform-random feasibility probability of approximately **0.244%**.

Both QAOA variants sampled the physically optimal schedule, with CVaR-QAOA assigning greater observed sampling frequency to it in this run.

The 12-qubit experiment illustrates an important distinction between minimising approximate QUBO energy and minimising actual physical operating cost. The exact QUBO optimum had a 3.09% physical-cost gap, while the trained CVaR-QAOA circuit most frequently sampled the physically optimal configuration.

The conditional subproblem's cost is not directly comparable to the total 24-period MILP objective because the subproblem fixes the remaining thermal generators.

## 6. Repository Structure

```text
data/
    2020-01-27.json
    reduced_instance.json
    multi_period_instance.json

notebooks/
    01_data_exploration.ipynb
    02_multi_period_exploration.ipynb

src/
    multi_period_milp.py
    temporal_qaoa.py

results/
    classical_best.json
    classical_results.csv
    quantum_benchmark.csv
    controlled_comparison.csv
    measurement_counts.json
    qaoa_comparison.png
    cvar_distribution.png

    multi_period/
        period_summary.csv
        thermal_schedule.csv
        renewable_schedule.csv
        solver_summary.json
        dispatch_plot.png

    temporal_quantum/
        classical_enumeration.csv
        quantum_results.csv
        standard_counts.json
        cvar_counts.json
        summary.json
```

## 7. Installation and Reproduction

### Requirements

- Python
- NumPy
- Pandas
- SciPy
- Matplotlib
- Qiskit
- Qiskit Aer
- Jupyter Notebook or VS Code with Jupyter support

Install the project's dependencies using:

```bash
python -m pip install -r requirements.txt
```

### Running the notebooks

Open the `notebooks/` directory in VS Code or Jupyter.

**Notebook 1:** `01_data_exploration.ipynb`

Runs the single-period experiment, including the classical benchmark, QUBO construction, QAOA optimisation, finite-shot evaluation, and result plots.

**Notebook 2:** `02_multi_period_exploration.ipynb`

Explores the full benchmark, constructs the 24-period instance, and executes the larger classical and quantum experiments.

### Running the Python implementations

From the repository root, execute:

```bash
python src/multi_period_milp.py
python src/temporal_qaoa.py
```

Alternatively, from the second Jupyter notebook, use:

```python
%run ../src/multi_period_milp.py
%run ../src/temporal_qaoa.py
```

Run the classical MILP first because the temporal quantum experiment uses the saved thermal schedule as input.

Results are automatically exported to the corresponding directories under `results/`.

## 8. Limitations and Future Work

The current work is a proof of concept with several limitations:

- **No demonstrated quantum speedup:** Classical optimisation and enumeration remain computationally efficient for the tested instances.
- **Approximate quantum objective:** The QUBO does not exactly represent all physical dispatch costs and constraints.
- **Conditional quantum scheduling:** The 12-qubit experiment optimises only four generators over three periods, with other thermal generators fixed.
- **Simplified grid representation:** The model uses system-wide power balance rather than network-constrained power flow.
- **Simplified startup costs:** The multi-period model uses the cheapest available startup-cost tier instead of downtime-dependent startup categories.
- **No renewable uncertainty modelling:** Renewable availability is treated as known.
- **No hardware noise evaluation:** Quantum results come from ideal statevector and shot-based simulation.

Future extensions could investigate larger quantum subproblems, improved constraint-preserving formulations, noise-aware QAOA, more realistic startup-cost modelling, rolling-horizon scheduling, and integrated classical–quantum decomposition.

## 9. Research References

1. **PGLib-UC Benchmark Library**, IEEE PES Task Force on Benchmarks for Validation of Emerging Power System Algorithms. [GitHub repository](https://github.com/power-grid-lib/pglib-uc).

2. **Koretsky et al. (2021)**, *Adapting Quantum Approximation Optimization Algorithm (QAOA) for Unit Commitment*, IEEE International Conference on Quantum Computing and Engineering. [DOI: 10.1109/QCE52317.2021.00035](https://doi.org/10.1109/QCE52317.2021.00035).

3. **Moncayo-Martínez and He (2026)**, *Quantum optimisation for supply chain: QUBO formulations and QAOA solutions for facility location and load balancing*, Results in Engineering. [DOI: 10.1016/j.rineng.2025.108373](https://doi.org/10.1016/j.rineng.2025.108373).

4. **Barkoutsos et al. (2020)**, *Improving Variational Quantum Optimization using CVaR*, Quantum 4, 256. [DOI: 10.22331/q-2020-04-20-256](https://doi.org/10.22331/q-2020-04-20-256).

---

**Project summary:** This project demonstrates a reproducible hybrid quantum–classical approach to renewable-aware unit commitment. It combines classical MILP scheduling, QUBO-based quantum optimisation, and physical feasibility verification, while explicitly examining the limitations of approximate quantum formulations and the absence of computational quantum advantage in the tested instances.