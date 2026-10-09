### 1. Novelty

We developed a hybrid quantum–classical approach for renewable-aware unit commitment. A 24-period classical scheduling model is combined with a 12-qubit QAOA subproblem spanning three periods. Our approach uses slack-free capacity penalties and separate physical feasibility checks, allowing us to investigate differences between approximate QUBO optimisation and actual grid operating costs.

### 2. Level of Qiskit Programming

We implemented QAOA directly using Qiskit, including QUBO-to-Ising conversion, parameterised RZ/RZZ cost layers, RX mixers, and COBYLA optimisation. We also implemented CVaR-QAOA, used statevector simulation for parameter optimisation, and performed shot-based measurements using Qiskit Aer. Our experiments use six-qubit and twelve-qubit circuits.

### 3. Measurable Results / Classical Benchmarking

The 24-period classical MILP achieved an operating cost of 166,249.36. For the 12-qubit subproblem, exact classical enumeration identified 10 feasible schedules among 4,096 possibilities. Standard QAOA and CVaR-QAOA sampled feasible schedules at 9.69% and 13.01%, respectively, and both recovered the classical optimum of 40,103.16.

### 4. Technical Quantum Advantage

Our slack-free QUBO formulation avoids additional qubits for capacity-adequacy slack variables. CVaR-QAOA achieved a 13.01% feasible-sampling rate, compared with 0.244% for uniform random sampling, demonstrating probability concentration on feasible schedules. However, we do not claim computational quantum speedup, since classical optimisation remains efficient for these instances.