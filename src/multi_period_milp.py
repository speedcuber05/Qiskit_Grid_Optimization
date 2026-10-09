"""24-period renewable-aware unit commitment baseline using SciPy/HiGHS MILP.

Run from repository root: python src/multi_period_milp.py
Or from notebooks/: %run ../src/multi_period_milp.py

Input: data/multi_period_instance.json (built by 02_multi_period_exploration.ipynb)
Outputs: results/multi_period/*

Simplifications: single-bus power balance; thermal spinning reserve; cheapest startup
cost tier; no emissions, grid network, renewables costs, or terminal look-ahead.
"""

import json
from pathlib import Path
from time import perf_counter

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import Bounds, LinearConstraint, milp
from scipy.sparse import coo_matrix


def solve_project(project_root: Path, time_limit: float = 90, rel_gap: float = 0.001):
    source = project_root / "data" / "multi_period_instance.json"
    if not source.is_file():
        raise FileNotFoundError(f"Missing {source}. First run the dataset construction cells.")
    with source.open(encoding="utf-8") as file:
        case = json.load(file)

    T = int(case["time_periods"])
    thermal = case["thermal_generators"]
    renewable = case["renewable_generators"]
    G = list(thermal)
    W = list(renewable)
    D = np.asarray(case["demand"], dtype=float)
    R = np.asarray(case["reserves"], dtype=float)
    assert len(D) == len(R) == T

    names = []
    lower, upper, c, integer = [], [], [], []

    def var(name, lb=0.0, ub=np.inf, cost=0.0, binary=False):
        j = len(names)
        names.append(name)
        lower.append(lb)
        upper.append(ub)
        c.append(cost)
        integer.append(1 if binary else 0)
        return j

    y, u, v, p, segments, w = {}, {}, {}, {}, {}, {}

    for i in G:
        g = thermal[i]
        pmin, pmax = float(g["power_output_minimum"]), float(g["power_output_maximum"])
        points = g["piecewise_production"]
        assert len(points) >= 2
        assert abs(points[0]["mw"] - pmin) < 1e-5
        assert abs(points[-1]["mw"] - pmax) < 1e-5
        slopes = np.array([
            (b["cost"] - a["cost"]) / (b["mw"] - a["mw"])
            for a, b in zip(points[:-1], points[1:])
        ])
        if np.any(np.diff(slopes) < -1e-6):
            raise ValueError(f"Non-convex piecewise cost: {i}; requires segment-order binaries")

        # Single cheapest startup tier for every start (simplification).
        startup_price = min(float(x["cost"]) for x in g["startup"])
        for t in range(T):
            y[i, t] = var(f"on/{i}/{t}", ub=1, cost=float(points[0]["cost"]), binary=True)
            u[i, t] = var(f"start/{i}/{t}", ub=1, cost=startup_price, binary=True)
            v[i, t] = var(f"stop/{i}/{t}", ub=1, binary=True)
            p[i, t] = var(f"power/{i}/{t}", ub=pmax)
            for k, (a, b) in enumerate(zip(points[:-1], points[1:])):
                segments[i, t, k] = var(
                    f"segment/{i}/{t}/{k}",
                    ub=float(b["mw"] - a["mw"]),
                    cost=float(slopes[k]),
                )

    for j in W:
        g = renewable[j]
        for t in range(T):
            lo, hi = float(g["power_output_minimum"][t]), float(g["power_output_maximum"][t])
            w[j, t] = var(f"renewable/{j}/{t}", lb=lo, ub=hi)

    # Construct sparse mixed-integer linear constraints, without dense matrices.
    ridx, cidx, values, lhs, rhs = [], [], [], [], []

    def add(row, lo=-np.inf, hi=np.inf):
        r = len(lhs)
        for j, coeff in row.items():
            if abs(coeff) > 1e-14:
                ridx.append(r)
                cidx.append(j)
                values.append(float(coeff))
        lhs.append(float(lo))
        rhs.append(float(hi))

    for t in range(T):
        # Nodal power balance: sum thermal + sum renewable = load.
        add({**{p[i, t]: 1 for i in G}, **{w[j, t]: 1 for j in W}}, D[t], D[t])

        # Spinning reserve is unused capacity of committed thermal generators.
        add({**{y[i, t]: float(thermal[i]["power_output_maximum"]) for i in G},
             **{p[i, t]: -1 for i in G}}, lo=R[t])

        for i in G:
            g = thermal[i]
            pmin = float(g["power_output_minimum"])
            pmax = float(g["power_output_maximum"])
            ups = int(g["time_up_minimum"])
            downs = int(g["time_down_minimum"])

            if g.get("must_run", 0):
                add({y[i, t]: 1}, 1, 1)

            # Power = minimum operating output + segment increments, if ON.
            eq = {p[i, t]: 1, y[i, t]: -pmin}
            for k in range(len(g["piecewise_production"]) - 1):
                eq[segments[i, t, k]] = -1
            add(eq, 0, 0)
            add({p[i, t]: 1, y[i, t]: -pmax}, hi=0)

            # Commitment transition: on_t - on_(t-1) = start_t - stop_t.
            if t == 0:
                prev_on = int(g["unit_on_t0"])
                prev_power = float(g["power_output_t0"])
                add({y[i, t]: 1, u[i, t]: -1, v[i, t]: 1}, prev_on, prev_on)
            else:
                add({y[i, t]: 1, y[i, t - 1]: -1,
                     u[i, t]: -1, v[i, t]: 1}, 0, 0)
            add({u[i, t]: 1, v[i, t]: 1}, hi=1)

            # Ramping with explicit startup/shutdown allowances.
            RU = float(g["ramp_up_limit"])
            RD = float(g["ramp_down_limit"])
            SU = float(g["ramp_startup_limit"])
            SD = float(g["ramp_shutdown_limit"])
            if t == 0:
                add({p[i, t]: 1, u[i, t]: -SU}, hi=prev_power + RU * prev_on)
                add({p[i, t]: -1, y[i, t]: -RD, v[i, t]: -SD}, hi=-prev_power)
            else:
                add({p[i, t]: 1, p[i, t - 1]: -1,
                     y[i, t - 1]: -RU, u[i, t]: -SU}, hi=0)
                add({p[i, t - 1]: 1, p[i, t]: -1,
                     y[i, t]: -RD, v[i, t]: -SD}, hi=0)

            # Minimum up/down time constraints.
            add({**{u[i, k]: 1 for k in range(max(0, t - ups + 1), t + 1)},
                 y[i, t]: -1}, hi=0)
            add({**{v[i, k]: 1 for k in range(max(0, t - downs + 1), t + 1)},
                 y[i, t]: 1}, hi=1)

            # Units that haven't yet satisfied initial min up/down duration.
            if int(g["unit_on_t0"]) and t < max(0, ups - int(g["time_up_t0"])):
                add({y[i, t]: 1}, 1, 1)
            if not int(g["unit_on_t0"]) and t < max(0, downs - int(g["time_down_t0"])):
                add({y[i, t]: 1}, 0, 0)

    A = coo_matrix((values, (ridx, cidx)), shape=(len(lhs), len(names))).tocsr()
    constraints = LinearConstraint(A, np.asarray(lhs), np.asarray(rhs))
    print(f"Model: {len(G)} thermal, {len(W)} renewable, {T} periods")
    print(f"MILP variables: {len(names)}; binary: {sum(integer)}; constraints: {len(lhs)}")
    print("Solving with HiGHS...")
    tic = perf_counter()
    result = milp(
        c=np.asarray(c),
        integrality=np.asarray(integer),
        bounds=Bounds(lower, upper),
        constraints=constraints,
        options={"time_limit": time_limit, "mip_rel_gap": rel_gap, "disp": False},
    )
    runtime = perf_counter() - tic
    print(f"Solver status: {result.message}")
    if result.x is None:
        raise RuntimeError("No feasible schedule found; inspect model or increase solver time limit")

    x = result.x
    records = []
    thermal_rows = []
    wind_rows = []
    for t in range(T):
        tp = sum(x[p[i, t]] for i in G)
        rp = sum(x[w[j, t]] for j in W)
        ravail = sum(renewable[j]["power_output_maximum"][t] for j in W)
        headroom = sum(thermal[i]["power_output_maximum"] * round(x[y[i, t]])
                       - x[p[i, t]] for i in G)
        records.append({
            "period": t + 1, "demand_mw": D[t], "thermal_mw": tp,
            "renewable_used_mw": rp, "renewable_available_mw": ravail,
            "curtailed_mw": max(0, ravail - rp),
            "reserve_required_mw": R[t], "reserve_headroom_mw": headroom,
            "generators_on": sum(round(x[y[i, t]]) for i in G),
            "starts": sum(round(x[u[i, t]]) for i in G),
            "stops": sum(round(x[v[i, t]]) for i in G),
        })
        for i in G:
            thermal_rows.append({"period": t + 1, "generator": i,
                                 "on": round(x[y[i, t]]), "startup": round(x[u[i, t]]),
                                 "shutdown": round(x[v[i, t]]), "power_mw": x[p[i, t]]})
        for j in W:
            wind_rows.append({"period": t + 1, "generator": j, "power_mw": x[w[j, t]],
                              "available_mw": renewable[j]["power_output_maximum"][t]})

    profile = pd.DataFrame(records)
    commitments = pd.DataFrame(thermal_rows)
    renewables = pd.DataFrame(wind_rows)
    assert np.allclose(profile["thermal_mw"] + profile["renewable_used_mw"], D, atol=1e-4)
    assert (profile["reserve_headroom_mw"] + 1e-4 >= R).all()

    out = project_root / "results" / "multi_period"
    out.mkdir(parents=True, exist_ok=True)
    profile.to_csv(out / "period_summary.csv", index=False)
    commitments.to_csv(out / "thermal_schedule.csv", index=False)
    renewables.to_csv(out / "renewable_schedule.csv", index=False)
    metadata = {
        "status_code": int(result.status), "status": str(result.message),
        "objective": float(result.fun), "runtime_seconds": runtime,
        "mip_gap": float(result.mip_gap) if result.mip_gap is not None else None,
        "time_periods": T, "thermal_count": len(G), "renewable_count": len(W),
        "total_curtailed_mwh_if_1hour_periods": float(profile["curtailed_mw"].sum()),
        "startup_cost_model": "cheapest listed startup tier for every event",
    }
    with (out / "solver_summary.json").open("w", encoding="utf-8") as file:
        json.dump(metadata, file, indent=2)

    fig, ax = plt.subplots(figsize=(11, 4.6))
    ax.plot(profile.period, profile.demand_mw, label="Demand", linewidth=2)
    ax.plot(profile.period, profile.thermal_mw, label="Thermal dispatch")
    ax.plot(profile.period, profile.renewable_used_mw, label="Renewables dispatched")
    ax.plot(profile.period, profile.renewable_available_mw, label="Renewables available", linestyle="--")
    ax.set(xlabel="Period", ylabel="Power (MW)", title="24-period unit commitment and dispatch")
    ax.grid(alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(out / "dispatch_plot.png", dpi=170)
    plt.close(fig)

    print(f"Objective: {result.fun:,.2f} benchmark cost units")
    print(f"Solve time: {runtime:.2f}s; reported MIP gap: {metadata['mip_gap']}")
    print(f"Total curtailed power summed over periods: {profile.curtailed_mw.sum():,.2f} MW-periods")
    print(f"Saved results to {out}")
    print(profile[["period", "demand_mw", "thermal_mw", "renewable_used_mw",
                   "curtailed_mw", "generators_on"]].head(8).round(2).to_string(index=False))
    return result, profile, commitments, renewables


if __name__ == "__main__":
    root = Path(__file__).resolve().parent.parent
    solve_project(root)
