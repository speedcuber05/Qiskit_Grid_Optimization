"""12-qubit, 3-period conditional unit commitment experiment (RTS-GMLC).

Usage from the project root:
    python src/temporal_qaoa.py
    python src/temporal_qaoa.py --classical-only

Required: scipy, numpy, pandas, matplotlib, qiskit, qiskit-aer.

This is a *conditional* subproblem: nonselected thermal generation is fixed
at the 24-period MILP solution over periods 16-18. Selected units and all
renewable dispatch are optimised. QAOA minimises a quadratic *surrogate*;
physical/economic feasibility is evaluated separately by an LP with ramps.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
from time import perf_counter
import numpy as np
import pandas as pd
from scipy.optimize import linprog, minimize

SELECTED = ["202_STEAM_3", "101_STEAM_4", "115_STEAM_3", "315_CT_7"]
PERIODS = (16, 17, 18)  # 1-indexed periods
CVaR_FRACTION = 0.2
SEED = 43


def locate_root():
    script_dir = Path(__file__).resolve().parent
    if (script_dir.parent / 'data' / 'multi_period_instance.json').exists():
        return script_dir.parent
    if (Path.cwd() / 'data' / 'multi_period_instance.json').exists():
        return Path.cwd()
    raise FileNotFoundError('Run the 24-period MILP first, from the project root')


def build_context(root):
    case = json.loads((root/'data'/'multi_period_instance.json').read_text())
    sched = pd.read_csv(root/'results'/'multi_period'/'thermal_schedule.csv')
    assert set(SELECTED) <= set(case['thermal_generators'])
    t_indices = [t-1 for t in PERIODS]
    G = len(SELECTED)
    L = len(PERIODS)
    units = [case['thermal_generators'][x] for x in SELECTED]
    previous = {}
    for name, unit in zip(SELECTED,units):
        prev = sched[(sched['period']==PERIODS[0]-1)&(sched['generator']==name)].iloc[0]
        state = int(prev['on'])
        streak = 0
        for tp in range(PERIODS[0]-1,0,-1):
            row = sched[(sched['period']==tp)&(sched['generator']==name)].iloc[0]
            if int(row['on']) != state: break
            streak += 1
        previous[name] = {'on':state,'mw':float(prev['power_mw']),'streak':streak}
    fixed_mw = []
    fixed_headroom = []
    renewable_max = []
    renewable_min = []
    for t in PERIODS:
        background = sched[(sched.period==t)&(~sched.generator.isin(SELECTED))]
        fixed_mw.append(float(background.power_mw.sum()))
        fixed_headroom.append(float(sum(
            case['thermal_generators'][r.generator]['power_output_maximum']*int(r.on)-r.power_mw
            for r in background.itertuples())))
        renewable_max.append(float(sum(v['power_output_maximum'][t-1] for v in case['renewable_generators'].values())))
        renewable_min.append(float(sum(v['power_output_minimum'][t-1] for v in case['renewable_generators'].values())))
    demand = np.array([case['demand'][j] for j in t_indices],float)
    reserve = np.array([case['reserves'][j] for j in t_indices],float)
    return dict(case=case,sched=sched,units=units,previous=previous,L=L,G=G,
                demand=demand,reserve=reserve,fixed=np.array(fixed_mw),
                headroom=np.array(fixed_headroom),
                renew_max=np.array(renewable_max),renew_min=np.array(renewable_min))


def bits_from_index(idx,G=4,L=3):
    # Generator-major: x_(i,t) occupies qubit i*L+t.
    return np.array([(idx>>j)&1 for j in range(G*L)],dtype=int).reshape(G,L)


def solve_physical(ctx, y):
    """LP dispatch for a commitment; infeasible returns None.
    Uses continuous piecewise segments, renewable aggregate, reserve and ramps;
    checks min up/down inside the 3-period horizon. Background dispatch fixed.
    """
    G,L=ctx['G'],ctx['L']
    U=ctx['units']
    starts=np.zeros((G,L),dtype=int); stops=np.zeros((G,L),dtype=int)
    for i, unit in enumerate(U):
        prev=ctx['previous'][SELECTED[i]]['on']
        m_up=int(unit['time_up_minimum']);m_down=int(unit['time_down_minimum'])
        prev_streak=ctx['previous'][SELECTED[i]]['streak']
        # Respect minimum remaining time from before block (if any).
        if prev and prev_streak < m_up and not np.all(y[i,:min(L,m_up-prev_streak)]==1):return None
        if not prev and prev_streak < m_down and not np.all(y[i,:min(L,m_down-prev_streak)]==0):return None
        for t in range(L):
            old=prev if t==0 else int(y[i,t-1]); new=int(y[i,t])
            starts[i,t]=int(old==0 and new==1)
            stops[i,t]=int(old==1 and new==0)
            if starts[i,t] and any(y[i,t:min(L,t+m_up)]==0):return None
            if stops[i,t] and any(y[i,t:min(L,t+m_down)]==1):return None
    # Variables: incremental cost-curve segments (one set per generator-period)
    # plus one total renewable-use variable per period.
    segment_indices={};bounds=[];cost=[]
    pmin=np.array([float(g['power_output_minimum']) for g in U]);
    pmax=np.array([float(g['power_output_maximum']) for g in U]);
    production_base=0.0
    startup_cost=0.0
    for i,g in enumerate(U):
        points=g['piecewise_production']
        startup_price=min(float(s['cost']) for s in g['startup'])
        for t in range(L):
            production_base+=points[0]['cost']*y[i,t]
            startup_cost+=startup_price*starts[i,t]
            segment_indices[i,t]=[]
            for a,b in zip(points[:-1],points[1:]):
                j=len(cost); segment_indices[i,t].append(j)
                bounds.append((0.0,(b['mw']-a['mw'])*y[i,t]))
                cost.append((b['cost']-a['cost'])/(b['mw']-a['mw']))
    renew_indices=[]
    for t in range(L):
        renew_indices.append(len(cost));cost.append(0.0)
        bounds.append((ctx['renew_min'][t],ctx['renew_max'][t]))
    N=len(cost);eq_rows=[];eq_values=[];ub_rows=[];ub_values=[]
    def expr(i,t,sign=1):
        v=np.zeros(N)
        for j in segment_indices[i,t]: v[j]+=sign
        return v
    for t in range(L):
        # Σ thermal P + renewable generation = residual load.
        row=np.zeros(N)
        row[renew_indices[t]]=1.0
        for i in range(G):row+=expr(i,t)
        eq_rows.append(row)
        eq_values.append(ctx['demand'][t]-ctx['fixed'][t]-float(pmin@y[:,t]))
        # thermal unused reserve; background includes its original headroom
        row=np.zeros(N)
        for i in range(G):row+=expr(i,t)
        ub_rows.append(row)
        ub_values.append(float((pmax-pmin)@y[:,t]+ctx['headroom'][t]-ctx['reserve'][t]))
        for i,g in enumerate(U):
            if t==0:
                old=ctx['previous'][SELECTED[i]]['on']
                old_power=ctx['previous'][SELECTED[i]]['mw']
                old_expr=np.zeros(N)
            else:
                old=int(y[i,t-1]);old_power=pmin[i]*old
                old_expr=expr(i,t-1)
            RU=float(g['ramp_up_limit']);RD=float(g['ramp_down_limit'])
            SU=float(g['ramp_startup_limit']);SD=float(g['ramp_shutdown_limit'])
            ub_rows.append(expr(i,t)-old_expr)
            ub_values.append(old_power + RU*old + SU*starts[i,t]-pmin[i]*y[i,t])
            ub_rows.append(old_expr-expr(i,t))
            ub_values.append(RD*y[i,t]+SD*stops[i,t]+pmin[i]*y[i,t]-old_power)
    sol=linprog(np.asarray(cost),A_ub=ub_rows,b_ub=ub_values,
                A_eq=eq_rows,b_eq=eq_values,bounds=bounds,method='highs')
    if not sol.success:return None
    return dict(cost=float(sol.fun+production_base+startup_cost),
                startup_cost=startup_cost, production_cost=float(sol.fun+production_base),
                y=y.copy())


def build_qubo(ctx,alpha=160.0,beta=0.3):
    """Construct QUBO from run costs, transition startup and adequacy surrogate.

    QUBO is an *approximation* of physical dispatch. It does not encode ramp
    rates, min-up/down, dispatch cost curvature or exact feasibility.
    """
    G,L=ctx['G'],ctx['L'];n=G*L
    linear=np.zeros(n);quad=np.zeros((n,n));constant=0.
    pmax=np.array([g['power_output_maximum'] for g in ctx['units']])
    fixed_min=np.array([g['piecewise_production'][0]['cost'] for g in ctx['units']])
    startups=np.array([min(s['cost'] for s in g['startup']) for g in ctx['units']])
    for i in range(G):
        for t in range(L):
            k=i*L+t
            linear[k]+=fixed_min[i]+startups[i]
            if t==0:
                linear[k]-=startups[i]*ctx['previous'][SELECTED[i]]['on']
            else:
                q=i*L+t-1
                quad[q,k]-=startups[i]
    # Need selected thermal capacity after max renewable + fixed thermal,
    # accounting for background thermal reserve headroom.
    for t in range(L):
        target=(ctx['demand'][t]-ctx['fixed'][t]-ctx['renew_max'][t]
                +max(0.,ctx['reserve'][t]-ctx['headroom'][t]))
        # alpha*(target - ΣPmax*y)+beta*(target - ΣPmax*y)^2
        constant+=alpha*target+beta*target**2
        for i in range(G):
            k=i*L+t;cap=pmax[i]
            linear[k]+=-alpha*cap+beta*(cap**2-2*target*cap)
            for j in range(i+1,G):
                m=j*L+t
                quad[k,m]+=2*beta*cap*pmax[j]
    return constant,linear,quad


def all_qubo_energies(c,linear,quad):
    n=len(linear)
    # basis labels: bit i corresponds to qubit i, same convention throughout
    states=((np.arange(2**n)[:,None]>>np.arange(n))&1).astype(float)
    E=c+states@linear+np.einsum('ni,ij,nj->n',states,quad,states)
    return states,E


def run_quantum(linear,quad,E,exact_feasible,physical_cost,out,depth=1,shots=4096):
    from qiskit import QuantumCircuit,transpile
    from qiskit.quantum_info import Statevector
    from qiskit_aer import AerSimulator
    n=len(linear)
    # x_i=(1-Z_i)/2; E=E0+Σh_i Z_i+ΣJ_ij Z_i Z_j
    J=quad/4.0
    h=-linear/2.0-(quad.sum(axis=0)+quad.sum(axis=1))/4.0
    # Verify Pauli mapping, ignoring the physically irrelevant constant shift.
    z=1-2*((np.arange(2**n)[:,None]>>np.arange(n))&1)
    ising=z@h+np.einsum('ni,ij,nj->n',z,J,z)
    assert np.allclose(E-E[0],ising-ising[0],atol=1e-5)
    scale=max(1000.,np.std(E))
    rng=np.random.default_rng(SEED)
    order=np.argsort(E)
    def circuit(params,meas=False):
        qc=QuantumCircuit(n)
        qc.h(range(n))
        for layer in range(depth):
            gamma=params[layer];b=params[depth+layer]
            for i in range(n):
                if abs(h[i])>1e-12:qc.rz(2*gamma*h[i]/scale,i)
            for i in range(n):
                for j in range(i+1,n):
                    if abs(J[i,j])>1e-12:qc.rzz(2*gamma*J[i,j]/scale,i,j)
            for i in range(n):qc.rx(2*b,i)
        if meas:qc.measure_all()
        return qc
    def probs(params):
        sv=Statevector.from_instruction(circuit(params))
        return np.abs(sv.data)**2
    def obj(params,fract):
        p=probs(params)
        if fract==1.:return float(p@E)/scale
        cum=0.;acc=0.
        for idx in order:
            v=min(fract-cum,p[idx]);acc+=v*E[idx];cum+=v
            if cum>=fract-1e-12:break
        return acc/fract/scale
    initials=[np.r_[rng.uniform(0,2*np.pi,depth),rng.uniform(0,np.pi/2,depth)] for _ in range(3)]
    simulator=AerSimulator(method='statevector')
    results=[]
    for name,frac in [('Standard QAOA',1.),('CVaR-QAOA',CVaR_FRACTION)]:
        start=perf_counter();fits=[]
        for x0 in initials:
            fits.append(minimize(obj,x0,args=(frac,),method='COBYLA',
                                 options={'maxiter':130,'rhobeg':0.4,'tol':1e-3}))
        fit=min(fits,key=lambda x:x.fun)
        duration=perf_counter()-start
        qc=circuit(fit.x,meas=True)
        backend_circuit=transpile(qc,simulator)
        # qiskit reports reversed bit order; integer index avoids confusion.
        counts=simulator.run(backend_circuit,shots=shots,seed_simulator=SEED).result().get_counts()
        count_by_idx={int(s.replace(' ',''),2):cnt for s,cnt in counts.items()}
        optimal_idx=int(np.argmin(E))
        p=probs(fit.x)
        result={
            'method':name,'depth':depth,'qubits':n,
            'optimization_seconds':duration,'function_evaluations':sum(x.nfev for x in fits),
            'probability_exact_qubo_optimum':float(p[optimal_idx]),
            'probability_physical_optimum':float(p[np.nanargmin(physical_cost)]),
            'sampled_exact_qubo_optimum_pct':100*count_by_idx.get(optimal_idx,0)/shots,
            'physically_feasible_pct':100*sum(c for idx,c in count_by_idx.items() if exact_feasible[idx])/shots,
            'physical_optimum_sampled_pct':100*count_by_idx.get(int(np.nanargmin(physical_cost)),0)/shots,
            'best_sampled_physical_cost':min((physical_cost[idx] for idx in count_by_idx if exact_feasible[idx]),default=None),
            'most_probable_bitstring':format(max(count_by_idx,key=count_by_idx.get),'012b')[::-1],
            'circuit_depth':qc.depth(),'rzz_count':qc.count_ops().get('rzz',0),
        }
        results.append(result)
        print(name, '\n', json.dumps(
    result,
    indent=2,
    default=lambda x: x.item() if isinstance(x, np.generic) else str(x)
))
        (out/('standard_counts.json' if frac==1. else 'cvar_counts.json')).write_text(json.dumps(count_by_idx,indent=2))
    return results


def main(classical_only=False):
    root=locate_root();out=root/'results'/'temporal_quantum';out.mkdir(parents=True,exist_ok=True)
    ctx=build_context(root);G,L=ctx['G'],ctx['L'];n=G*L
    print('Conditional quantum problem: 4 thermal units x 3 periods = 12 bits')
    print('Selected thermal:',SELECTED,'; periods:',PERIODS)
    print(pd.DataFrame({'period':PERIODS,'demand_mw':ctx['demand'],
          'renewable_available_mw':ctx['renew_max'],
          'fixed_other_thermal_mw':ctx['fixed'],
          'reserve_mw':ctx['reserve']}))
    c,linear,quad=build_qubo(ctx)
    states,E=all_qubo_energies(c,linear,quad)
    qubo_idx=int(np.argmin(E))
    print('QUBO optimum:',format(qubo_idx,f'0{n}b')[::-1],f'energy={E[qubo_idx]:,.1f}')
    # Enumerate exact physical dispatch LPs of 4096 binary schedules.
    feasible=np.zeros(2**n,dtype=bool);physical_cost=np.full(2**n,np.nan)
    tic=perf_counter()
    for idx in range(2**n):
        ans=solve_physical(ctx,bits_from_index(idx,G,L))
        if ans is not None:
            feasible[idx]=True
            physical_cost[idx]=ans['cost']
    assert feasible.any(), 'No physically feasible conditional schedules!'
    valid_idx=np.flatnonzero(feasible)
    classical_best_idx=int(valid_idx[np.argmin(physical_cost[valid_idx])])
    print('LP enumeration time:',round(perf_counter()-tic,2),'s')
    print('Physically feasible schedules:',len(valid_idx),'/',2**n)
    print('Exact classical optimum:',format(classical_best_idx,f'0{n}b')[::-1],
          'cost:',round(physical_cost[classical_best_idx],2))
    print('QUBO optimum physically feasible:',bool(feasible[qubo_idx]))
    if feasible[qubo_idx]:
        gap=100*(physical_cost[qubo_idx]/physical_cost[classical_best_idx]-1)
        print(f'QUBO optimum physical cost: {physical_cost[qubo_idx]:.2f}; gap: {gap:.2f}%')
    table=pd.DataFrame({'bitmask':np.arange(2**n),
       'bitstring_qubit_order':[format(i,f'0{n}b')[::-1] for i in range(2**n)],
       'qubo_energy':E,'physical_feasible':feasible,'physical_cost':physical_cost})
    table.to_csv(out/'classical_enumeration.csv',index=False)
    metadata={'selected_generators':SELECTED,'periods':PERIODS,
      'qubits':n,'qubo_optimum_bitstring':format(qubo_idx,f'0{n}b')[::-1],
      'qubo_optimum_energy':float(E[qubo_idx]),
      'physical_optimum_bitstring':format(classical_best_idx,f'0{n}b')[::-1],
      'physical_optimum_cost':float(physical_cost[classical_best_idx]),
      'feasible_schedules':int(feasible.sum()),'total_schedules':int(2**n),
      'conditioning':'16 nonselected thermal generators frozen to 24-period MILP dispatch',
      'limitations':'QUBO is a surrogate, omits physical dispatch/ramp/min-up/min-down. Physical check uses exact LP with frozen thermal background; no quantum advantage claim.'}
    (out/'summary.json').write_text(json.dumps(metadata,indent=2))
    if not classical_only:
        quantum_results=run_quantum(linear,quad,E,feasible,physical_cost,out)
        pd.DataFrame(quantum_results).to_csv(out/'quantum_results.csv',index=False)
    print('Saved:',out)

if __name__=='__main__':
    ap=argparse.ArgumentParser()
    ap.add_argument('--classical-only',action='store_true')
    args=ap.parse_args()
    main(args.classical_only)
