#!/usr/bin/env python3
import os
import sys
import time
import shutil
import numpy as np
import pandas as pd
import plotly.graph_objs as go
import plotly.io as pio

# Set up paths
currentdir = os.path.dirname(os.path.abspath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir)

import pyLLE

# Configuration
JULIA_BIN = "julia"  # Assumes julia is in your PATH.
# If julia is at a specific location, set it here, e.g.:
# JULIA_BIN = "/Applications/Julia-1.8.app/Contents/Resources/julia/bin/julia"

def check_julia():
    if shutil.which(JULIA_BIN) is None:
        print(f"WARNING: '{JULIA_BIN}' not found in PATH. Simulation may fail.")
        print("Please ensure Julia is installed and added to PATH, or update JULIA_BIN in this script.")

def run_benchmark():
    check_julia()
    
    # -------------------- Simulation 1 --------------------
    print(f"\n{'='*20} Running Simulation 1 {'='*20}")
    
    res = dict(
        R=23e-6, 
        Qi=1e6, 
        Qc=1e6, 
        γ=3.2, 
        dispfile=os.path.join(currentdir, "RW1000_H430.csv")
    )

    sim = dict(
        Pin=[140e-3], 
        f_pmp=[283e12],
        φ_pmp=[0], 
        δω=[None], 
        Tscan=0.7e6,
        μ_sim=[-220, 220],
        δω_init= 1e9 * 2 * np.pi,
        δω_end= -6.5e9 * 2 * np.pi,
        num_probe = 5000, 
    )

    solver = pyLLE.LLEsolver(sim=sim, res=res, debug=False)

    # Analyze (Plotting skipped or non-blocking)
    # fig = solver.Analyze(plot=True)
    # fig.update_yaxes(range = [-50, 50], title = 'D<sub>int</sub> (GHz)')
    # fig.show()

    solver.Setup(verbose=False)
    
    print("Starting temporal solver...")
    start_time = time.time()
    try:
        solver.SolveTemporal(bin=JULIA_BIN)
        solver.RetrieveData()
        end_time = time.time()
        print(f"Simulation 1 completed in {end_time - start_time:.2f} seconds")
    except Exception as e:
        print(f"Simulation 1 failed: {e}")

    # Plot results (optional)
    # fig = go.Figure()
    # tr = go.Scatter(y=solver.sol.Pcomb * 1e3)
    # fig.add_trace(tr)
    # fig.show()

    # -------------------- Simulation 2 --------------------
    print(f"\n{'='*20} Running Simulation 2 {'='*20}")

    dks_init_path = os.path.join(currentdir, 'DKS_init.csv')
    if os.path.exists(dks_init_path):
        df = pd.read_csv(dks_init_path)
        df["DKS_init"] = df.DKS_init.apply(lambda val: complex(val.strip('()')))
        δω = df.det.values[0] * 2*np.pi
        DKS_init = df.DKS_init.values
        
        sim = dict(
            Pin=[160e-3], 
            f_pmp=[283e12], 
            φ_pmp=[0], 
            δω=[None], 
            Tscan=0.2e6,
            μ_sim=[-220, 220],
            δω_end = δω, δω_init = δω, 
            DKS_init =  DKS_init, 
        )
        solver = pyLLE.LLEsolver(sim=sim, res=res, debug=False) # reusing res from sim 1

        solver.Setup(verbose=False)
        print("Starting temporal solver...")
        start_time = time.time()
        try:
            solver.SolveTemporal(bin=JULIA_BIN)
            solver.RetrieveData()
            end_time = time.time()
            print(f"Simulation 2 completed in {end_time - start_time:.2f} seconds")
        except Exception as e:
             print(f"Simulation 2 failed: {e}")
    else:
        print(f"Skipping Simulation 2: {dks_init_path} not found")


    # -------------------- Simulation 3 --------------------
    print(f"\n{'='*20} Running Simulation 3 {'='*20}")
    
    # Re-reading csv just to be safe as per notebook flow
    if os.path.exists(dks_init_path):
        df = pd.read_csv(dks_init_path)
        df["DKS_init"] = df.DKS_init.apply(lambda val: complex(val.strip('()')))
        δω = df.det.values[0] * 2*np.pi
        DKS_init = df.DKS_init.values
        
        sim = dict(
            Pin=[160e-3], 
            f_pmp=[283e12], 
            φ_pmp=[0], 
            δω=[None], 
            Tscan=0.2e6,
            μ_sim=[-220, 220],
            δω_end = δω, δω_init = δω, 
            DKS_init =  DKS_init, 
            D1_manual = 982.2557817183881 * 1e9 * 2*np.pi 
        )
        solver = pyLLE.LLEsolver(sim=sim, res=res, debug=False)

        solver.Setup(verbose=False)
        print("Starting temporal solver...")
        start_time = time.time()
        try:
            solver.SolveTemporal(bin=JULIA_BIN)
            solver.RetrieveData()
            end_time = time.time()
            print(f"Simulation 3 completed in {end_time - start_time:.2f} seconds")
        except Exception as e:
             print(f"Simulation 3 failed: {e}")

    # -------------------- Simulation 4 --------------------
    print(f"\n{'='*20} Running Simulation 4 {'='*20}")

    if os.path.exists(dks_init_path):
        df = pd.read_csv(dks_init_path)
        df["DKS_init"] = df.DKS_init.apply(lambda val: complex(val.strip('()')))
        δω = df.det.values[0] * 2*np.pi
        DKS_init = df.DKS_init.values

        sim = dict(
            Pin=[160e-3, 10e-3], 
            f_pmp=[283e12, 240e12], 
            φ_pmp=[0, 0], 
            δω=[None, 0], 
            Tscan=0.2e6,
            μ_sim=[-220, 220],
            δω_end = δω, δω_init = δω, 
            DKS_init =  DKS_init, 
            D1_manual = 982.2557817183881 * 1e9 * 2*np.pi
        )
        solver = pyLLE.LLEsolver(sim=sim, res=res, debug=False)

        # fig = solver.Analyze(plot=False)
        solver.Setup(verbose=False)
        print("Starting temporal solver...")
        start_time = time.time()
        try:
            solver.SolveTemporal(bin=JULIA_BIN)
            solver.RetrieveData()
            end_time = time.time()
            print(f"Simulation 4 completed in {end_time - start_time:.2f} seconds")
        except Exception as e:
            print(f"Simulation 4 failed: {e}")

        # fig = solver.PlotCombSpectra(4000)
        # fig.show()


if __name__ == "__main__":
    run_benchmark()
