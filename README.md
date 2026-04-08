# DEM Particle Simulator – HPSC 2026 Assignment 1

A three-dimensional **Discrete Element Method (DEM)** solver for spherical particles, implemented in **C++** (driver) and **Fortran 90** (compute kernels), parallelised with **OpenMP**, and post-processed with **Python 3**.
# Optimised build (default)
make

# Debug build (bounds checking, no optimisation)
make debug


The binary `dem_sim` is placed in the project root.

---

## Running

### Quick run (N=200 particles, 4 threads)

make run


### Custom particle count / thread count

make scale N=1000 THREADS=8


### Custom config file

./dem_sim configs/default.cfg


Config files use a simple `key = value` format – see `configs/default.cfg` for all options.

---

## Verification Tests

Three analytical tests verify the numerical implementation:


make verify          # runs all three tests below

make test_freefall   # z(t) = z0 - 0.5*g*t^2
make test_bounce     # particle bouncing with damping
make test_constvel   # constant-velocity (g=0) check


Each test runs the simulator and calls the corresponding Python verification script, saving plots to `results/figures/`.

---

## Scaling Study (Strong Scaling)

Runs N=1000 and N=5000 with 1, 2, 4, 8 threads and plots speedup / efficiency:


make scaling_study
# Output: results/figures/scaling_study.pdf


---

## Bonus Experiments

make bonus_damping    # compare gamma_n = 0, 10, 50, 100, 200
make bonus_timestep   # dt convergence study (1e-3 to 1e-6)
make bonus_settling   # N=500 particle cloud settling

---

## Generate All Plots

After any simulation run the bash command : make plot

---

## File Structure

```
dem_project/
├── Makefile                    --single entry point for all targets
├── src/
│   ├── main.cpp                --C++ simulation driver
│   └── dem_physics.f90         --Fortran 90 OpenMP compute kernels
├── include/
│   ├── dem_params.h            --parameter struct + config parser
│   ├── timer.h                 --wall-clock timers
│   └── dem_io.h                --I/O helpers
├── configs/
│   ├── default.cfg             --N=200 multi-particle run
│   ├── freefall.cfg            --free-fall verification
│   ├── bounce.cfg              --bounce verification
│   └── constvel.cfg            --constant-velocity verification
├── scripts/
│   ├── plot_all.py             --master plot script
│   ├── verify_freefall.py      --free-fall analytical comparison
│   ├── verify_bounce.py        --bounce verification plots
│   ├── verify_constvel.py      --constant-velocity verification
│   ├── plot_scaling.py         --speedup / efficiency plots
│   ├── plot_damping.py         --damping sensitivity figures
│   ├── plot_timestep_sensitivity.py -- convergence study figures
│   ├── plot_settling.py        --particle settling figures
│   ├── gen_config.py           --programmatic config generation
│   └── extract_timing.py       -- timing extraction for scaling study
├── report/
│   └── main.tex                -- IEEE-format LaTeX white paper
└── results/                    -- generated outputs (gitignored)
```

---

## LLM Usage Disclosure

AI assistance was used to debugg the initial codes.

---

## Clean Up

```bash
make clean      # remove compiled objects and binary
make cleanall   # also remove results/
```
