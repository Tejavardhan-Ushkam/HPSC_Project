
#!/usr/bin/env python3
"""
plot_timestep_sensitivity.py
Compares free-fall numerical accuracy across multiple dt values.
Reads from results/dt_<val>/ sub-folders.
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.rcParams.update({'font.size': 12, 'axes.labelsize': 12,
                            'axes.titlesize': 12, 'legend.fontsize': 10,
                            'xtick.labelsize': 11, 'ytick.labelsize': 11,
                            'axes.grid': True, 'grid.alpha': 0.3})
import matplotlib.pyplot as plt
import os

FIGDIR = 'results/figures'
os.makedirs(FIGDIR, exist_ok=True)

G   = 9.81
Z0  = 0.9       # initial z (mode 2 freefall)
DT_VALS = [1e-3, 1e-4, 1e-5, 1e-6]
COLOURS = ['crimson', 'darkorange', 'steelblue', 'seagreen']

def dt_label(dt):
    exp = int(round(np.log10(dt)))
    return f'Δt=10^{{{exp}}}'

def folder_name(dt):
    exp = int(round(np.log10(dt)))
    return f'results/dt_1e{exp}'

def analytical_z(t):
    return Z0 - 0.5 * G * t**2

# ── Figure 1: trajectories vs analytical ─────────────────────────────────────
fig1, ax1 = plt.subplots(figsize=(8, 4))

max_errors = []
dt_actuals = []

for dt, c in zip(DT_VALS, COLOURS):
    fdir  = folder_name(dt)
    fpath = f'{fdir}/trajectory.csv'
    if not os.path.exists(fpath):
        print(f"  [warn] {fpath} not found – skipping")
        continue
    df = pd.read_csv(fpath)
    df = df[df['pid'] == 0].copy()
    t  = df['time'].values
    z  = df['z'].values

    # Only pre-contact phase
    mask = z > 0.025
    t_c, z_c = t[mask], z[mask]
    z_a = analytical_z(t_c)
    err = np.abs(z_c - z_a)

    ax1.plot(t_c, z_c, color=c, lw=1, alpha=0.85, label=f'$\\Delta t=10^{{{int(round(np.log10(dt)))}}}$')
    max_errors.append(err.max())
    dt_actuals.append(dt)

# Analytical reference
t_ref = np.linspace(0, 0.4, 500)
z_ref = analytical_z(t_ref)
mask_ref = z_ref > 0.025
ax1.plot(t_ref[mask_ref], z_ref[mask_ref], 'k--', lw=1.5, label='Analytical')

ax1.set_xlabel('Time (s)')
ax1.set_ylabel('z position (m)')
ax1.set_title('Free-Fall: Numerical vs Analytical for Multiple Timesteps')
ax1.legend(loc='lower left')
plt.tight_layout()
plt.savefig(f'{FIGDIR}/timestep_trajectories.pdf', dpi=200, bbox_inches='tight')
plt.savefig(f'{FIGDIR}/timestep_trajectories.png', dpi=200, bbox_inches='tight')
plt.close()
print("  saved timestep_trajectories")

# ── Figure 2: max error vs dt (convergence plot) ─────────────────────────────
if len(max_errors) >= 2:
    fig2, ax2 = plt.subplots(figsize=(6, 4))
    dt_arr  = np.array(dt_actuals)
    err_arr = np.array(max_errors)

    ax2.loglog(dt_arr, err_arr, 'ko-', ms=7, lw=1.5, label='Max |error|')

    # Reference O(dt) line
    ref_x = np.array([dt_arr.min(), dt_arr.max()])
    ref_y = err_arr[-1] * (ref_x / dt_arr[-1])      # passes through finest point
    ax2.loglog(ref_x, ref_y, 'r--', lw=1.2, label='O(Δt) reference')

    ax2.set_xlabel('Timestep Δt (s)')
    ax2.set_ylabel('Max |z_num − z_ana| (m)')
    ax2.set_title('Convergence: Position Error vs Timestep')
    ax2.legend()
    plt.tight_layout()
    plt.savefig(f'{FIGDIR}/timestep_error_vs_dt.pdf', dpi=200, bbox_inches='tight')
    plt.savefig(f'{FIGDIR}/timestep_error_vs_dt.png', dpi=200, bbox_inches='tight')
    plt.close()
    print("  saved timestep_error_vs_dt")

    # Print measured convergence order
    orders = np.diff(np.log(err_arr)) / np.diff(np.log(dt_arr))
    print(f"  Measured convergence orders: {np.round(orders, 2)}")

print(f"\n[plot_timestep_sensitivity] Done – figures in {FIGDIR}/")
