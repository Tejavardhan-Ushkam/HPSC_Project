#!/usr/bin/env python3
"""
plot_all.py
Master plotting script – generates all publication-quality figures
for the white paper from whatever results/ files are present.
Run after: make run  (or make test_freefall, test_bounce, etc.)
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.rcParams.update({
    'font.size': 12,
    'axes.labelsize': 12,
    'axes.titlesize': 12,
    'legend.fontsize': 10,
    'xtick.labelsize': 11,
    'ytick.labelsize': 11,
    'figure.dpi': 150,
    'lines.linewidth': 1.5,
    'axes.grid': True,
    'grid.alpha': 0.3,
})
import matplotlib.pyplot as plt
import os, sys

FIGDIR = 'results/figures'
os.makedirs(FIGDIR, exist_ok=True)

def save(name):
    plt.tight_layout()
    plt.savefig(f'{FIGDIR}/{name}.pdf', dpi=200, bbox_inches='tight')
    plt.savefig(f'{FIGDIR}/{name}.png', dpi=200, bbox_inches='tight')
    plt.close()
    print(f"  saved {FIGDIR}/{name}.pdf")

# ══════════════════════════════════════════════════════════════════════════════
# 1. Kinetic Energy vs Time  (multi-particle run)
# ══════════════════════════════════════════════════════════════════════════════
ke_csv = 'results/kinetic_energy.csv'
if os.path.exists(ke_csv):
    ke = pd.read_csv(ke_csv)
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(ke['time'], ke['kinetic_energy'], 'steelblue')
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('Total Kinetic Energy (J)')
    ax.set_title('Kinetic Energy vs Time')
    save('kinetic_energy')

    # Contact count
    if 'n_contacts' in ke.columns:
        fig, ax = plt.subplots(figsize=(7, 4))
        ax.plot(ke['time'], ke['n_contacts'], 'darkorange')
        ax.set_xlabel('Time (s)')
        ax.set_ylabel('Number of Contacts')
        ax.set_title('Active Particle Contacts vs Time')
        save('contact_count')

# ══════════════════════════════════════════════════════════════════════════════
# 2. Timing / profiling pie chart
# ══════════════════════════════════════════════════════════════════════════════
timing_csv = 'results/timing.csv'
if os.path.exists(timing_csv):
    tm = pd.read_csv(timing_csv)
    tm = tm[tm['phase'] != 'total']
    fig, ax = plt.subplots(figsize=(5, 5))
    ax.pie(tm['time_s'], labels=tm['phase'], autopct='%1.1f%%',
           startangle=90, textprops={'fontsize': 11})
    ax.set_title('Runtime Distribution by Phase')
    save('timing_pie')

    # Bar chart of timings
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.bar(tm['phase'], tm['time_s'], color=['steelblue','tomato','seagreen'])
    ax.set_xlabel('Phase')
    ax.set_ylabel('Time (s)')
    ax.set_title('Runtime per Phase')
    save('timing_bar')

# ══════════════════════════════════════════════════════════════════════════════
# 3. Trajectory snapshots (x-z plane, first 20 particles)
# ══════════════════════════════════════════════════════════════════════════════
traj_csv = 'results/trajectory.csv'
if os.path.exists(traj_csv):
    traj = pd.read_csv(traj_csv)
    times_all = sorted(traj['time'].unique())

    # Pick 4 snapshot times
    snap_times = [times_all[int(len(times_all)*f)]
                  for f in [0.0, 0.25, 0.5, 1.0]]

    fig, axes = plt.subplots(1, 4, figsize=(14, 4), sharex=True, sharey=True)
    for ax, st in zip(axes, snap_times):
        snap = traj[np.isclose(traj['time'], st, atol=1e-8)]
        ax.scatter(snap['x'], snap['z'], s=10, alpha=0.7, c='steelblue')
        ax.set_xlim(0, 1); ax.set_ylim(0, 1)
        ax.set_xlabel('x (m)')
        ax.set_title(f't = {st:.3f} s')
        ax.set_aspect('equal')
    axes[0].set_ylabel('z (m)')
    fig.suptitle('Particle Configuration Snapshots (x-z plane)', fontsize=12)
    save('snapshots_xz')

    # 3D scatter at final time (matplotlib)
    tf = times_all[-1]
    snap_f = traj[np.isclose(traj['time'], tf, atol=1e-8)]
    fig = plt.figure(figsize=(6, 6))
    ax3d = fig.add_subplot(111, projection='3d')
    ax3d.scatter(snap_f['x'], snap_f['y'], snap_f['z'],
                 s=10, alpha=0.6, c='steelblue')
    ax3d.set_xlabel('x (m)'); ax3d.set_ylabel('y (m)'); ax3d.set_zlabel('z (m)')
    ax3d.set_title(f'Particle Configuration at t = {tf:.3f} s')
    save('snapshot_3d_final')

# ══════════════════════════════════════════════════════════════════════════════
# 4. Scaling study (if available)
# ══════════════════════════════════════════════════════════════════════════════
scale_csv = 'results/scaling_table.csv'
if os.path.exists(scale_csv):
    import subprocess
    subprocess.call(['python3', 'scripts/plot_scaling.py'])

print(f"\n[plot_all] All figures written to {FIGDIR}/")
