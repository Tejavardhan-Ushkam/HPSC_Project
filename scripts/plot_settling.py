#!/usr/bin/env python3
"""
plot_settling.py
Visualises the settling of a particle cloud under gravity.
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

def save(name):
    plt.tight_layout()
    plt.savefig(f'{FIGDIR}/{name}.pdf', dpi=200, bbox_inches='tight')
    plt.savefig(f'{FIGDIR}/{name}.png', dpi=200, bbox_inches='tight')
    plt.close()
    print(f"  saved {name}")

# ── Load data ─────────────────────────────────────────────────────────────────
traj_path = 'results/trajectory.csv'
ke_path   = 'results/kinetic_energy.csv'

if not os.path.exists(traj_path):
    print(f"[plot_settling] {traj_path} not found. Run: make run")
    raise SystemExit(1)

traj = pd.read_csv(traj_path)
times_all = sorted(traj['time'].unique())

# ── Figure 1: Centre-of-mass height vs time ──────────────────────────────────
com_z = traj.groupby('time')['z'].mean().reset_index()
com_z.columns = ['time', 'com_z']

fig, ax = plt.subplots(figsize=(7, 4))
ax.plot(com_z['time'], com_z['com_z'], 'steelblue', lw=1.5)
ax.set_xlabel('Time (s)')
ax.set_ylabel('Centre-of-mass z (m)')
ax.set_title('Particle Cloud: Centre-of-Mass Height vs Time')
save('settling_com_height')

# ── Figure 2: Kinetic energy ──────────────────────────────────────────────────
if os.path.exists(ke_path):
    ke = pd.read_csv(ke_path)
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(ke['time'], ke['kinetic_energy'], 'tomato', lw=1.5)
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('Total Kinetic Energy (J)')
    ax.set_title('Particle Cloud: Kinetic Energy During Settling')
    save('settling_ke')

# ── Figure 3: 4-panel snapshots ───────────────────────────────────────────────
snap_fracs = [0.0, 0.25, 0.5, 1.0]
snap_times = [times_all[min(int(f * (len(times_all)-1)), len(times_all)-1)]
              for f in snap_fracs]

fig, axes = plt.subplots(1, 4, figsize=(14, 4), sharex=True, sharey=True)
for ax, st in zip(axes, snap_times):
    # get closest time
    closest = times_all[np.argmin(np.abs(np.array(times_all) - st))]
    snap = traj[np.isclose(traj['time'], closest, atol=1e-9)]
    sc = ax.scatter(snap['x'], snap['z'], s=8, alpha=0.65,
                    c=snap['z'], cmap='viridis_r', vmin=0, vmax=1)
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.set_xlabel('x (m)')
    ax.set_title(f't = {closest:.3f} s')
    ax.set_aspect('equal')

axes[0].set_ylabel('z (m)')
plt.colorbar(sc, ax=axes[-1], label='z (m)', shrink=0.85)
fig.suptitle('Particle Configuration Snapshots (x–z plane)', fontsize=13, y=1.02)
save('settling_snapshots')

# ── Figure 4: Velocity magnitude distribution at final time ──────────────────
tf   = times_all[-1]
snap_f = traj[np.isclose(traj['time'], tf, atol=1e-9)].copy()
snap_f['speed'] = np.sqrt(snap_f['vx']**2 + snap_f['vy']**2 + snap_f['vz']**2)

fig, ax = plt.subplots(figsize=(6, 4))
ax.hist(snap_f['speed'], bins=30, color='steelblue', edgecolor='white', lw=0.5)
ax.set_xlabel('Particle speed (m/s)')
ax.set_ylabel('Count')
ax.set_title(f'Speed Distribution at t = {tf:.3f} s')
save('settling_speed_hist')

print(f"\n[plot_settling] All figures saved to {FIGDIR}/")
