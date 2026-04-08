#!/usr/bin/env python3
"""
plot_damping.py
Compares bouncing particle behaviour across different damping coefficients.
Reads from results/damping_gN/ sub-folders.
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.rcParams.update({'font.size': 12, 'axes.labelsize': 12,
                            'axes.titlesize': 12, 'legend.fontsize': 10,
                            'xtick.labelsize': 11, 'ytick.labelsize': 11,
                            'axes.grid': True, 'grid.alpha': 0.3})
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
import os

FIGDIR = 'results/figures'
os.makedirs(FIGDIR, exist_ok=True)

GAMMA_VALS = [0, 10, 50, 100, 200]
COLOURS    = ['steelblue', 'darkorange', 'seagreen', 'crimson', 'purple']

def load_traj(gamma):
    path = f'results/damping_g{gamma}/trajectory.csv'
    if not os.path.exists(path):
        print(f"  [warn] {path} not found – skipping gamma={gamma}")
        return None
    df = pd.read_csv(path)
    return df[df['pid'] == 0].copy()

def load_ke(gamma):
    path = f'results/damping_g{gamma}/kinetic_energy.csv'
    if not os.path.exists(path):
        return None
    return pd.read_csv(path)

# ── Figure 1: z(t) for all gamma values ─────────────────────────────────────
fig1, ax1 = plt.subplots(figsize=(8, 4))
for g, c in zip(GAMMA_VALS, COLOURS):
    df = load_traj(g)
    if df is None:
        continue
    ax1.plot(df['time'], df['z'], color=c, lw=1.2, label=f'γₙ = {g}')
ax1.set_xlabel('Time (s)')
ax1.set_ylabel('z position (m)')
ax1.set_title('Bouncing Height vs Time for Various Damping Coefficients')
ax1.legend(loc='upper right')
plt.tight_layout()
plt.savefig(f'{FIGDIR}/damping_height_vs_time.pdf', dpi=200, bbox_inches='tight')
plt.savefig(f'{FIGDIR}/damping_height_vs_time.png', dpi=200, bbox_inches='tight')
plt.close()
print(f"  saved damping_height_vs_time")

# ── Figure 2: Rebound heights vs bounce number ───────────────────────────────
fig2, ax2 = plt.subplots(figsize=(7, 4))
for g, c in zip(GAMMA_VALS, COLOURS):
    df = load_traj(g)
    if df is None:
        continue
    z = df['z'].values
    peaks, _ = find_peaks(z, height=0.03, distance=30)
    if len(peaks) < 2:
        continue
    heights = z[peaks]
    ax2.plot(np.arange(1, len(heights)+1), heights, 'o-',
             color=c, ms=5, lw=1.2, label=f'γₙ = {g}')
ax2.set_xlabel('Bounce Number')
ax2.set_ylabel('Rebound Height (m)')
ax2.set_title('Rebound Height vs Bounce Number')
ax2.legend()
plt.tight_layout()
plt.savefig(f'{FIGDIR}/damping_rebound_heights.pdf', dpi=200, bbox_inches='tight')
plt.savefig(f'{FIGDIR}/damping_rebound_heights.png', dpi=200, bbox_inches='tight')
plt.close()
print(f"  saved damping_rebound_heights")

# ── Figure 3: KE vs time ─────────────────────────────────────────────────────
fig3, ax3 = plt.subplots(figsize=(8, 4))
for g, c in zip(GAMMA_VALS, COLOURS):
    ke = load_ke(g)
    if ke is None:
        continue
    ax3.plot(ke['time'], ke['kinetic_energy'], color=c, lw=1.2, label=f'γₙ = {g}')
ax3.set_xlabel('Time (s)')
ax3.set_ylabel('Kinetic Energy (J)')
ax3.set_title('Kinetic Energy vs Time for Various Damping Coefficients')
ax3.legend()
plt.tight_layout()
plt.savefig(f'{FIGDIR}/damping_ke_comparison.pdf', dpi=200, bbox_inches='tight')
plt.savefig(f'{FIGDIR}/damping_ke_comparison.png', dpi=200, bbox_inches='tight')
plt.close()
print(f"  saved damping_ke_comparison")

print(f"\n[plot_damping] All figures saved to {FIGDIR}/")
