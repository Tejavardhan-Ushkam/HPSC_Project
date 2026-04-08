#!/usr/bin/env python3
"""
verify_bounce.py
Reads results/trajectory.csv (bounce.cfg run) and plots:
  1. z(t) – bouncing trajectory
  2. Rebound height vs bounce number
  3. Kinetic energy vs time
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

FIGDIR = 'results/figures'
os.makedirs(FIGDIR, exist_ok=True)

df = pd.read_csv('results/trajectory.csv')
df = df[df['pid'] == 0].copy()
t  = df['time'].values
z  = df['z'].values
vz = df['vz'].values

# ── Plot 1: bouncing height vs time ──────────────────────────────────────────
fig, ax = plt.subplots(figsize=(8, 4))
ax.plot(t, z, 'b-', lw=0.8)
ax.set_xlabel('Time (s)', fontsize=12)
ax.set_ylabel('z position (m)', fontsize=12)
ax.set_title('Bouncing Particle: Height vs Time', fontsize=12)
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(f'{FIGDIR}/bounce_trajectory.pdf', dpi=200)
plt.savefig(f'{FIGDIR}/bounce_trajectory.png', dpi=200)
plt.close()

# ── Extract rebound peaks ─────────────────────────────────────────────────────
from scipy.signal import find_peaks
peaks, _ = find_peaks(z, height=0.03, distance=50)
if len(peaks) > 1:
    fig2, ax2 = plt.subplots(figsize=(6, 4))
    bounce_nums   = np.arange(1, len(peaks) + 1)
    rebound_heights = z[peaks]
    ax2.plot(bounce_nums, rebound_heights, 'ro-', ms=6, lw=1.5)
    ax2.set_xlabel('Bounce number', fontsize=12)
    ax2.set_ylabel('Rebound height (m)', fontsize=12)
    ax2.set_title('Rebound Height vs Bounce Number', fontsize=12)
    ax2.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(f'{FIGDIR}/bounce_rebound_heights.pdf', dpi=200)
    plt.savefig(f'{FIGDIR}/bounce_rebound_heights.png', dpi=200)
    plt.close()
    print(f"[verify_bounce] Found {len(peaks)} bounces. Max rebound = {rebound_heights[0]:.3f} m")
else:
    print("[verify_bounce] Too few peaks detected – try longer simulation time.")

# ── KE from KE file ──────────────────────────────────────────────────────────
ke_file = 'results/kinetic_energy.csv'
if os.path.exists(ke_file):
    ke_df = pd.read_csv(ke_file)
    fig3, ax3 = plt.subplots(figsize=(7, 4))
    ax3.plot(ke_df['time'], ke_df['kinetic_energy'], 'g-', lw=1)
    ax3.set_xlabel('Time (s)', fontsize=12)
    ax3.set_ylabel('Kinetic Energy (J)', fontsize=12)
    ax3.set_title('Kinetic Energy vs Time (Bounce Test)', fontsize=12)
    ax3.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(f'{FIGDIR}/bounce_kinetic_energy.pdf', dpi=200)
    plt.savefig(f'{FIGDIR}/bounce_kinetic_energy.png', dpi=200)
    plt.close()

print("[verify_bounce] Done.")
