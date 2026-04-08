#!/usr/bin/env python3
"""
verify_freefall.py
Reads results/trajectory.csv (from freefall.cfg run) and compares
numerical z(t) with the analytical solution z(t) = z0 - 0.5*g*t^2.
Also plots error vs timestep for multiple dt values if available.
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

FIGDIR = 'results/figures'
os.makedirs(FIGDIR, exist_ok=True)

G  = 9.81
Z0 = 0.9   # initial z (fraction of Lz=1.0)  –  see init_particles mode 2

def analytical_z(t):
    return Z0 - 0.5 * G * t**2

def analytical_vz(t):
    return -G * t

df = pd.read_csv('results/trajectory.csv')
# Keep only particle 0
df = df[df['pid'] == 0].copy()
t  = df['time'].values
z_num = df['z'].values
vz_num = df['vz'].values

z_ana  = analytical_z(t)
vz_ana = analytical_vz(t)

# Only compare while particle hasn't hit floor
mask = z_num > 0.02   # above floor (radius = 0.02)
t_c    = t[mask]
z_c    = z_num[mask]
z_a_c  = z_ana[mask]

# ── Plot 1: trajectory comparison ────────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(10, 4))

ax = axes[0]
ax.plot(t_c, z_a_c,  'k--', lw=1.5, label='Analytical')
ax.plot(t_c, z_c,    'r-',  lw=1,   label='Numerical', alpha=0.8)
ax.set_xlabel('Time (s)', fontsize=12)
ax.set_ylabel('z position (m)', fontsize=12)
ax.set_title('Free-Fall: z(t)', fontsize=12)
ax.legend(fontsize=11)
ax.grid(True, alpha=0.3)

# ── Plot 2: error vs time ─────────────────────────────────────────────────────
err = np.abs(z_c - z_a_c)
ax2 = axes[1]
ax2.semilogy(t_c, err + 1e-16, 'b-', lw=1.5)
ax2.set_xlabel('Time (s)', fontsize=12)
ax2.set_ylabel('|z_num - z_ana| (m)', fontsize=12)
ax2.set_title('Free-Fall: Position Error', fontsize=12)
ax2.grid(True, alpha=0.3, which='both')

plt.tight_layout()
plt.savefig(f'{FIGDIR}/freefall_comparison.pdf', dpi=200)
plt.savefig(f'{FIGDIR}/freefall_comparison.png', dpi=200)
plt.close()

print(f"[verify_freefall] Max z error = {err.max():.3e} m")
print(f"[verify_freefall] Figures saved to {FIGDIR}/")

# ── Velocity comparison ───────────────────────────────────────────────────────
fig2, ax3 = plt.subplots(figsize=(6, 4))
ax3.plot(t_c, analytical_vz(t_c), 'k--', lw=1.5, label='Analytical vz')
ax3.plot(t_c, vz_num[mask],        'r-',  lw=1,   label='Numerical vz', alpha=0.8)
ax3.set_xlabel('Time (s)', fontsize=12)
ax3.set_ylabel('vz (m/s)', fontsize=12)
ax3.set_title('Free-Fall: Velocity', fontsize=12)
ax3.legend(fontsize=11)
ax3.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(f'{FIGDIR}/freefall_velocity.pdf', dpi=200)
plt.savefig(f'{FIGDIR}/freefall_velocity.png', dpi=200)
plt.close()

print("[verify_freefall] Done.")
