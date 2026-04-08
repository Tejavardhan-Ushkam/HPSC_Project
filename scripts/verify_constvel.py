#!/usr/bin/env python3
"""
verify_constvel.py
Verifies that a particle with zero gravity travels at constant velocity.
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

FIGDIR = 'results/figures'
os.makedirs(FIGDIR, exist_ok=True)

VZ0 = 2.0   # initial vz set in mode 4
Z0  = 0.1   # initial z

df = pd.read_csv('results/trajectory.csv')
df = df[df['pid'] == 0].copy()
t  = df['time'].values
z  = df['z'].values
vz = df['vz'].values

z_analytical = Z0 + VZ0 * t

fig, axes = plt.subplots(1, 2, figsize=(10, 4))

ax = axes[0]
ax.plot(t, z_analytical, 'k--', lw=1.5, label='Analytical z = z0 + v0*t')
ax.plot(t, z,            'r-',  lw=1,   label='Numerical', alpha=0.8)
ax.set_xlabel('Time (s)', fontsize=12)
ax.set_ylabel('z (m)', fontsize=12)
ax.set_title('Constant-Velocity: Position', fontsize=12)
ax.legend(fontsize=11)
ax.grid(True, alpha=0.3)

# velocity drift
ax2 = axes[1]
v_err = np.abs(vz - VZ0)
ax2.semilogy(t, v_err + 1e-16, 'b-', lw=1.5)
ax2.set_xlabel('Time (s)', fontsize=12)
ax2.set_ylabel('|vz - v0| (m/s)', fontsize=12)
ax2.set_title('Velocity Error', fontsize=12)
ax2.grid(True, alpha=0.3, which='both')

plt.tight_layout()
plt.savefig(f'{FIGDIR}/constvel_verification.pdf', dpi=200)
plt.savefig(f'{FIGDIR}/constvel_verification.png', dpi=200)
plt.close()

max_drift = v_err.max()
print(f"[verify_constvel] Max velocity drift = {max_drift:.3e} m/s")
print("[verify_constvel] Done.")
