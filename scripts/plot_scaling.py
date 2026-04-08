#!/usr/bin/env python3
"""
plot_scaling.py
Reads results/scaling_table.csv and produces:
  - Speedup S(p) vs threads
  - Efficiency E(p) vs threads
  - Runtime vs threads (log scale)
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

FIGDIR = 'results/figures'
os.makedirs(FIGDIR, exist_ok=True)

csv_path = 'results/scaling_table.csv'
if not os.path.exists(csv_path):
    print(f"[plot_scaling] {csv_path} not found. Run 'make scaling_study' first.")
    raise SystemExit(1)

df = pd.read_csv(csv_path)
n_vals = sorted(df['n_particles'].unique())

fig, axes = plt.subplots(1, 3, figsize=(14, 4))

for n_val in n_vals:
    sub = df[df['n_particles'] == n_val].sort_values('threads')
    threads = sub['threads'].values
    times   = sub['total_s'].values
    T1      = times[0]         # serial time (1 thread)
    speedup = T1 / times
    eff     = speedup / threads

    axes[0].plot(threads, times,   'o-', lw=1.5, ms=5, label=f'N={n_val}')
    axes[1].plot(threads, speedup, 'o-', lw=1.5, ms=5, label=f'N={n_val}')
    axes[2].plot(threads, eff*100, 'o-', lw=1.5, ms=5, label=f'N={n_val}')

# Ideal lines
max_t = max(df['threads'])
t_range = np.array(sorted(df['threads'].unique()))
axes[1].plot(t_range, t_range, 'k--', lw=1, label='Ideal')
axes[2].axhline(100, color='k', ls='--', lw=1, label='Ideal 100%')

for ax, ylabel, title in zip(
        axes,
        ['Total runtime (s)', 'Speedup S(p)', 'Efficiency E(p) (%)'],
        ['Runtime vs Threads', 'Speedup vs Threads', 'Efficiency vs Threads']):
    ax.set_xlabel('Number of threads', fontsize=12)
    ax.set_ylabel(ylabel, fontsize=12)
    ax.set_title(title, fontsize=12)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.set_xticks(t_range)

plt.tight_layout()
plt.savefig(f'{FIGDIR}/scaling_study.pdf', dpi=200)
plt.savefig(f'{FIGDIR}/scaling_study.png', dpi=200)
plt.close()
print(f"[plot_scaling] Saved {FIGDIR}/scaling_study.pdf")
