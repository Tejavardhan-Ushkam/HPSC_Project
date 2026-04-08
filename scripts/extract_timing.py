#!/usr/bin/env python3
"""
extract_timing.py – reads results/timing.csv and appends one row to
results/scaling_table.csv  (created by scaling_study Makefile target).
"""
import argparse, csv, os

def main():
    p = argparse.ArgumentParser()
    p.add_argument('--n',       type=int, required=True)
    p.add_argument('--threads', type=int, required=True)
    args = p.parse_args()

    timing_file = 'results/timing.csv'
    out_file    = 'results/scaling_table.csv'

    if not os.path.exists(timing_file):
        print(f"[extract_timing] {timing_file} not found, skipping.")
        return

    data = {}
    with open(timing_file) as f:
        reader = csv.DictReader(f)
        for row in reader:
            data[row['phase']] = float(row['time_s'])

    total  = data.get('total',  0.0)
    forces = data.get('forces', 0.0)

    header = not os.path.exists(out_file)
    with open(out_file, 'a') as f:
        if header:
            f.write('n_particles,threads,total_s,forces_s\n')
        f.write(f"{args.n},{args.threads},{total:.4f},{forces:.4f}\n")

    print(f"[extract_timing] N={args.n} T={args.threads} total={total:.4f}s")

if __name__ == '__main__':
    main()
