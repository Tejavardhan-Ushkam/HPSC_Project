#!/usr/bin/env python3
"""Append one row to the scaling table CSV."""
import argparse, csv, os

def main():
    p = argparse.ArgumentParser()
    p.add_argument('--n',       type=int,   required=True)
    p.add_argument('--threads', type=int,   required=True)
    p.add_argument('--timing',  type=str,   default='results/timing.csv')
    p.add_argument('--out',     type=str,   default='results/scaling/scaling_table.csv')
    args = p.parse_args()

    data = {}
    with open(args.timing) as f:
        for row in csv.DictReader(f):
            data[row['phase']] = float(row['time_s'])

    with open(args.out, 'a') as f:
        f.write(f"{args.n},{args.threads},{data.get('total',0):.6f},{data.get('forces',0):.6f}\n")

    print(f"  [scaling] N={args.n} T={args.threads} total={data.get('total',0):.3f}s")

if __name__ == '__main__':
    main()
