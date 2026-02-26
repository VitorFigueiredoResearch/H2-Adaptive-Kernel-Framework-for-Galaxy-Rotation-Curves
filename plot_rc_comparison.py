#!/usr/bin/env python3
"""
Rotation Curve Comparison Plotter
==================================
Plots V_obs vs V_H1 vs V_H2 to visually assess adaptive kernel performance.

Usage:
    python plot_rc_comparison.py --galaxy NGC3198
"""

import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def load_observed_data(galaxy: str):
    """Load SPARC observed rotation curve."""
    rotmod_file = Path('data') / 'sparc' / f'{galaxy}_rotmod.dat'
    
    if not rotmod_file.exists():
        raise FileNotFoundError(f"SPARC data not found: {rotmod_file}")
    
    R_obs = []
    V_obs = []
    V_err = []
    
    with open(rotmod_file) as f:
        for line in f:
            if line.strip() and not line.startswith('#'):
                parts = line.split()
                if len(parts) >= 3:
                    try:
                        R_obs.append(float(parts[0]))
                        V_obs.append(float(parts[1]))
                        V_err.append(float(parts[2]))
                    except ValueError:
                        continue
    
    return np.array(R_obs), np.array(V_obs), np.array(V_err)


def main():
    parser = argparse.ArgumentParser(description="Plot rotation curve comparison")
    parser.add_argument("--galaxy", required=True, help="Galaxy name (e.g., NGC3198)")
    parser.add_argument("--output", help="Output filename (default: auto)")
    args = parser.parse_args()
    
    galaxy = args.galaxy
    
    print(f"\n{'='*60}")
    print(f"ROTATION CURVE COMPARISON: {galaxy}")
    print(f"{'='*60}\n")
    
    # Load observed data
    try:
        R_obs, V_obs, V_err = load_observed_data(galaxy)
        print(f"✓ Loaded observed data: {len(R_obs)} points")
    except FileNotFoundError as e:
        print(f"✗ Error: {e}")
        return 1
    
    # Load H1 frozen
    h1_file = Path('data') / 'h1_frozen' / 'per_galaxy' / f'rc_decomp_{galaxy}_best.csv'
    if not h1_file.exists():
        print(f"✗ Error: H1 frozen not found: {h1_file}")
        return 1
    
    h1 = pd.read_csv(h1_file)
    print(f"✓ Loaded H1 frozen: {len(h1)} points")
    
    # Load H2 adaptive
    h2_file = Path('data') / 'derived' / 'phase4' / 'h2_outputs' / f'rc_decomp_{galaxy}_H2_adaptive.csv'
    if not h2_file.exists():
        print(f"✗ Error: H2 adaptive not found: {h2_file}")
        print(f"  Did you run Phase-4 for this galaxy?")
        return 1
    
    h2 = pd.read_csv(h2_file)
    print(f"✓ Loaded H2 adaptive: {len(h2)} points")
    
    # Extract velocities
    # H2 file contains both H1 and H2 curves
    R_h2 = h2['R_kpc'].values
    V_baryon = h2['V_baryon'].values
    V_total_h1 = h2['V_total_H1'].values
    V_total_h2 = h2['V_total_H2'].values
    deltaV = h2['dV_H2_minus_H1'].values
    
    # Use H2 file for everything (it has all data)
    R_h1 = R_h2
    R_delta = R_h2
    
    # Create figure with 2 subplots
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10), sharex=True)
    
    # --- Top panel: Rotation curves ---
    ax1.errorbar(R_obs, V_obs, yerr=V_err, fmt='o', color='black', 
                 markersize=4, capsize=3, alpha=0.7, label='V_obs (SPARC)')
    
    ax1.plot(R_h1, V_baryon, ':', color='green', linewidth=2, 
             label='V_baryon (baryons only)')
    
    ax1.plot(R_h1, V_total_h1, '-', color='blue', linewidth=2.5, 
             label='V_H1 (frozen kernel)')
    
    ax1.plot(R_h2, V_total_h2, '-', color='orange', linewidth=2.5, 
             label='V_H2 (adaptive kernel)')
    
    ax1.set_ylabel('V (km/s)', fontsize=12)
    ax1.set_title(f'Rotation Curve Comparison: {galaxy}', fontsize=14, fontweight='bold')
    ax1.legend(loc='best', fontsize=10)
    ax1.grid(True, alpha=0.3)
    
    # --- Bottom panel: ΔV = V_H2 - V_H1 ---
    ax2.plot(R_delta, deltaV, '-', color='red', linewidth=2)
    ax2.axhline(0, color='black', linestyle='--', linewidth=1, alpha=0.5)
    ax2.fill_between(R_delta, 0, deltaV, alpha=0.3, color='red')
    
    ax2.set_xlabel('R (kpc)', fontsize=12)
    ax2.set_ylabel('ΔV = V_H2 - V_H1 (km/s)', fontsize=12)
    ax2.set_title('Adaptive Kernel Effect', fontsize=12, fontweight='bold')
    ax2.grid(True, alpha=0.3)
    
    # Add statistics text
    deltaV_max = np.abs(deltaV).max()
    deltaV_mean = np.mean(deltaV)
    
    stats_text = f'max|ΔV| = {deltaV_max:.2f} km/s\nmean ΔV = {deltaV_mean:.2f} km/s'
    ax2.text(0.02, 0.98, stats_text, transform=ax2.transAxes,
             verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5),
             fontsize=10, family='monospace')
    
    plt.tight_layout()
    
    # Save figure
    if args.output:
        output_file = Path(args.output)
    else:
        output_dir = Path('data') / 'derived' / 'phase4' / 'h2_outputs'
        output_dir.mkdir(parents=True, exist_ok=True)
        output_file = output_dir / f'rc_comparison_{galaxy}.png'
    
    plt.savefig(output_file, dpi=200, bbox_inches='tight')
    print(f"\n✓ Saved plot: {output_file}")
    
    # Print summary statistics
    print(f"\n{'='*60}")
    print(f"SUMMARY STATISTICS")
    print(f"{'='*60}")
    print(f"Observed data points: {len(R_obs)}")
    print(f"R range: [{R_obs.min():.2f}, {R_obs.max():.2f}] kpc")
    print(f"V_obs range: [{V_obs.min():.2f}, {V_obs.max():.2f}] km/s")
    print(f"")
    print(f"ΔV statistics (V_H2 - V_H1):")
    print(f"  max|ΔV| = {deltaV_max:.2f} km/s")
    print(f"  mean ΔV = {deltaV_mean:.2f} km/s")
    print(f"  std ΔV  = {np.std(deltaV):.2f} km/s")
    print(f"{'='*60}\n")
    
    plt.show()
    
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())
