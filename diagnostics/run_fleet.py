#!/usr/bin/env python3
"""
H2 Fleet Runner
===============
Automated pipeline for running H2 diagnostics on multiple galaxies.

Usage:
    # Single galaxy
    python -m diagnostics.run_fleet --galaxies NGC3198
    
    # Multiple galaxies
    python -m diagnostics.run_fleet --galaxies NGC3198,IC2574,NGC0891
    
    # All available galaxies
    python -m diagnostics.run_fleet --all
    
    # Custom alpha parameter
    python -m diagnostics.run_fleet --galaxies NGC3198,IC2574 --alpha 2.0
"""

import argparse
import subprocess
import pandas as pd
from pathlib import Path
import sys


def check_prerequisites(galaxy):
    """Check if required input files exist for a galaxy."""
    rotmod = Path(f'data/sparc/{galaxy}_rotmod.dat')
    h1_frozen = Path(f'data/h1_frozen/per_galaxy/rc_decomp_{galaxy}_best.csv')
    
    missing = []
    if not rotmod.exists():
        missing.append(f'SPARC rotmod: {rotmod}')
    if not h1_frozen.exists():
        missing.append(f'H1 frozen: {h1_frozen}')
    
    return missing


def run_command(cmd, galaxy, step_name):
    """Run a command and return status."""
    print(f"\n[{galaxy}] Running {step_name}...")
    print(f"  Command: {' '.join(cmd)}")
    
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=False
        )
        
        if result.returncode != 0:
            print(f"  ✗ FAILED (exit code {result.returncode})")
            if result.stderr:
                print(f"  Error: {result.stderr[:500]}")
            return False
        else:
            print(f"  ✓ Success")
            return True
            
    except Exception as e:
        print(f"  ✗ Exception: {e}")
        return False


def extract_test2_results(galaxy):
    """Extract Test-2 results by running the fixed script."""
    cmd = [
        sys.executable, '-m', 'diagnostics.test2_chi_correlation',
        '--galaxy', galaxy
    ]
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=False)
        if result.returncode != 0:
            return None, None, None
        
        # Parse output
        lines = result.stdout.strip().split('\n')
        r, N, max_dV = None, None, None
        
        for line in lines:
            if 'Pearson r:' in line:
                r = float(line.split(':')[1].strip())
            elif 'N points:' in line:
                N = int(line.split(':')[1].strip())
            elif 'max|ΔV|' in line:
                parts = line.split(':')[1].strip().split()
                max_dV = float(parts[0])
        
        return r, N, max_dV
        
    except Exception as e:
        print(f"  Warning: Could not extract Test-2 results: {e}")
        return None, None, None


def extract_leff_metrics(galaxy):
    """Extract L_eff summary metrics from Phase-3 leff CSV."""
    leff_path = Path(f'data/derived/phase3/leff_{galaxy}.csv')
    try:
        df = pd.read_csv(leff_path)
        col = df['L_eff_kpc']
        return (
            round(float(col.mean()), 4),
            round(float(col.min()),  4),
            round(float(col.max()),  4),
            round(float(col.iloc[-1]), 4),   # L_eff at outermost radius
        )
    except Exception as e:
        print(f"  Warning: Could not extract L_eff metrics for {galaxy}: {e}")
        return None, None, None, None


def extract_test3_results(galaxy):
    """Extract Test-3 results by running the script."""
    cmd = [
        sys.executable, '-m', 'diagnostics.test3_inner_scatter',
        '--galaxy', galaxy
    ]
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=False)
        if result.returncode != 0:
            return None, None, None
        
        # Parse output
        lines = result.stdout.strip().split('\n')
        sigma_h1, sigma_h2, delta_sigma = None, None, None
        
        for line in lines:
            if 'sigma_inner(H1)' in line:
                sigma_h1 = float(line.split('=')[1].strip().split()[0])
            elif 'sigma_inner(H2)' in line:
                sigma_h2 = float(line.split('=')[1].strip().split()[0])
            elif 'Delta sigma' in line:
                delta_sigma = float(line.split('=')[1].strip().split()[0])
        
        return sigma_h1, sigma_h2, delta_sigma
        
    except Exception as e:
        print(f"  Warning: Could not extract Test-3 results: {e}")
        return None, None, None


def process_galaxy(galaxy, alpha, sigma_idx, taper):
    """Run full H2 pipeline for one galaxy."""
    
    print(f"\n{'='*70}")
    print(f"Processing: {galaxy}")
    print(f"{'='*70}")
    
    # Check prerequisites
    missing = check_prerequisites(galaxy)
    if missing:
        print(f"\n✗ SKIPPED: Missing prerequisites:")
        for m in missing:
            print(f"    - {m}")
        return {
            'galaxy': galaxy,
            'status': 'MISSING_PREREQUISITES',
            'test1_pass': None,
            'test2_r': None,
            'test2_N': None,
            'test2_max_dV_kms': None,
            'test3_sigma_h1_dex': None,
            'test3_sigma_h2_dex': None,
            'test3_delta_sigma_dex': None,
            'mean_L_eff_kpc': None,
            'min_L_eff_kpc': None,
            'max_L_eff_kpc': None,
            'L_eff_at_Rmax_kpc': None,
        }
    
    # Phase-3: L_eff computation
    cmd_phase3 = [
        sys.executable, '-m', 'diagnostics.phase3_leff_ngc3198',
        '--galaxy', galaxy,
        '--alpha', str(alpha),
        '--sigma_idx', str(sigma_idx)
    ]
    if taper:
        cmd_phase3.append('--taper')
    
    if not run_command(cmd_phase3, galaxy, 'Phase-3'):
        return {'galaxy': galaxy, 'status': 'FAILED_PHASE3'}
    
    # Phase-4: Adaptive convolution
    cmd_phase4 = [
        sys.executable, '-m', 'diagnostics.phase4_adaptive_convolution',
        '--galaxy', galaxy
    ]
    
    if not run_command(cmd_phase4, galaxy, 'Phase-4'):
        return {'galaxy': galaxy, 'status': 'FAILED_PHASE4'}
    
    # Extract Test-1 result from Phase-4 output
    try:
        h2_file = Path(f'data/derived/phase4/h2_outputs/rc_decomp_{galaxy}_H2_adaptive.csv')
        h2_data = pd.read_csv(h2_file)
        test1_pass = h2_data['outer_pass'].iloc[0] if 'outer_pass' in h2_data.columns else None
    except:
        test1_pass = None
    
    # L_eff metrics from Phase-3 output
    mean_leff, min_leff, max_leff, leff_at_rmax = extract_leff_metrics(galaxy)

    # Test-2: Chi correlation
    print(f"\n[{galaxy}] Running Test-2...")
    test2_r, test2_N, test2_max_dV = extract_test2_results(galaxy)

    # Test-3: Scatter comparison
    print(f"\n[{galaxy}] Running Test-3...")
    test3_sigma_h1, test3_sigma_h2, test3_delta_sigma = extract_test3_results(galaxy)

    # Compile results
    results = {
        'galaxy': galaxy,
        'status': 'OK',
        'test1_pass': test1_pass,
        'test2_r': test2_r,
        'test2_N': test2_N,
        'test2_max_dV_kms': test2_max_dV,
        'test3_sigma_h1_dex': test3_sigma_h1,
        'test3_sigma_h2_dex': test3_sigma_h2,
        'test3_delta_sigma_dex': test3_delta_sigma,
        'mean_L_eff_kpc': mean_leff,
        'min_L_eff_kpc': min_leff,
        'max_L_eff_kpc': max_leff,
        'L_eff_at_Rmax_kpc': leff_at_rmax,
    }
    
    print(f"\n✓ {galaxy} completed successfully")
    
    return results


def main():
    parser = argparse.ArgumentParser(
        description="Run H2 pipeline on multiple galaxies",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python -m diagnostics.run_fleet --galaxies NGC3198,IC2574
  python -m diagnostics.run_fleet --all --alpha 2.0 --taper
        """
    )
    
    parser.add_argument(
        '--galaxies',
        help='Comma-separated list of galaxy names (e.g., NGC3198,IC2574)'
    )
    parser.add_argument(
        '--all',
        action='store_true',
        help='Run on all galaxies found in data/sparc/'
    )
    parser.add_argument(
        '--alpha',
        type=float,
        default=2.0,
        help='Alpha parameter for deformation (default: 2.0)'
    )
    parser.add_argument(
        '--sigma_idx',
        type=float,
        default=1.0,
        help='Sigma index for smoothing (default: 1.0)'
    )
    parser.add_argument(
        '--taper',
        action='store_true',
        help='Enable tapering in Phase-3'
    )
    parser.add_argument(
        '--output',
        default='data/derived/fleet/fleet_summary.csv',
        help='Output CSV file for results summary'
    )
    
    args = parser.parse_args()
    
    # Determine galaxy list
    if args.all:
        sparc_dir = Path('data/sparc')
        if not sparc_dir.exists():
            print(f"✗ Error: SPARC directory not found: {sparc_dir}")
            return 1
        
        galaxy_files = list(sparc_dir.glob('*_rotmod.dat'))
        galaxies = [f.stem.replace('_rotmod', '') for f in galaxy_files]
        print(f"Found {len(galaxies)} galaxies in {sparc_dir}")
        
    elif args.galaxies:
        galaxies = [g.strip() for g in args.galaxies.split(',')]
        
    else:
        print("✗ Error: Must specify either --galaxies or --all")
        parser.print_help()
        return 1
    
    print(f"\n{'='*70}")
    print(f"H2 FLEET RUNNER")
    print(f"{'='*70}")
    print(f"Galaxies to process: {len(galaxies)}")
    print(f"Alpha: {args.alpha}")
    print(f"Sigma index: {args.sigma_idx}")
    print(f"Taper: {args.taper}")
    print(f"Output: {args.output}")
    print(f"{'='*70}\n")
    
    # Process each galaxy
    all_results = []
    
    for i, galaxy in enumerate(galaxies, 1):
        print(f"\n[{i}/{len(galaxies)}] Galaxy: {galaxy}")
        
        result = process_galaxy(
            galaxy=galaxy,
            alpha=args.alpha,
            sigma_idx=args.sigma_idx,
            taper=args.taper
        )
        
        all_results.append(result)
    
    # Save results
    df = pd.DataFrame(all_results)
    
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    df.to_csv(output_path, index=False)
    
    print(f"\n{'='*70}")
    print(f"FLEET SUMMARY")
    print(f"{'='*70}")
    print(f"Total galaxies: {len(galaxies)}")
    print(f"Successful: {(df['status'] == 'OK').sum()}")
    print(f"Failed: {(df['status'] != 'OK').sum()}")
    print(f"\nResults saved to: {output_path}")
    print(f"{'='*70}\n")
    
    # Print summary table
    print("\nQUICK RESULTS:")
    print(df[['galaxy', 'status', 'test2_r', 'test3_delta_sigma_dex']].to_string(index=False))
    
    return 0


if __name__ == "__main__":
    sys.exit(main())
