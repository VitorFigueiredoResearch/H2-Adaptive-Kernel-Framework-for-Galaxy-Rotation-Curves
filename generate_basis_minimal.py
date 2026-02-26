#!/usr/bin/env python3
"""
Minimal H2 Basis Generator
===========================
Directly calls H1's predict_rc_for_params() to generate basis files
without any fitting logic.

Usage:
    python generate_basis_minimal.py --L_values 10 30 50 80
"""

import sys
import os
import argparse
from pathlib import Path
import pandas as pd

# Add H1 source to path
H2_ROOT = Path(__file__).resolve().parent
H1_SRC = H2_ROOT / "vendor" / "h1_src"
sys.path.insert(0, str(H1_SRC))

# Import H1 core function
from run_sparc_lite import predict_rc_for_params

def main():
    parser = argparse.ArgumentParser(description="Generate H2 basis files")
    parser.add_argument("--galaxy", default="NGC3198", help="Galaxy name")
    parser.add_argument("--L_values", nargs="+", type=float, required=True,
                       help="L values in kpc (e.g., 10 30 50 80)")
    parser.add_argument("--mu", type=float, default=10.0, 
                       help="Coupling parameter (default: 10.0)")
    parser.add_argument("--kernel", default="ananta-hybrid",
                       help="Kernel type (default: ananta-hybrid)")
    parser.add_argument("--beta", type=float, default=1.15,
                       help="Kernel beta parameter (default: 1.15)")
    args = parser.parse_args()
    
    print(f"\n{'='*60}")
    print(f"MINIMAL H2 BASIS GENERATION")
    print(f"{'='*60}")
    print(f"Galaxy: {args.galaxy}")
    print(f"L values: {args.L_values} kpc")
    print(f"mu: {args.mu}")
    print(f"kernel: {args.kernel}")
    print(f"beta: {args.beta}")
    print(f"{'='*60}\n")
    
    # Hardcode NGC3198 parameters (from SPARC database)
    if args.galaxy != "NGC3198":
        print(f"✗ Error: Only NGC3198 is currently supported")
        print(f"  (Other galaxies require CSV parsing fix)")
        return 1
    
    gal = {
        "name": "NGC3198",
        "Rd_star": 3.14,      # kpc
        "Mstar": 1.91e10,     # Msun
        "hz_star": 0.628,     # kpc
        "Rd_gas": 5.65,       # kpc  
        "Mgas": 1.45e10,      # Msun
        "hz_gas": 0.15        # kpc (assumed)
    }
    
    print(f"✓ Found galaxy: {gal['name']}")
    print(f"  Mstar: {gal['Mstar']:.2e} Msun")
    print(f"  Mgas: {gal['Mgas']:.2e} Msun")
    print(f"  Rd_star: {gal['Rd_star']:.3f} kpc")
    
    # Create output directory
    output_dir = H2_ROOT / "data" / "derived" / "phase4" / "basis_h1"
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Generate basis files
    for L in args.L_values:
        print(f"\n--- Generating L = {L} kpc ---")
        
        try:
            # Call H1's core function directly
            result = predict_rc_for_params(
                gal=gal,
                L=L,
                mu=args.mu,
                kernel=args.kernel,
                beta=args.beta
            )
            
            if result is None:
                print(f"✗ Error: predict_rc_for_params returned None")
                continue
            
            R_pred, V_b, V_k, V_pred = result
            
            # Verify non-zero velocities
            if V_b.max() == 0.0 or V_pred.max() == 0.0:
                print(f"✗ Error: Generated velocities are all zero")
                print(f"  V_baryon max: {V_b.max()}")
                print(f"  V_total max: {V_pred.max()}")
                continue
            
            # Create output filename
            L_tag = f"{int(round(L))}kpc" if abs(L - round(L)) < 1e-9 else f"{L:.3f}kpc"
            filename = f"rc_decomp_{args.galaxy}_L{L_tag}.csv"
            filepath = output_dir / filename
            
            # Save to CSV (match H1 format exactly)
            df = pd.DataFrame({
                'R_kpc': R_pred,
                'V_baryon': V_b,
                'V_kernel': V_k,
                'V_total': V_pred
            })
            df.to_csv(filepath, index=False)
            
            print(f"✓ Saved: {filename}")
            print(f"  R range: [{R_pred.min():.2f}, {R_pred.max():.2f}] kpc")
            print(f"  V_baryon: [{V_b.min():.2f}, {V_b.max():.2f}] km/s")
            print(f"  V_total: [{V_pred.min():.2f}, {V_pred.max():.2f}] km/s")
            
        except Exception as e:
            print(f"✗ Error generating L={L}: {e}")
            import traceback
            traceback.print_exc()
            continue
    
    print(f"\n{'='*60}")
    print(f"GENERATION COMPLETE")
    print(f"{'='*60}\n")
    
    return 0


if __name__ == "__main__":
    sys.exit(main())
