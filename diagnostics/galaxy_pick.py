from __future__ import annotations

import argparse

from core.galaxy_io import list_galaxies, load_galaxy_rc, get_Rd_star_kpc


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--galaxy", help="e.g. NGC3198")
    ap.add_argument("--list", action="store_true", help="List available galaxy names")
    args = ap.parse_args()

    if args.list:
        gals = list_galaxies()
        print(f"[Galaxy] Found {len(gals)} galaxies")
        print("[Galaxy] First 20:", ", ".join(gals[:20]))
        return

    if not args.galaxy:
        raise SystemExit("Use --list or provide --galaxy NAME")

    g = load_galaxy_rc(args.galaxy)

    print(f"[Galaxy] {g.name}")
    print(f"[Galaxy] points = {len(g.R_kpc)} | Rmax = {g.R_kpc.max():.6g} kpc")
    print(f"[Galaxy] first row: R={g.R_kpc[0]:.6g}  Vb={g.V_baryon_kms[0]:.6g}  Vt={g.V_total_kms[0]:.6g}")

    # Rd_star is optional; only works if galaxies.csv is present at the expected path
    try:
        Rd = get_Rd_star_kpc(args.galaxy)
        print(f"[Galaxy] Rd_star = {Rd:.6g} kpc")
    except Exception as e:
        print(f"[Galaxy] Rd_star unavailable: {e}")


if __name__ == "__main__":
    main()
