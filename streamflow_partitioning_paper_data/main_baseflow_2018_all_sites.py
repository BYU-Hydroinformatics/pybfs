#!/usr/bin/env python3
"""
Example usage of BFS for baseflow separation and forecasting

This script demonstrates:
1. Loading streamflow data and site parameters
2. Generating a baseflow table
3. Running BFS for baseflow separation
4. Visualizing results
5. Creating forecasts
"""
import os
import argparse
import pandas as pd
import numpy as np
import pybfs


import os
import argparse
import pandas as pd
import numpy as np
import pybfs


def pick_site_col(params_df: pd.DataFrame) -> str:
    for c in ["tmp.site", "site_no", "SiteID", "siteid", "site_id", "site"]:
        if c in params_df.columns:
            return c
    raise KeyError(
        f"Could not find a site id column in bfs_params file. Found columns: {list(params_df.columns)}"
    )


def main():
    parser = argparse.ArgumentParser(
        description="Run PyBFS for all sites and save BASEFLOW ONLY as a wide CSV (Date + one column per site)."
    )
    parser.add_argument("--params", default="streamflow_partitioning_paper_data/bfs_params_python_all_sites.csv",
                        help="CSV with calibrated parameters for all sites")
    parser.add_argument("--streamflow", default="streamflow_partitioning_paper_data/streamflow_with_date.csv",
                        help="Wide streamflow table: Date column + site columns")
    parser.add_argument("--outdir", default="streamflow_partitioning_paper_data",
                        help="Output directory")
    parser.add_argument("--outfile", default="baseflow_only_2018_all_sites_copy.csv",
                        help="Output wide CSV filename")
    parser.add_argument("--start", default="2018-01-01", help="Start date (inclusive)")
    parser.add_argument("--end", default="2018-12-31", help="End date (inclusive)")
    parser.add_argument("--date_col", default="Date", help="Date column name in streamflow.csv")
    parser.add_argument("--skip_missing", action="store_true",
                        help="Skip sites missing in streamflow columns or failing BFS")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # ---- Load calibrated parameters ----
    bfs_params = pd.read_csv(args.params, encoding="utf-8-sig")
    site_col = pick_site_col(bfs_params)
    site_ids = bfs_params[site_col].astype(str).str.strip().unique().tolist()
    print(f"Loaded {len(site_ids)} sites from: {args.params} (site col='{site_col}')")

    # ---- Load wide streamflow table ----
    qdv = pd.read_csv(args.streamflow, encoding="utf-8-sig")
    if args.date_col not in qdv.columns:
        raise KeyError(f"'{args.date_col}' not found in {args.streamflow}. Columns: {list(qdv.columns)}")

    qdv[args.date_col] = pd.to_datetime(qdv[args.date_col], errors="coerce")
    if qdv[args.date_col].isna().any():
        bad_n = int(qdv[args.date_col].isna().sum())
        raise ValueError(f"{bad_n} Date values failed to parse in {args.streamflow}.")

    # Filter to 2018 window
    start = pd.to_datetime(args.start)
    end = pd.to_datetime(args.end)
    qdv_2018 = qdv[(qdv[args.date_col] >= start) & (qdv[args.date_col] <= end)].copy()

    if len(qdv_2018) == 0:
        raise ValueError(f"No rows remain after filtering {args.start}..{args.end} in {args.streamflow}.")

    # Normalize column names to strings (site ids as headers)
    qdv_2018.columns = qdv_2018.columns.astype(str)
    streamflow_cols = set(qdv_2018.columns)

    # Output frame: Date + baseflow columns (added as we compute them)
    out_df = pd.DataFrame({args.date_col: qdv_2018[args.date_col].values})

    failures = []

    # ---- Run BFS site-by-site and store Baseflow only ----
    for idx, site_id in enumerate(site_ids, start=1):
        site_id_str = str(site_id).strip()

        if site_id_str not in streamflow_cols:
            msg = f"Site {site_id_str} not found as a column in streamflow.csv"
            if args.skip_missing:
                print("SKIP -", msg)
                failures.append({"SiteID": site_id_str, "Reason": "Missing streamflow column"})
                continue
            raise KeyError(msg + " (use --skip_missing to continue)")

        print(f"\n[{idx}/{len(site_ids)}] Computing BASEFLOW for site {site_id_str}")

        try:
            # Build site-specific streamflow df for pybfs
            sf = pd.DataFrame({
                "Date": qdv_2018[args.date_col].values,
                "Streamflow": pd.to_numeric(qdv_2018[site_id_str], errors="coerce").values
            })

            # Get parameters for the site
            try:
                site_int = int(float(site_id_str))
                basin_char, gw_hyd, flow = pybfs.get_values_for_site(bfs_params, site_int)
            except Exception:
                basin_char, gw_hyd, flow = pybfs.get_values_for_site(bfs_params, site_id_str)

            # base_table args: (lb, x1, wb, beta, kb, streamflow_df, por)
            lb = basin_char[1]
            x1 = basin_char[2]
            wb = basin_char[3]
            por = basin_char[4]
            beta = gw_hyd[1]
            kb = gw_hyd[3]

            SBT = pybfs.base_table(lb, x1, wb, beta, kb, sf, por)
            result = pybfs.bfs(sf, SBT, basin_char, gw_hyd, flow)

            # Determine baseflow column name
            if "Baseflow" in result.columns:
                base = result["Baseflow"].values
            elif "Baseflow.L3" in result.columns:
                base = result["Baseflow.L3"].values
            else:
                raise KeyError(f"Baseflow column not found. Result columns: {list(result.columns)}")

            # Add as a new column in wide output
            out_df[site_id_str] = base
            print(f"OK   - {site_id_str}")

        except Exception as e:
            msg = str(e)
            print(f"FAIL - {site_id_str}: {msg}")
            failures.append({"SiteID": site_id_str, "Reason": msg})
            if not args.skip_missing:
                # If you prefer hard fail, remove this and just raise
                pass

    # ---- Save output ----
    out_path = os.path.join(args.outdir, args.outfile)
    out_df.to_csv(out_path, index=False, float_format="%.6g")
    print("\nSaved baseflow-only wide file:", out_path)

    if failures:
        fail_path = os.path.join(args.outdir, "baseflow_2018_failures.csv")
        pd.DataFrame(failures).to_csv(fail_path, index=False)
        print("Saved failures log:", fail_path)


if __name__ == "__main__":
    main()
