#!/usr/bin/env python3
import argparse
import gzip
import json
import sys
import math
from scipy.stats import norm

def compute_se(beta, p_value):
    """Computes standard error from BETA and P-value."""
    try:
        p = float(p_value)
        b = float(beta)
        if p <= 0 or p >= 1: return "NA"
        
        # isf is the Inverse Survival Function
        z = abs(norm.isf(p / 2)) 
        
        if z == 0: return "NA"
        se = abs(b / z)
        return f"{se:.6g}"
    except (ValueError, TypeError):
        return "NA"

def compute_p(beta, se):
    """Computes two-tailed P-value from BETA and SE."""
    try:
        b = float(beta)
        s = float(se)
        if s == 0: return "NA"
        
        z = abs(b / s)
        # Survival function (1 - cdf) * 2 for two-tailed p-value
        p = 2 * norm.sf(z)
        return f"{p:.6g}"
    except (ValueError, TypeError):
        return "NA"

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', required=True)
    parser.add_argument('--meta', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()

    with open(args.meta, 'r') as f:
        meta = json.load(f)

    delim = meta.get("original_delimiter", "\t")
    rename_map = meta.get("rename_map", {})
    needs_cpra = meta.get("needs_cpra_split", False)
    cpra_col_name = meta.get("cpra_column", None)
    needs_beta = meta.get("needs_beta", False)
    needs_se = meta.get("needs_se", False)
    needs_p = meta.get("needs_p", False)
    needs_snpid = meta.get("needs_snpid", False)

    open_in = gzip.open if args.input.endswith('.gz') else open
    open_out = gzip.open if args.output.endswith('.gz') else open

    with open_in(args.input, 'rt') as f_in, open_out(args.output, 'wt') as f_out:
        orig_headers = f_in.readline().strip().split(delim)
        
        internal_headers = []
        for h in orig_headers:
            clean_h = h.strip().upper()
            mapped_h = rename_map.get(clean_h, clean_h)
            internal_headers.append(mapped_h)
        
        CORE_ORDER = ['CHR', 'BP', 'SNPID', 'REF', 'ALT', 'BETA', 'SE', 'P']
        
        cols_to_drop = ['OR', 'Z']
        if needs_cpra and cpra_col_name:
            cols_to_drop.append(rename_map.get(cpra_col_name.upper(), cpra_col_name.upper()))

        extra_cols = [h for h in internal_headers if h not in CORE_ORDER and h not in cols_to_drop]
        final_headers = CORE_ORDER + extra_cols
        
        f_out.write('\t'.join(final_headers) + '\n')

        idx_map = {h: i for i, h in enumerate(internal_headers)}
        cpra_idx = orig_headers.index(cpra_col_name) if (needs_cpra and cpra_col_name in orig_headers) else -1

        for line in f_in:
            if not line.strip(): continue
            parts = line.strip().split(delim)
            
            row_data = {}
            
            # 1. Base Population
            for h in internal_headers:
                if h in idx_map and idx_map[h] < len(parts):
                    row_data[h] = parts[idx_map[h]].strip()
                else:
                    row_data[h] = "NA"

            # 2. Split CPRA if requested
            if needs_cpra and cpra_idx != -1 and cpra_idx < len(parts):
                cpra_val = parts[cpra_idx]
                sub_parts = cpra_val.split(':')
                if len(sub_parts) >= 4:
                    row_data['CHR'] = sub_parts[0]
                    row_data['BP'] = sub_parts[1]
                    row_data['REF'] = sub_parts[2]
                    row_data['ALT'] = sub_parts[3]

            # 3. Clean CHR values 
            if 'CHR' in row_data and row_data['CHR'] != "NA":
                clean_chr = row_data['CHR'].upper().replace('CHR', '')
                if clean_chr == '23': clean_chr = 'X'
                if clean_chr == '24': clean_chr = 'Y'
                if clean_chr == '26': clean_chr = 'MT'
                row_data['CHR'] = clean_chr

            # 4. Construct SNPID if missing
            if needs_snpid and row_data.get('CHR', 'NA') != "NA" and row_data.get('BP', 'NA') != "NA":
                row_data['SNPID'] = f"{row_data['CHR']}:{row_data['BP']}"

            # 5. Compute BETA from OR
            if needs_beta and 'OR' in row_data and row_data['OR'] not in ['NA', '']:
                try:
                    row_data['BETA'] = str(math.log(float(row_data['OR'])))
                except ValueError:
                    row_data['BETA'] = "NA"

            # 6. Compute SE from P and BETA
            if needs_se and row_data.get('BETA', 'NA') != "NA" and row_data.get('P', 'NA') != "NA":
                row_data['SE'] = compute_se(row_data['BETA'], row_data['P'])

            # 7. Compute P from BETA and SE
            if needs_p and row_data.get('BETA', 'NA') != "NA" and row_data.get('SE', 'NA') != "NA":
                row_data['P'] = compute_p(row_data['BETA'], row_data['SE'])

            # Build final output string
            final_row = [row_data.get(h, "NA") for h in final_headers]
            f_out.write('\t'.join(final_row) + '\n')

if __name__ == "__main__":
    main()