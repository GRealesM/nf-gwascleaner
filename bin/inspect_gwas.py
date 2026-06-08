#!/usr/bin/env python3
import argparse
import gzip
import json
import sys
import re
import os

ALIAS_MAP = {
    "CHR": ["CHROMOSOME", "CHROM", "CHR_ID", "HG18CHR", "#CHROM", "#CHR", "CHRLOC"],
    "BP": ["POS", "BASE_PAIR_LOCATION", "BP(HG19)", "POSITION", "CHR_POSITION", "POSITION(HG19)", "BP_HG19", "COORDINATE"],
    "SNPID": ["ID", "VARIANT_ID", "MARKERNAME", "SNP", "RSID", "RS_ID", "RSIDS", "SNP_NAME", "SNP_ID", "#SNPID", "RS_NUMBER", "DB_SNP_RS_ID/MARKER", "DBSNP_RS_ID", "VARIANT", "ÍD"],
    "REF": ["OTHERALLELE", "REFERENCE_ALLELE", "REF_ALLELE", "OTHER_ALLELE", "A2_OTHER", "NEA", "REF", "ALLELE1"],
    "ALT": ["EFFECT_ALLELE", "EFFECTALLELE", "A1_EFFECT", "RISK_ALLELE", "EA", "ALT", "ALLELE2"],
    "BETA": ["EFFECT", "BETA_SNP_ADD", "EFFECT_ALT", "EFFB", "ALL_INV_VAR_META_BETA"],
    "SE": ["STANDARD_ERROR", "STDERR", "SEBETA_SNP_ADD", "SEBETA", "SE_EFFB", "ALL_INV_VAR_META_SEBETA", "LOG(OR)_SE"],
    "OR": ["ODDS_RATIO", "ODDSRATIO", "OR(A1)", "ORX"],
    "P": ["P_VALUE", "P.VALUE", "PVALUE", "P-VALUE", "PVAL", "P-VAL", "ALL.P.VALUE", "GC-ADJUSTED_P_", "CHI-SQUARED__P", "P1DF", "ALL_INV_VAR_META_P"],
    "LOG10P": ["LOG10P"],
    "-LOG10P": ["_-LOG10_P-VALUE"],
    "ALT_FREQ": ["EFFECT_ALLELE_FREQUENCY", "MAF"],
    "EMP_BETA": ["EMP_BETA"],
    "EMP_P": ["EMP1"],
    "EMP_SE": ["EMP_SE"],
    "hm_ALT": ["HM_EFFECT_ALLELE"],
    "hm_BETA": ["HM_BETA"],
    "hm_BP": ["HM_POS"],
    "hm_CHR": ["HM_CHROM"],
    "hm_OR": ["HM_ODDS_RATIO"],
    "hm_REF": ["HM_OTHER_ALLELE"],
    "hm_SNPID": ["HM_RSID"],
    "N": ["N"],
    "RSQ": ["RSQ"],
    "Z": ["ZSCORE", "Z_STAT"]
}

LOOKUP_DICT = {}
for std_name, aliases in ALIAS_MAP.items():
    LOOKUP_DICT[std_name] = std_name
    for alias in aliases:
        LOOKUP_DICT[alias] = std_name

CPRA_REGEX = re.compile(r'^[0-9XYxyMTmt]+:\d+:[A-Za-z]+:[A-Za-z]+$')
NA_VALUES = {'', 'NA', 'NAN', 'NULL', '.', 'N/A'}

def load_manifest(manifest_path):
    hg18, hg19, hg38 = set(), set(), set()
    try:
        with open(manifest_path, 'r') as f:
            next(f)
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 7:
                    hg18.add(f"{parts[1].replace('chr', '')}:{parts[2]}")
                    hg19.add(f"{parts[3].replace('chr', '')}:{parts[4]}")
                    hg38.add(f"{parts[5].replace('chr', '')}:{parts[6]}")
    except Exception as e:
        print(f"Warning: Could not load manifest: {e}", file=sys.stderr)
    return hg18, hg19, hg38

def sniff_format_and_build(file_path, manifest_path):
    open_func = gzip.open if file_path.endswith('.gz') else open
    hg18_set, hg19_set, hg38_set = load_manifest(manifest_path)
    
    try:
        with open_func(file_path, 'rt') as f:
            header_line = f.readline().strip()
            
            if '\t' in header_line: delim = '\t'
            elif ',' in header_line: delim = ','
            elif ';' in header_line: delim = ';'
            else: delim = ' '
            
            original_headers = header_line.split(delim)
            translated_headers = []
            rename_map = {}

            for orig in original_headers:
                clean_orig = orig.strip().upper()
                std_header = LOOKUP_DICT.get(clean_orig, clean_orig)
                translated_headers.append(std_header)
                if clean_orig != std_header:
                    rename_map[clean_orig] = std_header

            idx_map = {h: i for i, h in enumerate(translated_headers)}
            
            cpra_col_orig_name = None
            build_scores = {"hg18": 0, "hg19": 0, "hg38": 0}
            detected_build = "unknown"
            
            lines_checked = 0
            na_counts = {'P': 0, 'BETA': 0, 'OR': 0, 'SE': 0}
            
            while lines_checked < 500000:
                data_line = f.readline().strip()
                if not data_line: break 
                
                row_values = data_line.split(delim)
                if len(row_values) != len(original_headers):
                    continue

                for col in na_counts.keys():
                    if col in idx_map:
                        val = row_values[idx_map[col]].strip().upper()
                        if val in NA_VALUES:
                            na_counts[col] += 1

                if not cpra_col_orig_name and lines_checked < 5:
                    for i, val in enumerate(row_values):
                        if CPRA_REGEX.match(val.strip()):
                            cpra_col_orig_name = original_headers[i]
                            break

                q_chr, q_bp = None, None
                if 'CHR' in idx_map and 'BP' in idx_map:
                    q_chr = str(row_values[idx_map['CHR']]).upper().replace('CHR', '')
                    q_bp = str(row_values[idx_map['BP']])
                elif cpra_col_orig_name:
                    cpra_idx = original_headers.index(cpra_col_orig_name)
                    parts = row_values[cpra_idx].split(':')
                    if len(parts) >= 2:
                        q_chr, q_bp = parts[0].upper().replace('CHR', ''), parts[1]

                if q_chr and q_bp:
                    coord = f"{q_chr}:{q_bp}"
                    if coord in hg38_set: build_scores["hg38"] += 1
                    if coord in hg19_set: build_scores["hg19"] += 1
                    if coord in hg18_set: build_scores["hg18"] += 1

                if build_scores["hg38"] >= 5: detected_build = "hg38"; break
                if build_scores["hg19"] >= 5: detected_build = "hg19"; break
                if build_scores["hg18"] >= 5: detected_build = "hg18"; break

                lines_checked += 1

            if detected_build == "unknown" and max(build_scores.values()) > 0:
                detected_build = max(build_scores, key=build_scores.get)
                
            na_ratios = {k: (v / lines_checked) if lines_checked > 0 else 0 for k, v in na_counts.items()}

            return translated_headers, original_headers, delim, rename_map, cpra_col_orig_name, detected_build, na_ratios

    except Exception as e:
        print(f"Error reading file: {e}", file=sys.stderr)
        sys.exit(1)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', required=True, help="Path to input GWAS file")
    parser.add_argument('--manifest', required=True, help="Path to Manifest_build_translator.tsv")
    parser.add_argument('--output', default="meta.json", help="Output JSON file name")
    parser.add_argument('--study_id', required=False, help="Identifier")
    args = parser.parse_args()

    study_id = args.study_id
    if not study_id:
        base = os.path.basename(args.input)
        if base.lower().endswith('.gz'): base = base[:-3]
        if base.lower().endswith('.tsv') or base.lower().endswith('.txt') or base.lower().endswith('.csv'): 
            base = base[:-4]
        study_id = base

    translated_headers, orig_headers, original_delim, rename_map, cpra_col, build, na_ratios = sniff_format_and_build(args.input, args.manifest)
    
    meta = {
        "study_id": study_id,
        "fatal_error": False,
        "error_message": "",
        "original_delimiter": original_delim,
        "rename_map": rename_map,
        "has_cpra": bool(cpra_col),
        "needs_cpra_split": False,
        "cpra_column": cpra_col,
        "needs_beta": False,
        "needs_se": False,
        "needs_p": False,
        "needs_snpid": False,
        "build": build
    }

    # Removed P from core_cols
    core_cols = ['CHR', 'BP', 'REF', 'ALT']
    missing_core = [col for col in core_cols if col not in translated_headers]
    
    if meta["has_cpra"] and ('CHR' in missing_core or 'BP' in missing_core):
        meta["needs_cpra_split"] = True
        for c in ['CHR', 'BP', 'REF', 'ALT']:
            if c in missing_core: missing_core.remove(c)

    if missing_core:
        meta["fatal_error"] = True
        meta["error_message"] = f"Missing core columns: {', '.join(missing_core)}"
        
    if 'SNPID' not in translated_headers:
        meta["needs_snpid"] = True
    
    # Logic for missingness and effect size
    if 'BETA' not in translated_headers:
        if 'OR' in translated_headers:
            meta["needs_beta"] = True
            if na_ratios['OR'] > 0.5:
                meta["fatal_error"] = True
                meta["error_message"] = "Over 50% missing data in OR column."
        else:
            meta["fatal_error"] = True
            meta["error_message"] = "Missing BETA and OR."
    else:
        if na_ratios['BETA'] > 0.5:
            meta["fatal_error"] = True
            meta["error_message"] = "Over 50% missing data in BETA column."

    if 'SE' not in translated_headers:
        meta["needs_se"] = True
    elif na_ratios['SE'] > 0.5:
        meta["fatal_error"] = True
        meta["error_message"] = "Over 50% missing data in SE column."
        
    # Check for P
    if 'P' not in translated_headers:
        meta["needs_p"] = True
        # To compute P, we MUST have an effect size and an SE
        has_effect = ('BETA' in translated_headers) or ('OR' in translated_headers)
        has_se = 'SE' in translated_headers
        if not (has_effect and has_se):
            meta["fatal_error"] = True
            meta["error_message"] = "Missing P-value, and lacks BETA/OR + SE required to compute it."
    elif na_ratios['P'] > 0.5:
        meta["fatal_error"] = True
        meta["error_message"] = "Over 50% missing data in P column."

    with open(args.output, 'w') as f:
        json.dump(meta, f, indent=4)

if __name__ == "__main__":
    main()