#!/usr/bin/env python3
import argparse
import gzip
import subprocess
import os
import sys

def main():
    parser = argparse.ArgumentParser(description="Wrapper to securely liftover GWAS files.")
    parser.add_argument('--input', required=True)
    parser.add_argument('--chain', required=True)
    parser.add_argument('--liftover', default='liftOver')
    parser.add_argument('--output', required=True)
    parser.add_argument('--log', required=True) # New required argument
    args = parser.parse_args()

    open_in = gzip.open if args.input.endswith('.gz') else open
    open_out = gzip.open if args.output.endswith('.gz') else open
    
    bed_file = "temp.bed"
    mapped_bed = "mapped.bed"
    unmapped_bed = "unmapped.bed"

    # 1. Create a safe BED file using row indices instead of SNPIDs
    with open_in(args.input, 'rt') as f_in, open(bed_file, 'w') as f_bed:
        headers = f_in.readline().strip().split('\t')
        chr_idx = headers.index('CHR')
        bp_idx = headers.index('BP')
        
        row_idx = 0
        for line in f_in:
            if not line.strip(): continue
            parts = line.strip().split('\t')
            chrom = parts[chr_idx]
            bp = parts[bp_idx]
            
            if chrom != "NA" and bp != "NA":
                chrom_out = chrom if chrom.startswith('chr') else f"chr{chrom}"
                try:
                    pos = int(bp)
                    f_bed.write(f"{chrom_out}\t{pos-1}\t{pos}\t{row_idx}\n")
                except ValueError:
                    pass
            row_idx += 1

    total_rows_processed = row_idx

    # 2. Execute UCSC liftOver
    cmd = [args.liftover, bed_file, args.chain, mapped_bed, unmapped_bed]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"liftOver failed: {result.stderr}", file=sys.stderr)
        sys.exit(1)

    # 3. Load successfully mapped coordinates into memory
    mapped_coords = {}
    with open(mapped_bed, 'r') as f_map:
        for line in f_map:
            parts = line.strip().split('\t')
            new_chr = parts[0].replace('chr', '')
            new_bp = parts[2]
            r_idx = int(parts[3])
            mapped_coords[r_idx] = (new_chr, new_bp)

    # 4. Stream original file, replacing CHR/BP for mapped rows, dropping unmapped
    with open_in(args.input, 'rt') as f_in, open_out(args.output, 'wt') as f_out:
        header_line = f_in.readline()
        f_out.write(header_line)
        
        row_idx = 0
        for line in f_in:
            if not line.strip(): continue
            if row_idx in mapped_coords:
                parts = line.strip().split('\t')
                new_chr, new_bp = mapped_coords[row_idx]
                
                parts[chr_idx] = new_chr
                parts[bp_idx] = new_bp
                
                f_out.write('\t'.join(parts) + '\n')
            row_idx += 1

    # 5. Calculate dropped SNPs and write to log
    dropped_count = total_rows_processed - len(mapped_coords)
    with open(args.log, 'w') as f_log:
        f_log.write(f"{dropped_count} SNPs could not be mapped to hg38 and were removed.")

    # Clean up temporary files
    os.remove(bed_file)
    os.remove(mapped_bed)
    os.remove(unmapped_bed)

if __name__ == "__main__":
    main()