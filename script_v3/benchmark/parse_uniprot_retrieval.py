import sys
import re

# Configure which column has the PDB info (0-based index)
# 8th column = index 7
TARGET_COL = 7 

def process_file(input_file, output_file=None):
    output=[]
    with open(input_file, 'r') as f:
        lines = f.readlines()
        for idx in range(1, len(lines)):
            line = lines[idx].strip()
            if not line: continue
            
            # Split line into columns (assuming tab-separated)
            cols = line.split('\t')
            uniprot_id = cols[0]
            review_status = cols[1]
            seq_length = int(cols[6])
            
            # Safety check: ensure the line actually has 8 columns
            if len(cols) <= TARGET_COL:
                print(line) # Print as-is if column missing
                continue

            # Get the PDB blob (8th column)
            pdb_blob = cols[TARGET_COL]
            
            # Regex explanation:
            # "      : Match a literal quote
            # (.*?)  : Capture everything inside non-greedily
            # "      : Match the closing literal quote
            pdb_entries = re.findall(r'"(.*?)"', pdb_blob)

            # If no PDB entries found, print original line
            if not pdb_entries:
                print(line)
                continue

            # Print a new line for each PDB entry
            for entry in pdb_entries:
                entry_items = [i.strip() for i in entry.split(";")]
                pdb_id = entry_items[0]
                chain_info = entry_items[3] if len(entry_items) > 3 else ""
                chain_info = chain_info.rstrip('.')
                chain_infos = chain_info.split(',') if chain_info else []

                for info in chain_infos:
                    info = info.strip()
                    if '=' not in info:
                        continue
                    
                    chain_str, range_str = info.split('=')
                    chains = chain_str.split('/')
                    cov_len = 0
                    coverage = 0.0
                    if len(range_str) > 3:
                        cov_len = int(range_str.split('-')[1]) - int(range_str.split('-')[0]) + 1
                        coverage = round(cov_len / seq_length,3)
                    for chain in chains:
                        chain_line=";".join([chain, range_str, str(cov_len), f"{coverage:.3f}"])
                        output.append((uniprot_id, review_status, seq_length, pdb_id, chain_line))
                        # print('\t'.join([uniprot_id, review_status, str(seq_length), pdb_id, chain, range_str, str(cov_len), f"{cov_len / seq_length:.3f}"]))

    if output_file:
        with open(output_file, 'w') as out_f:
            for record in output:
                out_f.write('\t'.join(map(str, record)) + '\n')


if __name__ == "__main__":
    if len(sys.argv) > 2:
        process_file(sys.argv[1], sys.argv[2])
    else:
        print("Usage: python parse_uniprot_retrieval.py input.txt output.txt")