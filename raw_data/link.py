import gzip
import json
import math
import time
from collections import Counter, defaultdict

import numpy as np
from rich.console import Console
from sklearn.cluster import KMeans

# Record the start time
start_time = time.time()

console = Console(stderr=True, highlight=True)

promoters = {line.strip() for line in open('promoters.fasta') if not line.startswith('>')}

oligos = {line.strip() for line in open('oligo_guides.fasta') if not line.startswith('>')}

import gzip
import json
import time
from collections import Counter, defaultdict

import numpy as np

# Assuming console and start_time are defined
# Assuming promoters and oligos are defined

counts = defaultdict(int)
barcode_to_promoter_oligo = defaultdict(lambda: defaultdict(int))

filtered_counts = 0
dark_promoters = 0
dark_oligos = 0
dark_barcodes = 0
oligo_offcount = 0
read_count = 0  # Replacing i
sgRNA_error_rate = 0

def calculate_min_count(barcode_to_promoter_oligo):
    # Extract barcode frequencies and reshape for K-means
    barcode_frequencies = np.array([sum(counts.values()) for counts in barcode_to_promoter_oligo.values()])
    barcode_frequencies_reshaped = barcode_frequencies.reshape(-1, 1)

    # Perform K-means clustering
    kmeans = KMeans(n_clusters=2, random_state=0, n_init = 1).fit(barcode_frequencies_reshaped)
    cluster_centers = sorted(kmeans.cluster_centers_.flatten())

    # Calculate threshold N
    threshold_N = np.mean(cluster_centers)
    
    return threshold_N


with console.status("[bold]Reading and filtering...[/bold]", spinner="dots") as status:
    for lane in ['L001']:
        with gzip.open(f'Undetermined_S0_{lane}_R1_001.fastq.gz', 'rt') as f1, \
             gzip.open(f'Undetermined_S0_{lane}_R2_001.fastq.gz', 'rt') as f2, \
             gzip.open(f'Undetermined_S0_{lane}_R3_001.fastq.gz', 'rt') as f3:

            for line1, line2, line3 in zip(f1, f2, f3):

                if read_count % 4 == 1:
                    stripped_line1, stripped_line2, stripped_line3 = line1.strip(), line2.strip(), line3.strip()

                    is_dark_promoter = all(c == 'G' for c in stripped_line1)
                    is_dark_barcode = all(c == 'G' for c in stripped_line2)
                    is_dark_oligo = all(c == 'G' for c in stripped_line3)

                    dark_promoters += is_dark_promoter
                    dark_barcodes += is_dark_barcode
                    dark_oligos += is_dark_oligo
                    
                    if not (is_dark_promoter or is_dark_barcode or is_dark_oligo):
                        if stripped_line3 not in oligos:
                            oligo_offcount += 1
                        elif stripped_line1 in promoters and stripped_line3 in oligos:
                            counts[(stripped_line1, stripped_line2, stripped_line3)] += 1
                            barcode_to_promoter_oligo[stripped_line2][(stripped_line1, stripped_line3)] += 1
                            filtered_counts += 1

                        if read_count // 4 % 100000 == 0 and read_count // 4 > 0:
                            sgRNA_error_rate = round(oligo_offcount / (read_count // 4 - dark_oligos), 4) if filtered_counts > 0 else 0
                            seconds_elapsed = round(time.time() - start_time, 2)
                            data = {
                                "processed_reads": read_count // 4,
                                "read_pass": filtered_counts,
                                "sgRNA_error_rate": sgRNA_error_rate,
                                "unique_barcodes": len(barcode_to_promoter_oligo),
                                "oligo_offcount": oligo_offcount,
                                "dark_barcodes": dark_barcodes,
                                "dark_oligos": dark_oligos,
                                "dark_promoters": dark_promoters,
                                "seconds_elapsed": seconds_elapsed,
                                "lane": lane,
                            }
                            console.print(json.dumps(data, indent=4))
                read_count += 1

pass_count = 0

# min_count = calculate_min_count(barcode_to_promoter_oligo)
min_count = 0

# Filtering logic
filtered_data = {}
console.log(f"[bold]Identifying best barcodes with at least {int(round(0, 0))} counts and 95% unique...[/bold]")

for barcode, promoter_oligo_counts in barcode_to_promoter_oligo.items():
    total = sum(promoter_oligo_counts.values())
    for (promoter, oligo), count in promoter_oligo_counts.items():
        proportion = count / total
        if count >= min_count and proportion >= 0.95 :  
            pass_count += 1
            filtered_data[(promoter, barcode, oligo)] = (count, proportion)



# Specify the output file name
output_file = "joined.tsv"
written_lines = 0
console.log(f"[bold]{pass_count:,} ({round(100 * pass_count/ len(barcode_to_promoter_oligo), 2)}%) barcodes passed the filter. Writing to {output_file}.[/bold]")

# Open the output file for writing
with open(output_file, "w") as outfile:
    # Write header to the file
    outfile.write("promoter\tbarcode\toligo\tcount\tproportion\n")
    
    # Iterate through filtered_data and write each entry to the file
    for (promoter, barcode, oligo), (count, proportion) in filtered_data.items():
        # Create a tab-separated line and write it to the file
        line = f'{promoter}\t{barcode}\t{oligo}\t{count}\t{proportion:.4f}\n'
        outfile.write(line)
        written_lines += 1

        if written_lines % 25000 == 0:
            console.log(f"{written_lines:,} lines written to file.")

