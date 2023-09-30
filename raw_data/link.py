import gzip
from collections import Counter, defaultdict
from rich.console import Console

console = Console(stderr=True)

promoters = {line.strip() for line in open('promoters.fasta') if not line.startswith('>')}

oligos = {line.strip() for line in open('oligo_guides.fasta') if not line.startswith('>')}

counts = Counter()
barcode_to_promoter_oligo = defaultdict(Counter)
filtered_counts = 0

with console.status("Processing...", spinner="dots") as status:
    for lane in ['L001', 'L002']:
        with gzip.open(f'Rci_lib_0421_S1_{lane}_R1_001.fastq.gz', 'rt') as f1, \
            gzip.open(f'Rci_lib_0421_S1_{lane}_R2_001.fastq.gz', 'rt') as f2, \
            gzip.open(f'Rci_lib_0421_S1_{lane}_R3_001.fastq.gz', 'rt') as f3:

            for i, (line1, line2, line3) in enumerate(zip(f1, f2, f3)):
                if i % 4 == 1:
                    stripped_line1, stripped_line2, stripped_line3 = line1.strip(), line2.strip(), line3.strip()
                    if stripped_line1 in promoters and stripped_line3 in oligos:
                        counts[(stripped_line1, stripped_line2, stripped_line3)] += 1
                        barcode_to_promoter_oligo[stripped_line2][(stripped_line1, stripped_line3)] += 1
                        filtered_counts += 1

                        if filtered_counts % 10000 == 0:
                            status.update(f"{filtered_counts:,} reads were found in defined promoters and oligos...")

# Filtering barcodes based on the 90% threshold and calculating proportion
filtered_data = {}
for barcode, promoter_oligo_counts in barcode_to_promoter_oligo.items():
    total = sum(promoter_oligo_counts.values())
    for (promoter, oligo), count in promoter_oligo_counts.items():
        proportion = count / total
        if proportion >= 0.9:
            filtered_data[(promoter, barcode, oligo)] = (count, proportion)

# Print filtered counts and proportions
for (promoter, barcode, oligo), (count, proportion) in filtered_data.items():
    print(f'{promoter}\t{barcode}\t{oligo}\t{count}\t{proportion:.4f}')
