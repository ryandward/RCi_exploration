# RCi_exploration
# Barcode Mapping and Threshold Determination with FDR

## Overview

1. **Data Collection**: Reads are extracted from gzipped FASTQ files.
2. **Threshold Determination**: A False Discovery Rate (FDR)-based approach is used.
3. **Data Filtering**: Barcodes meeting certain criteria are selected.

## Data Collection

The data collection process uses gzipped FASTQ files for each lane. The reads are filtered in real-time to exclude 'dark' promoters, barcodes, and oligos (those comprising a single nucleotide type, e.g., all 'G'). Counts are incremented for each unique combination of promoter, barcode, and oligo.

```python
for line1, line2, line3 in zip(f1, f2, f3):
    if read_count % 4 == 1:
        # ... (data collection and filtering logic here)
        counts[(stripped_line1, stripped_line2, stripped_line3)] += 1
```

## Threshold Determination Using FDR

We use an FDR-based approach to determine the minimum count \(N\) for considering a barcode as real. The expected number of false positives is calculated based on the sgRNA error rate.

$$
\text{Expected False Positives} = M \times \text{sgRNA error rate}^N
$$

Where:
* $M$ is the total number of unique barcodes
* $N$ is the count threshold for considering a barcode as real

$$
\text{sgRNA Error Rate} = \frac{\text{Number of off-target sgRNA reads}}{\text{Total number of reads} - \text{Number of dark oligos}}
$$


```python
def calculate_min_count(FDR_target, sgRNA_error_rate, barcode_to_promoter_oligo):
    # ... (FDR logic here)
    expected_false_positives = M * (sgRNA_error_rate ** N)
```

## Data Filtering

Barcodes with counts greater than or equal to \(N\) and representing more than 95% of the total counts for that barcode are considered "real" and are filtered into the final dataset.

```python
if count >= min_count and proportion >= 0.95:  
    pass_count += 1
    filtered_data[(promoter, barcode, oligo)] = (count, proportion)
```
