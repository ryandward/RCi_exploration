# # library imports
library(data.table)
library(tidyverse)
library(magrittr)
library(edgeR)
library(poolr)
library(lme4)

setDTthreads(12)

##########################################################################################

lane_spread_data <- fread("data_products/lane_spread_data.tsv.gz")

##########################################################################################

# get rid of barcodes and sum across promoter x guide combos
# data has now lost a level of complexity in favor of comparability across promoters
# data size is also smaller.


barcodes <- fread(
	"data_products//barcode.winners.tsv.gz", 
	header = F, 
	col.names = c("barcode", "promoter", "spacer", "count", "all_counts", "family_size")) %>% 
	select(barcode, promoter, spacer)

lane_spread_data <- barcodes[lane_spread_data, on = .(barcode)]

# get rid of barcodes and sum across promoter x guide combos
# data has now lost a level of complexity in favor of comparability across promoters
# data size is also smaller.

lane_spread_spacers <- lane_spread_data %>% 
	melt(id.vars = c("nucleotide", "barcode", "promoter", "spacer", "timing"), value.name = "count", variable.name = "replicate")
	
lane_spread_spacers <- lane_spread_data %>% 	
	group_by(promoter, spacer, nucleotide, timing, replicate) %>%
	summarise(count = sum(count, na.rm = TRUE)) %>%
	ungroup 

lane_spread_spacers <- lane_spread_spacers %>% 
pivot_wider(
		id_cols = c("spacer"), 
		values_from = "count", 
		names_from = c("promoter", "nucleotide", "timing", "replicate"),
		values_fill = 0)

lane_spread_spacers <- lane_spread_spacers %>% 
	data.table %>% 
	melt(id.vars = "spacer", value.name = "count") %>% 
	separate(variable, c("promoter", "nucleotide", "timing", "replicate", "lane"), sep = "_") %>% 
	pivot_wider(id_cols = c("promoter", "timing", "replicate", "lane", "nucleotide"), names_from = "spacer", values_from = "count") 

no_RNA_spacers <- lane_spread_spacers %>% 
	filter(nucleotide == "RNA") %>% 
	select(-promoter, -timing, -replicate, -lane, -nucleotide) %>%
	data.matrix() %>% t %>% rowSums() %>% enframe %>% filter(value == 0) %>% pull(name)


no_DNA_spacers <- lane_spread_spacers %>% 
	filter(nucleotide == "DNA") %>% 
	select(-promoter, -timing, -replicate, -lane, -nucleotide) %>%
	data.matrix() %>% t %>% rowSums() %>% enframe %>% filter(value == 0) %>% pull(name)

incomplete_spacers <- c(no_RNA_spacers, no_DNA_spacers)

columns_to_keep <- setdiff(colnames(lane_spread_spacers), incomplete_spacers)

lane_spread_spacers <- lane_spread_spacers[, columns_to_keep, with = FALSE]

lane_spread_spacers <- lane_spread_spacers %>% arrange(promoter, nucleotide, replicate, lane)

fwrite(lane_spread_spacers, "lane_spread_spacers.tsv.gz", sep = "\t")

##########################################################################################

