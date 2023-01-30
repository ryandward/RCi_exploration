library(pacman)

p_load(data.table, tidyverse, mpra, magrittr)

# load the barcode_stats_norm data
barcode_stats_norm <- fread("barcode_stats_norm.tsv")

# load the barcode_stats_norm data
barcode_stats <- fread("barcode_stats.tsv")

# load the guides data
guides <- fread("oligo_guides.tsv")

# identify unique essential guides
essentials <- guides %>% filter(type == "mismatch") %>% pull(locus_tag) %>% unique

#rename guide types
guides[locus_tag %in% essentials & type == "perfect", type := "perfect_essential"]

# select columns containing T
data_columns <- barcode_stats %>% 
	colnames %>% `[`(barcode_stats %>% colnames %>% grepl("^T", .))

spacer_design <- data.frame(
	intcpt = 1, 
	T2 = grepl("T2", data_columns)) %>%
	set_rownames(data_columns)

spacer_block <- data_columns %>% stringr::str_extract("(A|B|C|D)_[0-9]")

control_barcodes <- barcode_stats 
# %>%
# 	filter(spacer %in% (guides %>% filter(type == "perfect") %>% pull(spacer)))

# for each spacer
control_MPRAset <- MPRASet(
	RNA = filter(control_barcodes, nucleotide == "RNA") %>% select(all_of(data_columns)),
	DNA = filter(control_barcodes, nucleotide == "DNA") %>% select(all_of(data_columns)),
	eid = filter(control_barcodes, nucleotide == "RNA") %>% mutate(eid = paste(promoter, spacer, sep = "_")) %>% pull(eid),
	barcode = filter(control_barcodes, nucleotide == "RNA") %>% pull(barcode))

control_mpralm <- mpralm(
	control_MPRAset,
	design = spacer_design,
	block = spacer_block,
	model_type = "corr_groups",
	aggregate = "sum",
	plot = FALSE,
	normalize = TRUE)

control_results <- control_mpralm %>% 
	topTable(coef = 2, number = Inf) %>% 
	data.table(keep.rownames = "eid") %>% 
	separate(eid, c("promoter", "spacer"), sep = "_") %>% 
	inner_join(guides)

# control_results %>% ggplot(aes(x = logFC, y = -log10(adj.P.Val))) + geom_point(aes(colour = promoter))

control_results[, logFC.adj := logFC - mean(logFC[promoter == "j23119"]), by = spacer]

