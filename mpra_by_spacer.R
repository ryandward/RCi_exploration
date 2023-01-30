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
data_columns <- barcode_stats_norm %>% 
	colnames %>% `[`(barcode_stats_norm %>% colnames %>% grepl("^T", .))

spacer_design <- data.frame(
	intcpt = 1, 
	T2 = grepl("T2", data_columns)) %>%
	set_rownames(data_columns)

spacer_block <- data_columns %>% stringr::str_extract("(A|B|C|D)_[0-9]")


# for each spacer
spacer_results <- barcode_stats %>%
	filter(spacer %in% (guides %>% filter(type == "controls") %>% pull(spacer))) %>%
	group_by(spacer) %>% nest %>% 
	mutate(
		# create a new MPRASet with selected columns
		mpraset = map(data, ~MPRASet(
			RNA = .x %>% filter(nucleotide == "RNA") %>% select(all_of(data_columns)),
			DNA = .x %>% filter(nucleotide == "DNA") %>% select(all_of(data_columns)),
			eid = .x %>% filter(nucleotide == "RNA") %>% pull(promoter),
			barcode = .x %>% filter(nucleotide == "RNA") %>% pull(barcode)))) %>% 
	mutate(
		mpralm = map(mpraset, ~mpralm(
			.x,
			design = spacer_design,
			block = spacer_block,
			model_type = "corr_groups",
			aggregate = "sum",
			plot = FALSE,
			normalize = TRUE)),
		toptable = map(
			mpralm[1], ~topTable(.x, number = Inf, coef = 2) %>% 
				data.table(keep.rownames = "promoter"))) %>%
	select(spacer, toptable) %>%
	unnest(cols = c(toptable))

# for each promoter
promoter_results <- barcode_stats %>%
	filter(spacer %in% (guides %>% pull(spacer))) %>%
	group_by(promoter) %>% nest %>% 
	mutate(
		# create a new MPRASet with selected columns
		mpraset = map(data, ~MPRASet(
			RNA = .x %>% filter(nucleotide == "RNA") %>% select(all_of(data_columns)),
			DNA = .x %>% filter(nucleotide == "DNA") %>% select(all_of(data_columns)),
			eid = .x %>% filter(nucleotide == "RNA") %>% pull(spacer),
			barcode = .x %>% filter(nucleotide == "RNA") %>% pull(barcode)))) %>% 
	mutate(
		mpralm = map(mpraset, ~mpralm(
			.x,
			design = spacer_design,
			block = spacer_block,
			model_type = "corr_groups",
			aggregate = "sum",
			plot = FALSE,
			normalize = TRUE)),
		toptable = map(
			mpralm[1], ~topTable(.x, number = Inf, coef = 2) %>% 
				data.table(keep.rownames = "spacer"))) %>%
	select(promoter, toptable) %>%
	unnest(cols = c(toptable)) %>% 
	inner_join(guides)

promoter_results %>% ggplot(aes(x = AveExpr, y = -log10(adj.P.Val))) + geom_point(aes(colour = promoter))

