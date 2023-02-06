library(pacman)
p_load(data.table, tidyverse, magrittr, mpra, pheatmap)

# load the barcode_stats_norm data
barcode_stats <- fread("barcode_stats.tsv")


guides <- fread("oligo_guides.tsv")


promoters <- unique(barcode_stats$promoter)


promoters_interest <- c("j23119", "lacUV5", "acrAB", "bamE", "dapA", "hdeA", "letA", "lolB", "lptD", "marR", "micA", "mraZ", "pspA", "rcsA", "rpoE", "tolC", "waaQ", "ygiM") 


promoters_interest <- promoters_interest %>% factor(levels = promoters_interest)


spacer_stats_spread <- barcode_stats %>% 
	filter(promoter %in% promoters_interest) %>%
	mutate(promoter = factor(promoter, levels = promoters_interest)) %>%
	melt(id.vars = c("promoter", "spacer", "barcode", "nucleotide"), variable.name = "batch", value.name = "count") %>% 
	group_by(promoter, spacer, nucleotide, batch) %>% summarise(count = sum(count)) %>% 
	data.table %>% dcast(spacer + nucleotide ~ batch + promoter, value.var = "count") 


spread_data_columns <- spacer_stats_spread %>% 
	colnames %>% `[`(spacer_stats_spread %>% colnames %>% grepl("^T", .))


design_promoters <- spread_data_columns %>% stringr::str_extract("[^_]+$")
design_promoters <- gsub("j23119", "constitutive", design_promoters)
design_promoters <- gsub("lacUV5", "constitutive", design_promoters)


design_timing    <- spread_data_columns %>% stringr::str_extract("^T[02]")


spread_data_design <- model.matrix( ~ design_timing + design_promoters) %>% 
	set_colnames(c("intercept", (design_timing %>% unique)[-1], (design_promoters %>% unique)[-1])) %>%
	set_rownames(spread_data_columns)

spread_data_block <- spread_data_columns %>% stringr::str_extract("(A|B|C|D)_[0-9]")

spread_data_MPRAset <- MPRASet(
	RNA = filter(spacer_stats_spread, nucleotide == "RNA") %>% select(all_of(spread_data_columns)),
	DNA = filter(spacer_stats_spread, nucleotide == "DNA") %>% select(all_of(spread_data_columns)),
	eid = filter(spacer_stats_spread, nucleotide == "RNA") %>% mutate(spacer = as.character(spacer)) %>% pull(spacer))

spread_data_mpralm <- mpralm(
	spread_data_MPRAset,
	design = spread_data_design,
	block = spread_data_block,
	model_type = "corr_groups",
	aggregate = "none",
	plot = FALSE,
	normalize = TRUE)


### these are some data manipulation to get the results into a decent format

results <- spread_data_mpralm %>% topTable(n = Inf) %>% data.table(keep.rownames = "spacer") %>% select(spacer, T2, all_of(promoters_interest[-1:-2])) %>% inner_join(guides)

reshaped <- results %>% melt(id.vars = c("type", "gene", "locus_tag", "offset", "spacer"), variable.name = "promoter", value.name = "logFC") %>% group_by(type, gene, promoter) %>% summarise(mlogFC = mean(logFC)) %>% filter(type != "controls") %>% filter(type == "perfect") %>% ungroup %>% data.table %>% dcast(gene ~ promoter, value.var = "mlogFC")

spread_data_mpralm %>% topTable(n = Inf, coef = "pspA") %>% data.table(keep.rownames = "spacer") %>% inner_join(guides) %>% group_by(gene, locus_tag, type) %>% summarise(mLFC = mean(logFC)) %>% ungroup() %>% filter(type %like% "perfect") %>% select(locus_tag, mLFC) %>% clipr::write_clip()

###############################

results_colnames <- spread_data_mpralm %>% topTable() %>% data.table() %>% colnames


##########################################################################################

# Initialize data.table objects to store results
results_FDR <- guides %>% select(spacer)
results_LFC <- guides %>% select(spacer)

##########################################################################################

rm(overall_results)
predictors = c("intercept", "T2", design_promoters %>% unique %>% `[`(-1))

for (i in predictors) {
	
	# Perform a generalized linear model test using the contrast in aba_contrast[,i]
	one_result <- spread_data_mpralm %>% topTable(coef = i, n = Inf) %>% data.table(keep.rownames = "spacer")

	# Print a message indicating which contrast is being processed
	print(paste("Processing results for", i, "..."))

	if (!exists("overall_results")) {
		overall_results <- one_result %>% 
			select(spacer, logFC, adj.P.Val) %>% 
			mutate(predictor = i) }

	else {
		overall_results <- overall_results %>% rbind(one_result %>% 
			select(spacer, logFC, adj.P.Val) %>% 
			mutate(predictor = i)) }

}

overall_results_summary <- overall_results %>% 
	inner_join(guides) %>% 
	group_by(gene, locus_tag, type, predictor) %>% 
	summarise(mLFC = mean(logFC), FDR = poolr::stouffer(adj.P.Val)$p) %>%
	mutate(predictor = factor(predictor, levels = predictors))

overall_results_summary$predictor <- factor(overall_results_summary$predictor, levels = predictors)

overall_results_mLFC <- overall_results_summary %>% data.table %>% dcast(gene + type ~ predictor, value.var = "mLFC")
overall_results_FDR <- overall_results_summary %>% data.table %>% dcast(gene + type ~ predictor, value.var = "FDR")

