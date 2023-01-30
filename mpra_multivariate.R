spacer_stats_spread <- barcode_stats %>% 
	filter(promoter %in% c("lacUV5", "j23119")) %>%
	melt(id.vars = c("promoter", "spacer", "barcode", "nucleotide"), variable.name = "batch", value.name = "count") %>% 
	group_by(promoter, spacer, nucleotide, batch) %>% summarise(count = sum(count)) %>% 
	setDT %>% dcast(spacer + nucleotide ~ batch + promoter, value.var = "count", fill = 0) %>%
	arrange("promoter", "nucleotide")

spread_data_columns <- spacer_stats_spread %>% 
	colnames %>% `[`(spacer_stats_spread %>% colnames %>% grepl("^T", .))

spread_data_design <- data.frame(
	intcpt = 1, 
	T2 = grepl("T2", spread_data_columns),
	lacUV5 = grepl("lacUV5", spread_data_columns),
	acrAB = grepl("acrAB", spread_data_columns),
	lptD = grepl("lptD", spread_data_columns)) %>%
	set_rownames(spread_data_columns)

spread_data_design <- 
	model.matrix(~  spread_data_design, data = spacer_stats_spread %>% data.table) %>%
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
	model_type = "indep_groups",
	aggregate = "none",
	plot = FALSE,
	normalize = TRUE)

