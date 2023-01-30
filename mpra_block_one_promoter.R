# Read in gene data
gene_data <- fread("gene_to_b.tsv")

# read raw data
lane_spread_spacers <- fread("lane_spread_spacers.tsv.gz")



# Reshape data
promoter_data <- lane_spread_spacers %>% 
	melt(
		id.vars = c("promoter", "timing", "replicate", "nucleotide", "lane"),
		variable.name = "spacer", value.name = "count") %>%
	# mutate(count = as.numeric(count)) %>%
	# mutate(random_offset = runif(nrow(.), min = 0, max = 20)) %>%
	# mutate(count = count + random_offset) %>%
	# mutate(count = case_when(count == 0 ~ NA_real_, TRUE ~ count)) %>%
	filter(promoter == "lolB") %>%
	dcast(spacer + nucleotide ~ promoter + timing + replicate + lane, value.var = "count")

# Create block set
promoter_MPRASet <- MPRASet(
	RNA = promoter_data %>% filter(nucleotide == "RNA") %>% select(-spacer, -nucleotide),
	DNA = promoter_data %>% filter(nucleotide == "DNA") %>% select(-spacer, -nucleotide),
	eid = promoter_data %>% filter(nucleotide == "RNA") %>% pull(spacer) %>% as.character 
)

promoter_design <- data.frame(
	intcpt = 1, 
	T2 = grepl("T2", colnames(promoter_MPRASet))
) %>% 
	set_rownames(colnames(promoter_MPRASet))

promoter_block <- promoter_data %>% 
	colnames %>% `[`(promoter_data %>% colnames %>% grepl("_", .)) %>% 
	stringr::str_extract("(A|B|C|D)_[0-9]")

promoter_mpralm <- mpralm(
	promoter_MPRASet, 
	design = promoter_design,
	block = promoter_block,
	normalize = TRUE,
	model_type = "indep_groups", 
	aggregate = "none")

# Define function to handle missing values
handle_NA <- function(x) {
	if (!any(is.na(x))) {
		result <- poolr::stouffer(x)
		p_value <- result$p
	} else {
		p_value <- NA
	}
	return(p_value)
}

##########################################################################################

promoter_guide_results_raw <- topTable(promoter_mpralm, coef = 2, number = Inf) %>% 
	data.table(keep.rownames = "spacer") %>% 
	inner_join(guides) 

promoter_gene_results_raw <- promoter_guide_results_raw %>% 
	group_by(locus_tag, type) %>% 
	summarise(
		logFC = mean(logFC), 
		FDR = handle_NA(adj.P.Val))

promoter_gene_results_raw %>% ggplot(aes(x = logFC, y = -log10(FDR))) + geom_point()

