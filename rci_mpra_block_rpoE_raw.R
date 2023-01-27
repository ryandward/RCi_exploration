# Read in gene data
gene_data <- fread("gene_to_b.tsv")

# read raw data
lane_spread_spacers <- fread("lane_spread_spacers.tsv.gz")



# Reshape data
rpoE_data <- lane_spread_spacers %>% 
	melt(
		id.vars = c("promoter", "timing", "replicate", "nucleotide", "lane"),
		variable.name = "spacer", value.name = "count") %>%
	# mutate(count = as.numeric(count)) %>%
	# mutate(random_offset = runif(nrow(.), min = 0, max = 20)) %>%
	# mutate(count = count + random_offset) %>%
	# mutate(count = case_when(count == 0 ~ NA_real_, TRUE ~ count)) %>%
	filter(promoter == "rpoE") %>%
	dcast(spacer + nucleotide ~ promoter + timing + replicate + lane, value.var = "count")

# Create block set
rpoE_MPRASet <- MPRASet(
	RNA = rpoE_data %>% filter(nucleotide == "RNA") %>% select(-spacer, -nucleotide),
	DNA = rpoE_data %>% filter(nucleotide == "DNA") %>% select(-spacer, -nucleotide),
	eid = rpoE_data %>% filter(nucleotide == "RNA") %>% pull(spacer) %>% as.character 
)

rpoE_design <- data.frame(
	intcpt = 1, 
	T2 = grepl("T2", colnames(rpoE_MPRASet))
) %>% 
	set_rownames(colnames(rpoE_MPRASet))

rpoE_block <- rpoE_data %>% 
	colnames %>% `[`(rpoE_data %>% colnames %>% grepl("_", .)) %>% 
	stringr::str_extract("(A|B|C|D)_[0-9]")

rpoE_mpralm <- mpralm(
	rpoE_MPRASet, 
	design = rpoE_design,
	block = rpoE_block,
	model_type = "corr_groups", 
	aggregate = "mean")

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

rpoE_guide_results_raw <- topTable(rpoE_mpralm, coef = 2, number = Inf) %>% 
	data.table(keep.rownames = "spacer") %>% 
	inner_join(guides) %>% 
	arrange(adj.P.Val) %>% 
	inner_join(gene_data, by = c("locus_tag" = "b_name")) 


rpoE_gene_results_raw <- rpoE_guide_results_raw %>% 
	group_by(locus_tag, type) %>% 
	summarise(
		logFC = mean(logFC), 
		FDR = handle_NA(adj.P.Val))

rpoE_guide_results_raw %>% ggplot(aes(x = logFC, y = -log10(adj.P.Val))) + geom_point()

