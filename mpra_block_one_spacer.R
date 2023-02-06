# Read in gene data
gene_data <- fread("gene_to_b.tsv")

# read raw data
lane_spread_spacers <- fread("lane_spread_spacers.tsv.gz")

# Reshape data

spacer_data <- lane_spread_spacers %>%
	melt(
		id.vars = c("promoter", "timing", "replicate", "nucleotide", "lane"),
		variable.name = "spacer", value.name = "count") %>%
	mutate(count = as.numeric(count)) %>%
	mutate(count = case_when(count == 0 ~ NA_real_, TRUE ~ count)) %>%
	filter(spacer == "CGGAGATCCGTCAGGCGTTT") %>%
	dcast(promoter + nucleotide ~ timing + replicate + lane, value.var = "count")

# spacer_data <- lane_spacers_norm_complete_DNA %>% 
# 	mutate(nucleotide = "DNA") %>%
# 	rbind(lane_spacers_norm_complete_RNA %>% mutate(nucleotide = "RNA")) %>% 
# 	filter(spacer == "TTTCTTCGCGGTATATTCAC") %>% 
# 		dcast(promoter + nucleotide ~ spacer + timing + replicate + lane, value.var = "imputed_count")
	
# Create block set
spacer_MPRASet <- MPRASet(
	RNA = spacer_data %>% filter(nucleotide == "RNA") %>% select(-promoter, -nucleotide),
	DNA = spacer_data %>% filter(nucleotide == "DNA") %>% select(-promoter, -nucleotide),
	eid = spacer_data %>% filter(nucleotide == "RNA") %>% pull(promoter) %>% as.character 
)

spacer_design <- data.frame(
	intcpt = 1, 
	T2 = grepl("T2", colnames(spacer_MPRASet))
) %>% 
	set_rownames(colnames(spacer_MPRASet))

spacer_block <- spacer_data %>% 
	colnames %>% `[`(spacer_data %>% colnames %>% grepl("_", .)) %>% 
	stringr::str_extract("(A|B|C|D)_[0-9]")
	# stringr::str_extract("[0-9]$")
	# stringr::str_extract("_(A|B|C|D)_") %>% gsub("_", "", .)

spacer_mpralm <- mpralm(
	spacer_MPRASet, 
	design = spacer_design,
	block = spacer_block,
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

spacer_guide_results_raw <- topTable(spacer_mpralm, coef = 2, number = Inf) %>% 
	data.table(keep.rownames = "promoter") %>% 
	arrange(adj.P.Val) %>%
	mutate(logFC = round(logFC, 3))

spacer_guide_results_raw %>% ggplot(aes(x = logFC, y = -log10(adj.P.Val))) + geom_point()

