source("rci_alra_driver.R")

# Read in gene data
gene_data <- fread("gene_to_b.tsv")

# Reshape data
rpoE_data <- lane_spacers_norm_complete_DNA %>% 
	mutate(nucleotide = "DNA") %>%
	rbind(lane_spacers_norm_complete_RNA %>% mutate(nucleotide = "RNA")) %>% 
	filter(promoter == "rpoE") %>%
	dcast(spacer + nucleotide ~ promoter + timing + replicate + lane, value.var = "imputed_count")

# Create block set
rpoE_MPRASet <- MPRASet(
	RNA = rpoE_data %>% filter(nucleotide == "RNA") %>% select(-spacer, -nucleotide),
	DNA = rpoE_data %>% filter(nucleotide == "DNA") %>% select(-spacer, -nucleotide),
	eid = rpoE_data %>% filter(nucleotide == "RNA") %>% pull(spacer)
)

rpoE_design <- data.frame(
	intcpt = 1, 
	T2 = grepl("T2", colnames(rpoE_MPRASet))
) %>% 
	set_rownames(colnames(rpoE_MPRASet))

rpoE_block <- rpoE_data %>% 
	colnames %>% `[`(rpoE_data %>% colnames %>% grepl("_", .)) %>% 
	stringr::str_extract("(A|B)_[0-9]")

rpoE_mpralm <- mpralm(
	rpoE_MPRASet, 
	design = rpoE_design,
	block = rpoE_block,
	model_type = "corr_groups", 
	aggregate = "sum")

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

rpoE_guide_results_imputed <- topTable(rpoE_mpralm, coef = 2, number = Inf) %>% 
	data.table(keep.rownames = "spacer") %>% 
	inner_join(guides) %>% 
	arrange(adj.P.Val) %>% 
	inner_join(gene_data, by = c("locus_tag" = "b_name")) 


rpoE__results_imputed <- rpoE_results %>% 
	group_by(locus_tag, type) %>% 
	summarise(
		logFC = mean(logFC), 
		FDR = handle_NA(adj.P.Val))

rpoE_guide_results_imputed %>% ggplot(aes(x = logFC, y = -log10(adj.P.Val))) + geom_point()