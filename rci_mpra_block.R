genes <- fread("gene_to_b.tsv")

# reshape data
imputed_reshaped <- lane_spacers_norm_complete_DNA %>% mutate(nucleotide = "DNA") %>%
	rbind(lane_spacers_norm_complete_RNA %>% mutate(nucleotide = "RNA")) %>% 
	dcast(spacer + nucleotide ~ promoter + timing + replicate + lane, value.var = "imputed_count")

# make block set
imputed_reshaped_MPRASet <-  MPRASet(
	RNA = imputed_reshaped %>% filter(nucleotide == "RNA") %>% select(-spacer, -nucleotide),
	DNA = imputed_reshaped %>% filter(nucleotide == "DNA") %>% select(-spacer, -nucleotide),
	eid = imputed_reshaped %>% filter(nucleotide == "RNA") %>% pull(spacer))

imputed_reshaped_design <- data.frame(
	intcpt = 1, 
	T2 = grepl("T2", colnames(imputed_reshaped_MPRASet))) %>% 
	set_rownames(colnames(imputed_reshaped_MPRASet))

imputed_reshaped_block <- imputed_reshaped %>% 
	colnames %>% `[`(imputed_reshaped %>% colnames %>% grepl("_", .)) %>% 
	gsub("_T[0-9]_[AB]_[1-4]", "", .)

imputed_reshaped_mpralm <- mpralm(
	imputed_reshaped_MPRASet, 
	design = imputed_reshaped_design,
	block = imputed_reshaped_block,
	# design = imputed_reshaped_block,
	# block = imputed_reshaped_design,
		model_type = "corr_groups", 
	aggregate = "sum")

stouffer_try_catch <- function(x) {
	if (!any(is.na(x))) {
		result <- poolr::stouffer(x)
		p_value <- result$p
	} else {
		p_value <- NA
	}
	return(p_value)
}

results_block <- topTable(imputed_reshaped_mpralm, coef = 2, number = Inf) %>% 
	data.table(keep.rownames = "spacer") %>% 
	inner_join(guides) %>% 
	arrange(adj.P.Val) %>% 
	inner_join(genes, by = c("locus_tag" = "b_name")) %>% 
	group_by(locus_tag, type) %>% 
	summarise(
		logFC = mean(logFC), 
		FDR = stouffer_try_catch(adj.P.Val)) 

results_block %>% filter(type == "perfect" | type == "mismatch") %>% 
	ungroup() %>% select(locus_tag, logFC) %>% clipr::write_clip()

topTable(wide_data_block_model.spacer, coef = 2, number = Inf) %>% 
	data.table(keep.rownames = "spacer") %>% 
	inner_join(guides) %>% 
	arrange(adj.P.Val) %>% 
	inner_join(genes, by = c("locus_tag" = "b_name")) %>% 
	group_by(locus_tag, type) %>% filter(gene == 'topA') %>% ungroup %>% select(spacer) %>% inner_join(wide_data.spacer) %>% pivot_longer(-c(spacer, nucleotide)) %>% separate(name, c("promoter", "timing", "replicate")) %>% pivot_wider(id_cols = c(spacer, promoter, replicate), names_from = c(timing, nucleotide), values_from = value) %>% View

