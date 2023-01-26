# Library imports
pacman::p_load(
	data.table, edgeR, poolr, lme4, mpra,
	tidyverse, magrittr, broom, limma, ggallin, mice, Rtsne)

# rejected_promoters <- c("mdtA", "ydiF", "lacUV5", "j23119", "mraZ", "hdeA")
# rejected_promoters <- c("mdtA", "ydiF")
rejected_promoters <- c()

lane_spread_spacers <- fread("lane_spread_spacers.tsv.gz") %>%
	filter(!promoter %in% rejected_promoters)

sample_identifiers <- lane_spread_spacers %>% 
	filter(nucleotide == "RNA") %>% 
	select(promoter, timing, replicate, lane) %>% 
	unite(identifier, c("promoter", "timing", "replicate", "lane"), remove = FALSE)

##########################################################################################
# RNA
lane_spread_matrix_RNA <- lane_spread_spacers %>% 
	filter(nucleotide == "RNA") %>% 
	select(-promoter, -timing, -replicate, -lane, -nucleotide) %>%
	data.matrix() %>% 
	set_rownames(sample_identifiers$identifier)

# RNA MDS
lane_spread_matrix_RNA_MDS <- lane_spread_matrix_RNA %>% 
	dist(method = "canberra") %>% 
	cmdscale %>%
	data.table %>% cbind(sample_identifiers)

lane_spread_matrix_RNA_MDS %>% 
	ggplot(aes(x = V1, y = V2)) + geom_point(aes(colour = promoter)) + facet_grid(~timing) +
	ggtitle("RNA tSNE") + ggtitle("RNA MDS by Canberra") 

#################################################################################

# RNA tSNE
lane_spread_tsne_RNA <- lane_spread_matrix_RNA %>% Rtsne

lane_spread_tsne_RNA_DT <- lane_spread_tsne_RNA$Y %>% data.table

lane_spread_tsne_RNA_DT <- lane_spread_spacers %>% 
	filter(nucleotide == "RNA") %>% 
	select(promoter, timing, replicate, lane, nucleotide) %>% 
	cbind(lane_spread_tsne_RNA_DT)

lane_spread_tsne_RNA_DT %>% 
	as.data.table %>% 
	ggplot(aes(x = V1, y = V2)) + geom_point(aes(colour = promoter)) + facet_grid(~lane) +
	ggtitle("RNA tSNE")

##########################################################################################
# DNA
lane_spread_matrix_DNA <- lane_spread_spacers %>% 
	filter(nucleotide == "DNA") %>% 
	select(-promoter, -timing, -replicate, -lane, -nucleotide) %>%
	data.matrix() %>% 
	set_rownames(sample_identifiers$identifier)

# DNA MDS
lane_spread_matrix_DNA_MDS <- lane_spread_matrix_DNA %>% 
	dist(method = "canberra") %>% 
	cmdscale %>%
	data.table %>% cbind(sample_identifiers)

lane_spread_matrix_DNA_MDS %>% 
	ggplot(aes(x = V1, y = V2)) + geom_point(aes(colour = promoter)) + facet_grid(~timing) +
	ggtitle("DNA tSNE") + ggtitle("DNA MDS by Canberra") 

#################################################################################

# DNA tSNE
lane_spread_tsne_DNA <- lane_spread_matrix_DNA %>% Rtsne

lane_spread_tsne_DNA_DT <- lane_spread_tsne_DNA$Y %>% data.table

lane_spread_tsne_DNA_DT <- lane_spread_spacers %>% 
	filter(nucleotide == "DNA") %>% 
	select(promoter, timing, replicate, lane, nucleotide) %>% 
	cbind(lane_spread_tsne_DNA_DT)

lane_spread_tsne_DNA_DT %>% 
	as.data.table %>% 
	ggplot(aes(x = V1, y = V2)) + geom_point(aes(colour = promoter)) + facet_grid(~lane) +
	ggtitle("DNA tSNE")



