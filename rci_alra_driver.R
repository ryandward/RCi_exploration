source("rci_signal_and_noise.R")
source("rci_alra_rna.R")
source("rci_alra_dna.R")

lane_spacers_complete_DNA %>% 
	mutate(nucleotide = "DNA") %>% 
	rbind(lane_spacers_complete_RNA %>% mutate(nucleotide = "RNA")) %>% 
	inner_join(guides) %>% 
	ggplot(aes(x = type, y = count_norm)) + 
	geom_boxplot() + 
	facet_grid(nucleotide~promoter) + 
	scale_x_discrete(guide = guide_axis(angle = 90))

complete_normalized <- lane_spacers_complete_DNA %>% 
	rename("count_norm" = "count_norm_DNA") %>% 
	inner_join(
		lane_spacers_complete_RNA %>% 
			rename("count_norm" = "count_norm_RNA")) %>% 
	mutate(ratio = count_norm_RNA/count_norm_DNA)


fwrite(complete_normalized, "complete_normalized.tsv.gz", sep = "\t")
complete_normalized <- fread("complete_normalized.tsv.gz")

complete_normalized %>% 
	inner_join(guides) %>% 
	ggplot(aes(x = type, y = log(ratio))) + 
	geom_boxplot() + 
	facet_grid(~promoter) + 
	scale_x_discrete(guide = guide_axis(angle = 90))


################################################################################

lane_spacers_norm_complete_tsne_RNA <- lane_spacers_norm_complete_RNA %>% Rtsne

lane_spacers_norm_complete_tsne_RNA_DT <- lane_spacers_norm_complete_tsne_RNA$Y %>% 
	data.table %>% 
	cbind(sample_identifiers)

lane_spacers_norm_complete_tsne_RNA_DT %>% 
	ggplot(aes(x = V1, y = V2)) + geom_point(aes(colour = promoter)) +
	ggtitle("RNA tSNE")


lane_spacers_norm_complete_tsne_DNA <- lane_spacers_norm_complete_DNA %>% Rtsne

lane_spacers_norm_complete_tsne_DNA_DT <- lane_spacers_norm_complete_tsne_DNA$Y %>% 
	data.table %>% 
	cbind(sample_identifiers)

lane_spacers_norm_complete_tsne_DNA_DT %>% 
	ggplot(aes(x = V1, y = V2)) + geom_point(aes(colour = promoter)) +
	ggtitle("DNA tSNE")


