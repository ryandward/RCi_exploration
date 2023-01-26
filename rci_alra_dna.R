library(pacman)

p_load_current_gh("rcannood/ALRA")

# Library imports
pacman::p_load(
	data.table, edgeR, poolr, lme4, mpra,
	tidyverse, magrittr, broom, limma, ggallin, mice)

lane_spread_matrix_DNA <- lane_spread_spacers %>% 
	filter(nucleotide == "DNA") %>% 
	select(-promoter, -timing, -replicate, -lane, -nucleotide) %>%
	data.matrix() %>%
	set_rownames(sample_identifiers$identifier)

guides <- fread("oligo_guides.tsv")

essential_guides <- guides %>% 
	filter(type == "mismatch") %>% 
	select(locus_tag) %>% unique() %>% 
	inner_join(guides) %>% 
	filter(type == "perfect")

guides[spacer %in% essential_guides$spacer, type := "perfect (essential)"]

lane_spread_matrix_DNA[is.na(lane_spread_matrix_DNA)] <- 0

gc()

lane_spread_matrix_DNA_norm <- normalize_data(lane_spread_matrix_DNA)


##########################################################################################

k_choice <- choose_k(lane_spread_matrix_DNA_norm, noise_start = 80)

# honestly I don't get this
data.table(
	x = 1:length(k_choice$d), 
	y = k_choice$d) %>% 
	ggplot(aes(x = x, y = y)) + 
	geom_line(size = 0.5) + 
	geom_vline(xintercept = k_choice$k) + 
	geom_point() + 
	xlab("Index") +
	scale_x_continuous(breaks = seq(10, 100, 10)) + 
	ylab('s_i') + 
	ggtitle('Singular values')

data.table(
	x = 2:length(k_choice$d), 
	y = diff(k_choice$d)) %>% 
	ggplot(aes(x = x, y = y)) + 
	geom_line(size = 0.5) + 
	geom_vline(xintercept = k_choice$k) + 
	geom_point() + 
	xlab("Index") +
	scale_x_continuous(breaks = seq(10, 100, 10)) + 
	ylab('s_{i} - s_{i-1}') + 
	ggtitle('Singular value spacings')

##########################################################################################

# lane_spacers_norm_complete_DNA <- alra(lane_spread_matrix_DNA_norm, k = k_choice$k)[[3]]
lane_spacers_norm_complete_DNA <- alra(lane_spread_matrix_DNA_norm, k = 5)[[3]] %>% 
	exp() %>%  `-`(1) %>% round(3) %>%
	set_rownames(sample_identifiers$identifier) %>% 
	t %>% 
	data.table(keep.rownames = "spacer") %>% 
	melt(id.vars = "spacer", value.name = "imputed_count") %>% 
	separate(variable, c("promoter", "timing", "replicate", "lane"))

