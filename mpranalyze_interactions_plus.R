library(pacman)
p_load(data.table, tidyverse, magrittr, MPRAnalyze, pheatmap)

# # Read barcode winners
# barcode_winners <- fread(
# 	'data_products/barcode.winners.tsv.gz',
# 	header = F,
# 	col.names = c(
# 		"barcode",
# 		"promoter",
# 		"spacer",
# 		"lib_count",
# 		"lib_total_count",
# 		"lib_size"))
# 
# barcode_stats <- fread("barcode_stats_lane_agg.tsv.gz")
# 
# lib_stats <- barcode_winners %>% 
# 	group_by(promoter, spacer) %>% 
# 	summarise(
# 		lib_count = sum(lib_count), 
# 		lib_total_count = sum(lib_total_count),  
# 		family_lib_size = sum(lib_size), 
# 		this_lib_size = n(),
# 		count_doubt = lib_total_count - lib_count,
# 		size_doubt = family_lib_size - this_lib_size,
# 		doubt_ratio = count_doubt/size_doubt)

# Read oligo guides data
guides <- fread("oligo_guides.tsv")

# Read promoter data and set constitutive promoters, move them to the beginning of the list
promoters <- c("j23119", "lacUV5", sort(setdiff(fread("promoters.tsv")$promoter, c("j23119", "lacUV5"))))

# Read spacer stats data
spacer_stats_spread <- fread("spacer_stats_spread_full.tsv.gz")

# Preprocess spacer stats data
spacer_stats <- spacer_stats_spread %>% 
	melt(id.vars = c("spacer", "nucleotide")) %>% 
	separate(variable, c("timing", "replicate", "lane", "promoter"), sep = "_") %>%
	group_by(spacer, nucleotide, timing, replicate, promoter) %>% 
	summarise(count = sum(value)) %>%
	mutate(annot = paste(timing, replicate, promoter, sep = ":")) %>% 
	arrange(annot, spacer, nucleotide)

# Create annotations data frame
annotations <- spacer_stats %>% ungroup() %>% 
	select(annot, timing, replicate, promoter) %>% 
	unique() %>% select(-annot) %>% data.frame() %>% 
	set_rownames(unique(spacer_stats$annot))

# Convert timing and replicate to factors
annotations$timing <- factor(annotations$timing)
annotations$replicate <- factor(annotations$replicate)

# Convert promoter names to "constitutive" for j23119 and lacUV5
annotations$promoter <- gsub("j23119|lacUV5", "constitutive", annotations$promoter)
annotations$promoter <- factor(annotations$promoter, levels = c("constitutive", sort(setdiff(promoters, c("j23119", "lacUV5")))))

# Join spacers data with oligo guides
spacers <- spacer_stats %>% ungroup %>% select(spacer) %>% unique() %>% inner_join(guides)

# Filter controls data
controls <- spacers %>% data.table %>% `[`( , type == "controls")

# Arrange RNA data to match spacers order
spacer_stats.RNA <- spacer_stats %>% 
	filter(nucleotide == "RNA") %>% ungroup %>% 
	select(spacer, annot, count) %>% data.table %>% 
	dcast(spacer ~ annot, value.var = "count")

# Convert RNA data to matrix
spacer_stats.RNA <- spacers %>% select(spacer) %>% left_join(spacer_stats.RNA)

# Turn RNA data into matrix
spacer_stats.RNA <- spacer_stats.RNA %>% 
	select(-spacer) %>% 
	data.matrix %>% 
	set_rownames(spacer_stats.RNA$spacer)

# Preprocess DNA data
spacer_stats.DNA <- spacer_stats %>% 
	filter(nucleotide == "DNA") %>% ungroup %>% 
	select(spacer, annot, count) %>% data.table %>% 
	dcast(spacer ~ annot, value.var = "count")

# Arrange DNA data to match spacers order
spacer_stats.DNA <- spacers %>% select(spacer) %>% left_join(spacer_stats.DNA)

# Convert DNA data to matrix
spacer_stats.DNA <- spacer_stats.DNA %>% 
	select(-spacer) %>% 
	data.matrix %>% 
	set_rownames(spacer_stats.DNA$spacer)

# Create MPRAnalyze object
rci_mpranalyze <- MpraObject(
	dnaCounts = spacer_stats.DNA, 
	rnaCounts = spacer_stats.RNA,
	dnaAnnot = annotations, 
	rnaAnnot = annotations,
	controls = controls)

# Estimate depth factors
rci_mpranalyze <- estimateDepthFactors(
	rci_mpranalyze, 
	lib.factor = c("timing", "promoter"),
	which.lib = "both",
	depth.estimator = "uq")

# Perform comparative analysis
rci_mpranalyze_intx_plus_comp <- analyzeComparative(
	rci_mpranalyze,
	dnaDesign = ~ timing + promoter,
	rnaDesign = ~ timing + promoter + timing:promoter,
	reducedDesign = ~ timing + promoter,
	fit.se = TRUE)

# Save comparative results
rci_mpranalyze_intx_plus_comp %>% saveRDS("processed_results/rci_mpranalyze_intx_plus_comp")

rci_logFC <- rci_mpranalyze_intx_plus_comp@modelFits$r.coef %>% data.table(keep.rownames = "spacer") %>% set_colnames(c("spacer", "disp", colnames(rci_mpranalyze_intx_plus_comp@designs@rnaFull)))

rci_se <- rci_mpranalyze_intx_plus_comp@modelFits$r.se  %>% data.table(keep.rownames = "spacer") %>% set_colnames(c("spacer", "disp", colnames(rci_mpranalyze_intx_plus_comp@designs@rnaFull)))

rci_results <- rci_logFC %>% 
	melt(id.vars = "spacer", variable.name = "predictor", value.name = "logFC") %>% 
	inner_join(rci_se %>% melt(id.vars = "spacer", variable.name = "predictor", value.name = "se")) %>% 
	mutate(stat_pred = (logFC / se) ^ 2 ) %>% 
	filter(predictor != "disp")


rci_results <- rci_results %>% 
	filter(logFC != 0) %>%
	inner_join(rci_mpranalyze_intx_plus_comp@modelFits$r.df %>% data.frame(r.df = .) %>% data.table(keep.rownames = "spacer")) %>%
	# mutate(stat_p.val = pchisq(stat_pred, df = r.df, lower.tail = FALSE)) %>% 
	mutate(stat_p.val = pchisq(stat_pred, df = 1, lower.tail = FALSE)) %>%
	group_by(predictor) %>% 
	mutate(fdr = p.adjust(stat_p.val, 'BH'))  %>%
	arrange(fdr)



rci_results %>% inner_join(rci_mpranalyze_intx_plus_comp@modelFits$r.df %>% data.frame(r.df = .) %>% data.table(keep.rownames = "spacer"))


rci_logFC %>% melt(id.vars = "spacer", variable.name = "predictor", value.name = "logFC") %>% inner_join(rci_se %>% melt(id.vars = "spacer", variable.name = "predictor", value.name = "se")) %>% mutate(statistic = (logFC / se) ^ 2 ) %>% arrange(desc(statistic)) %>% inner_join(guides) %>% filter(predictor != "disp") %>%  View


guides %>% inner_join(rci_mpranalyze_intx_plus_comp@modelFits$r.coef %>% data.table(keep.rownames = "spacer") %>% set_colnames(c("spacer", "disp", colnames(rci_mpranalyze_intx_plus_comp@designs@rnaFull)))) %>% melt(id.vars = c("gene", "type", "locus_tag", "offset", "spacer"), variable.name = "predictor", value.name = "logFC")  %>% inner_join(testLrt(rci_mpranalyze_intx_plus_comp) %>% data.table(keep.rownames = "spacer") %>% arrange(desc(statistic))) %>% View