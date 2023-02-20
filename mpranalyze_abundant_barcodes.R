library(pacman)
p_load(data.table, tidyverse, magrittr, MPRAnalyze, pheatmap)

barcode_stats <- fread("barcode_stats_lane_agg.tsv.gz")

low_counts <- c("mdtA", "ydiF")

# abundant_barcodes <- barcode_stats %>%
# 	filter(nucleotide == "RNA") %>%
# 	filter(!promoter %in% low_counts) %>%
# 	group_by(promoter, nucleotide, timing, replicate) %>%
# 	mutate(cpm = count/sum(count)*1e6) %>% 
# 	filter(cpm > 5) %>% filter(n() >= 2)  %>% arrange(barcode, replicate, timing)

# fwrite(abundant_barcodes, "abundant_barcodes.tsv.gz", sep = "\t")

abundant_barcodes <- fread("abundant_barcodes.tsv.gz")

# abundant_barcode_stats <- barcode_stats %>% filter(barcode %in% abundant_barcodes$barcode)
# 
# abundant_spacer_stats <- abundant_barcode_stats %>% group_by(promoter, spacer, timing, replicate, nucleotide) %>% summarise(count = sum(count))

# abundant_spacer_stats <- abundant_spacer_stats %>% 
# 	mutate(annot = paste(timing, replicate, promoter, sep = ":")) %>% 
# 	arrange(annot, spacer, nucleotide)
# 
# fwrite(abundant_spacer_stats, "abundant_spacer_stats.tsv.gz", sep = "\t")

abundant_spacer_stats <- fread("abundant_spacer_stats.tsv.gz")


# Read oligo guides data
guides <- fread("oligo_guides.tsv")

# Read promoter data and set constitutive promoters, move them to the beginning of the list
promoters <- c("j23119", "lacUV5", sort(setdiff(fread("promoters.tsv")$promoter, c("j23119", "lacUV5"))))
promoters <- setdiff(promoters, low_counts)

# Create annotations data frame
annotations <- abundant_spacer_stats %>% ungroup() %>% 
	select(annot, timing, replicate, promoter) %>% 
	unique() %>% select(-annot) %>% data.frame() %>% 
	set_rownames(unique(abundant_spacer_stats$annot))

# Convert timing and replicate to factors
annotations$timing <- factor(annotations$timing)
annotations$replicate <- factor(annotations$replicate)

# Convert promoter names to "constitutive" for j23119 and lacUV5
annotations$promoter <- gsub("j23119|lacUV5", "constitutive", annotations$promoter)
annotations$promoter <- factor(annotations$promoter, levels = c("constitutive", sort(setdiff(promoters, c("j23119", "lacUV5")))))

# Join spacers data with oligo guides
spacers <- abundant_spacer_stats %>% ungroup %>% select(spacer) %>% unique() %>% inner_join(guides)

# Filter controls data
controls <- spacers %>% data.table %>% `[`( , type == "controls")

# Arrange RNA data to match spacers order
abundant_spacer_stats.RNA <- abundant_spacer_stats %>% 
	filter(nucleotide == "RNA") %>% ungroup %>% 
	select(spacer, annot, count) %>% data.table %>% 
	dcast(spacer ~ annot, value.var = "count")

# Convert RNA data to matrix
abundant_spacer_stats.RNA <- spacers %>% select(spacer) %>% left_join(abundant_spacer_stats.RNA)

# Turn RNA data into matrix
abundant_spacer_stats.RNA <- abundant_spacer_stats.RNA %>% 
	select(-spacer) %>% 
	data.matrix %>% 
	set_rownames(abundant_spacer_stats.RNA$spacer)

# Preprocess DNA data
abundant_spacer_stats.DNA <- abundant_spacer_stats %>% 
	filter(nucleotide == "DNA") %>% ungroup %>% 
	select(spacer, annot, count) %>% data.table %>% 
	dcast(spacer ~ annot, value.var = "count")

# Arrange DNA data to match spacers order
abundant_spacer_stats.DNA <- spacers %>% select(spacer) %>% left_join(abundant_spacer_stats.DNA)

# Convert DNA data to matrix
abundant_spacer_stats.DNA <- abundant_spacer_stats.DNA %>% 
	select(-spacer) %>% 
	data.matrix %>% 
	set_rownames(abundant_spacer_stats.DNA$spacer)

abundant_spacer_stats.DNA[is.na(abundant_spacer_stats.DNA)] <- 0
abundant_spacer_stats.RNA[is.na(abundant_spacer_stats.RNA)] <- 0

# Create MPRAnalyze object
rci_mpranalyze <- MpraObject(
	dnaCounts = abundant_spacer_stats.DNA, 
	rnaCounts = abundant_spacer_stats.RNA,
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
abundant_rci_mpranalyze_intx_plus_comp <- analyzeComparative(
	rci_mpranalyze,
	dnaDesign = ~ timing + promoter,
	rnaDesign = ~ timing + promoter + timing:promoter,
	reducedDesign = ~ timing + promoter,
	fit.se = TRUE)

# Save comparative results
abundant_rci_mpranalyze_intx_plus_comp %>% saveRDS("processed_results/abundant_rci_mpranalyze_intx_plus_comp")

rci_logFC <- abundant_rci_mpranalyze_intx_plus_comp@modelFits$r.coef %>% data.table(keep.rownames = "spacer") %>% set_colnames(c("spacer", "disp", colnames(abundant_rci_mpranalyze_intx_plus_comp@designs@rnaFull)))

rci_se <- abundant_rci_mpranalyze_intx_plus_comp@modelFits$r.se  %>% data.table(keep.rownames = "spacer") %>% set_colnames(c("spacer", "disp", colnames(abundant_rci_mpranalyze_intx_plus_comp@designs@rnaFull)))

rci_results <- rci_logFC %>% 
	melt(id.vars = "spacer", variable.name = "predictor", value.name = "logFC") %>% 
	inner_join(rci_se %>% melt(id.vars = "spacer", variable.name = "predictor", value.name = "se")) %>% 
	mutate(stat_pred = (logFC / se) ^ 2 ) %>% 
	filter(predictor != "disp")


rci_results <- rci_results %>% 
	filter(logFC != 0) %>%
	inner_join(abundant_rci_mpranalyze_intx_plus_comp@modelFits$r.df %>% data.frame(r.df = .) %>% data.table(keep.rownames = "spacer")) %>%
	# mutate(stat_p.val = pchisq(stat_pred, df = r.df, lower.tail = FALSE)) %>% 
	mutate(stat_p.val = pchisq(stat_pred, df = 1, lower.tail = FALSE)) %>%
	group_by(predictor) %>% 
	mutate(fdr = p.adjust(stat_p.val, 'BH'))  %>%
	arrange(fdr)

