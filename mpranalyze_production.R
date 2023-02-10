library(pacman)
p_load(data.table, tidyverse, magrittr, MPRAnalyze, pheatmap)

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
annotations$promoter <- factor(annotations$promoter, levels = c("constitutive", unique(promoters)))

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

# Perform quantification analysis
rci_mpranalyze_quant <- analyzeQuantification(
	rci_mpranalyze,
	dnaDesign = ~ timing + promoter,
	rnaDesign = ~ timing + promoter)

# Save quantification results
rci_mpranalyze_quant %>% saveRDS("processed_results/rci_mpranalyze_quant")

# Perform comparative analysis
rci_mpranalyze_comp <- analyzeComparative(
	rci_mpranalyze,
	dnaDesign = ~ timing + promoter,
	rnaDesign = ~ timing + promoter,
	reducedDesign = ~ 1,
	fit.se = TRUE)

# Save comparative results
rci_mpranalyze_comp %>% saveRDS("processed_results/rci_mpranalyze_comp")


# Wrangle promoter results
promoter_results <- data.table()
for (i in levels(rci_mpranalyze_comp@dnaAnnot$promoter)[-1]) {
	one_result <- testCoefficient(rci_mpranalyze_comp, "promoter", i) %>% 
		data.table(keep.rownames = "spacer") %>% 
		inner_join(guides)
	setDT(one_result)
	set(one_result, j = "promoter", value = i)
	promoter_results <- rbind(promoter_results, one_result)
}

# Save promoter results
timing_results %>% fwrite("processed_results/timing_results.tsv", sep = "\t")

# Calculate timing results
timing_results <- testCoefficient(rci_mpranalyze_comp, "timing", "T2") %>%
	data.table(keep.rownames = "spacer") %>% 
	inner_join(guides) %>% arrange((fdr))

# Save timing results to file
promoter_results %>% fwrite("processed_results/promoter_results.tsv", sep = "\t")

# Create a column indicating if a locus_tag is unknown
unknown <- fread("yome/S3 No Information.tsv")
yome <- fread("yome/S1 y-ome Genes.tsv")

yome <- yome %>% mutate(unknown = case_when(locus_tag %in% unknown$locus_tag ~ TRUE, TRUE ~ FALSE))


# promoter_results is calculated by testing the coefficients of the promoter factor in the comparative analysis of the MPRA data.
promoter_results %>% inner_join(yome) %>% arrange(desc(statistic)) %>% group_by(locus_tag, promoter) %>% filter(statistic == max(statistic)) %>% View()

# getAlpha calculates the transcription rate for each spacer and associates it with the guide sequence information.
getAlpha(rci_mpranalyze_comp, by.factor = "promoter") %>% data.table(keep.rownames = "spacer") %>% inner_join(guides) %>%  inner_join(yome) %>% View

# testLrt performs a likelihood ratio test to determine if each spacer is significantly different from the constitutive promoter.
testLrt(rci_mpranalyze_comp) %>% data.table(keep.rownames = "spacer") %>% inner_join(guides) %>%  inner_join(yome) %>% View 

# testEmpirical performs an empirical test to determine if each spacer is significantly different from the non-targeting control.
testEmpirical(rci_mpranalyze_comp, twoSided = TRUE) %>% data.table(keep.rownames = "spacer") %>% inner_join(guides) %>%  inner_join(yome) %>% View 
