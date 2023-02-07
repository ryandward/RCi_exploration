library(pacman)
p_load(data.table, tidyverse, magrittr, mpra, pheatmap)

# load the barcode_stats_norm data
barcode_stats <- data

# barcode_stats <- fread("barcode_stats.tsv")


guides <- fread("oligo_guides.tsv")


# Select unique promoter names from the 'barcode_stats' data frame
promoters <- unique(barcode_stats$promoter)

# Sort the list of promoter names by placing 'j23119' and 'lacUV5' at the beginning
# and sorting the rest of the names in alphabetical order
promoters <- c("j23119", "lacUV5", sort(setdiff(promoters, c("j23119", "lacUV5"))))
# 
# data_summarize <- data %>% mutate(promoter = factor(promoter, levels = promoters)) %>% group_by(promoter, spacer, nucleotide, batch, timing) %>% summarise(count = sum(count))
# 
# data_summarize <- data_summarize %>% data.table
# 
# data_summarize[, batch := paste(timing, batch, sep = "_")]
# 
# spacer_stats_spread <- data_summarize %>% data.table %>% dcast(spacer + nucleotide ~ batch + promoter, value.var = "count", fill = 0)
# 
# spacer_stats_spread %>% fwrite("spacer_stats_spread_full.tsv.gz", sep = "\t")

spacer_stats_spread <- fread("spacer_stats_spread_full.tsv.gz")


spread_data_columns <- spacer_stats_spread %>% 
  colnames %>% `[`(spacer_stats_spread %>% colnames %>% grepl("^T", .))


design_promoters <- spread_data_columns %>% stringr::str_extract("[^_]+$")
design_promoters <- gsub("j23119", "constitutive", design_promoters)
design_promoters <- gsub("lacUV5", "constitutive", design_promoters)


design_timing <- spread_data_columns %>% stringr::str_extract("^T[02]")


spread_data_design <- model.matrix( ~ design_timing + design_promoters) %>% 
  set_colnames(c("intercept", (design_timing %>% unique)[-1], (design_promoters %>% unique)[-1])) %>%
  set_rownames(spread_data_columns)

spread_data_block <- spread_data_columns %>% stringr::str_extract("(A|B|C|D)_[0-9]")

spread_data_MPRAset <- MPRASet(
  RNA = filter(spacer_stats_spread, nucleotide == "RNA") %>% select(all_of(spread_data_columns)),
  DNA = filter(spacer_stats_spread, nucleotide == "DNA") %>% select(all_of(spread_data_columns)),
  eid = filter(spacer_stats_spread, nucleotide == "RNA") %>% mutate(spacer = as.character(spacer)) %>% pull(spacer))

spread_data_mpralm <- mpralm(
  spread_data_MPRAset,
  design = spread_data_design,
  block = spread_data_block,
  model_type = "corr_groups",
  aggregate = "none",
  plot = FALSE,
  normalize = TRUE)

##########################################################################################

# Initialize data.table objects to store results
results_FDR <- guides %>% select(spacer)
results_LFC <- guides %>% select(spacer)

##########################################################################################

rm(overall_results)
predictors = c("intercept", "T2", design_promoters %>% unique %>% `[`(-1))

for (i in predictors) {
  
  # Perform a generalized linear model test using the contrast in aba_contrast[,i]
  one_result <- spread_data_mpralm %>% topTable(coef = i, n = Inf) %>% data.table(keep.rownames = "spacer")
  
  # Print a message indicating which contrast is being processed
  print(paste("Processing results for", i, "..."))
  
  if (!exists("overall_results")) {
    overall_results <- one_result %>% 
      select(spacer, logFC, adj.P.Val) %>% 
      mutate(predictor = i) }
  
  else {
    overall_results <- overall_results %>% rbind(
      one_result %>% select(spacer, logFC, adj.P.Val) %>% mutate(predictor = i)) }
  
}

overall_results_summary <- overall_results %>% 
  inner_join(guides) %>% 
  group_by(gene, locus_tag, type, predictor) %>% 
  summarise(mLFC = mean(logFC), FDR = poolr::stouffer(adj.P.Val)$p) %>%
  mutate(predictor = factor(predictor, levels = predictors))

overall_results_summary$predictor <- factor(overall_results_summary$predictor, levels = predictors)

overall_results_mLFC <- overall_results_summary %>% data.table %>% dcast(gene + type ~ predictor, value.var = "mLFC")
overall_results_FDR <- overall_results_summary %>% data.table %>% dcast(gene + type ~ predictor, value.var = "FDR")

