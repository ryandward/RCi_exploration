abundant_rci_mpranalyze_intx_plus_comp <- readRDS("processed_results/abundant_rci_mpranalyze_intx_plus_comp")

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

testLrt(abundant_rci_mpranalyze_intx_plus_comp) %>% arrange(fdr) %>% 
	data.table(keep.rownames = "spacer") %>% 
	rename("lrt_stat" = "statistic") %>% 
	rename("lrt_p" = "pval") %>% 
	rename("lrt_fdr" = "fdr") %>% 
	select(spacer, lrt_stat, lrt_p, lrt_fdr) %>% 
	inner_join(rci_results) %>% inner_join(guides) %>%
	arrange(desc(lrt_stat), desc(stat_pred)) %>% 
	filter(lrt_fdr < 0.05 & fdr < 0.05) %>% group_by(spacer) %>% 
	filter(stat_pred == max(stat_pred)) %>% View

#check guides 
abundant_spacer_stats %>% group_by(timing, replicate, nucleotide, promoter) %>% mutate(count = 10^6*count/sum(count))  %>% data.table %>% inner_join(guides) %>% filter(spacer == "GTTCGATTGCCACCGCAATC") %>% dcast(gene + spacer + offset + nucleotide + timing ~ promoter + replicate, value.var = "count") %>% arrange(offset)


