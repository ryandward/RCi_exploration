# Read in gene data
gene_data <- fread("gene_to_b.tsv")

# Read in data about the guides
guides <- fread("oligo_guides.tsv")

# read raw data
lane_spread_spacers <- fread("lane_spread_spacers.tsv.gz")

# Reshape data

lane_spacers <- lane_spread_spacers %>% 
	melt(
		id.vars = c("promoter", "timing", "replicate", "nucleotide", "lane"),
		variable.name = "spacer", value.name = "count") %>%
	filter(replicate %in% c("A", "B"))

lane_spacers <- lane_spacers %>% 
	group_by(nucleotide, timing, replicate, lane) %>% 
	mutate(cpm = (count/sum(count) * 10^6)) %>% ungroup %>% data.table %>% 
	mutate(spacer = as.character(spacer)) 

lane_spread_spacers_norm <- lane_spacers %>% dcast(promoter + timing + replicate + lane + nucleotide ~ spacer, value.var = "cpm") 
	
lane_spread_spacers <- lane_spacers %>% dcast(promoter + timing + replicate + lane + nucleotide ~ spacer, value.var = "count") 

