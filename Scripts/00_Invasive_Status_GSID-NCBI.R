### Data summairies with gathered GSID/EASIN status 

source("Scripts/Functions/FUNC_process_metadata.R")
source("Scripts/Functions/FUNC_process_invasives.R")

bge12s_cs <- readRDS("Raw_data/12S_OTUs_NCBI.RDS")
bge18S_cs <- readRDS("Raw_data/18S_OTUs_NCBI.RDS")
bgeCOI_cs <- readRDS("Raw_data/COI_OTUs_NCBI.RDS")

# Steps 1–2 only
meta_12S <- process_metadata("12S", ps_input = bge12s_cs)
meta_18S <- process_metadata("18S", ps_input = bge18S_cs)
meta_COI <- process_metadata("COI", ps_input = bgeCOI_cs)

# Steps 3–4 (only when ready)
inv_12S <- process_invasives("12S", meta_12S$phyloseq, country_col = "Sampling.area.Country")
inv_18S <- process_invasives("18S", meta_18S$phyloseq, country_col = "Sampling.area.Country")
inv_COI <- process_invasives("COI", meta_COI$phyloseq, country_col = "Sampling.area.Country")



