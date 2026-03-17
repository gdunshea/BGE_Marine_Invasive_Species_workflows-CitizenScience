# Load the database query script...
source("Scripts/Functions/FUNC_query_invasive_status_databases.R")

# Process each marker's species list through WoRMS/GBIF
# (This takes 30-60 min per marker due to API rate limits)

### Can be loaded for this dataset from saved results below:
load( "Processed_data/post_01.Database_query_WORMS_NCBI.RData")

##if new data, proceed with the process_species_list function

verified_12S <- process_species_list(
  inv_12S$invasive_status, 
  output_file = "invasive_verified_12S.csv"
)

verified_18S <- process_species_list(
  inv_18S$invasive_status, 
  output_file = "invasive_verified_18S.csv"
)

verified_COI <- process_species_list(
  inv_COI$invasive_status, 
  output_file = "invasive_verified_COI.csv"
)


save.image(file = "Processed_data/post_01.Database_query_WORMS_NCBI.RData")
