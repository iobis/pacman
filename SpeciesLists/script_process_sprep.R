library(xlsx)
library(dplyr)
library(tidyr)

df <- read.xlsx2("SPREP_countries_brackish_marine_Joape.xlsx", sheetIndex = 1)

df <- df %>%
  mutate(IsInvasive = replace(IsInvasive, IsInvasive == "Null", NA)) %>%
  unite(., col = "taxonRemarks", country, establishmentMeans, IsInvasive, na.rm = TRUE, sep = ":") %>%
  group_by(scientificName) %>%
  summarize(taxonRemarks = paste0(taxonRemarks, collapse = ", ")) %>%
  mutate(references = "list by Shyama Pagad")

# get taxon matches in batches of 50
matches <- c(
  wm_records_taxamatch(df$scientificName[1:50]),
  wm_records_taxamatch(df$scientificName[51:nrow(df)])
)

# insert dummy tibbles
for (i in 1:length(matches)) {
  if (nrow(matches[[i]]) == 0) {
    matches[[i]] = tibble(scientificName = NA)
  } else if (nrow(matches[[i]]) > 1) {
    matches[[i]] <- matches[[i]][1,]
  }
}

taxa <- bind_rows(matches)
stopifnot(nrow(taxa) == nrow(df))

result <- df %>%
  mutate(scientificName_original = scientificName) %>%
  bind_cols(taxa) %>%
  select(references, scientificName_original, AphiaID, scientificName = scientificname, scientificNameAuthorship = authority, AphiaID_accepted = valid_AphiaID, scientificName_accepted = valid_name, scientificNameAuthorship_accepted = valid_authority, taxonRemarks)

write.xlsx(result, file = "sprep_list_processed.xlsx", row.names = FALSE)
