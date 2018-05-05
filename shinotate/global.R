CRAN_packages <- c("shiny",  "tidyverse", "sqldf", "dbplyr", "DT", "RSQLite")
pacman::p_load(char=CRAN_packages)

# Install and load bioconductor packages
bioc_packages <- c("Biostrings")
install.deps(bioc_packages, repo="bioc")

con <- DBI::dbConnect(RSQLite::SQLite(), "/srv/shiny-server/bioinfo03/data/Paspaley/Mantle/Trinity_assembly/Annotation/Trinotate_db/P_maxima_mantle_Trinotate_local.sqlite")
mantle_annot <- dplyr::tbl(src_dbi(con), "ORF_annotation") %>% 
  dplyr::collect(n=Inf)
