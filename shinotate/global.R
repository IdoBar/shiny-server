
# con <- DBI::dbConnect(RSQLite::SQLite(), "/srv/shiny-server/bioinfo03/data/Paspaley/Mantle/Trinity_assembly/Annotation/Trinotate_db/P_maxima_mantle_Trinotate_local.sqlite")
# mantle_annot <- dplyr::tbl(src_dbi(con), "ORF_annotation") %>% 
#   dplyr::collect(n=Inf)
onStop(function() dbDisconnect(con))
