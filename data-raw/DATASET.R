dlpfc <-
  list(spatial = readRDS('data-raw/LIBD_spatial_3639cells_2cols.rds'), count_exp = readRDS('data-raw/LIBD_mrna_count_18011genes_3639cells.rds'))

usethis::use_data(dlpfc)
