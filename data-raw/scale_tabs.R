## code to prepare `scale_tabs` dataset goes here
usethis::use_data_raw()


library(tidyverse)

# Read folder with all csv-scale tables and store as list
files <- list.files(path = "data-raw/scale_tabs/",
                    pattern = "*.csv",
                    full.names = T)
files

names <- list.files(path = "data/",
                    pattern = "*.csv") %>%
  gsub("*.csv$", "", .) %>%
  make.names()

scale_tabs <- lapply(setNames(files, names),
                     read.csv)
scale_tabs


# Export to plain text (obsolete)
dput(scale_tabs)


# Export to data-folder

#save(scale_tabs, file = "data/scale_tabs.rda")
#usethis::use_data(scale_tabs, overwrite = TRUE)  # Tidy version
