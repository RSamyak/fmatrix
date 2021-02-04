install.packages(path_to_package, repos = NULL, type = "source")

library(tidyverse)
library(fmatrix)

source("sa_chain_hetero.R")
source("sa_chain_helper.R")

path_to_files <- "./"

files_list <- list.files(path_to_files)
files_list <- files_list[grepl(".trees$", files_list)]

suffix_list <- files_list

do_everything <- function(filename, suffix) {
  assign(paste0("input_", suffix), list(filename = filename))


  assign(
    paste0("output_", suffix),
    find_sa_mean(
      get(paste0("input_", suffix)),
      preprocess = TRUE,
      track = TRUE,
      schedule = "exp",
      init = NULL,
      max.iter = 1000,
      init.temp = 1000,
      alpha = .99
    )
  )

  assign(paste0("mcctree_", suffix),
         phangorn::mcc(get(paste0("output_", suffix))$data_list$trees_list))

  saveRDS(get(paste0("output_", suffix)), file = paste0("output_", suffix, ".RDS"))
  saveRDS(get(paste0("mcctree_", suffix)), file = paste0("mcctree_", suffix, ".RDS"))
}

for(i in 1:length(files_list)){
  do_everything(file.path(path_to_files, files_list[i]), suffix_list[i])
}
