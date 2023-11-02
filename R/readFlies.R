require(magrittr)
# function for read multiple files in spicific path into a data table.
readFlies <- function(path,pattern="",header = F,drop = NULL, ...) {
  filenames = list.files(path, pattern,full.names = T)
  tables = lapply(filenames, function(i) {
    table = data.table::fread(i,header = header,drop = drop)
    table$filenames = i #增加一列标注来源
    table
  }) %>% data.table::rbindlist()
}