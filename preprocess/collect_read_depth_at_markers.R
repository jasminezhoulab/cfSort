library(data.table)

args <- commandArgs(trailingOnly = T)
input_file <- args[1]
output_file <- args[2]
marker_file <- args[3]
MAX_INDEX <- 1045098

marker <- fread(marker_file, data.table = F, header = T)$region_index
dt <- matrix(NA, ncol = length(marker), nrow = 1)
rownames(dt) <- input_file
cur_x <- fread(input_file, data.table = F, header = T)
cur_x <- cur_x[cur_x$marker_index %in% marker, ]
cur_x <- cur_x[! duplicated(cur_x$marker_index),]
rownames(cur_x) <- as.character(cur_x$marker_index)
dt[input_file,] <- cur_x[as.character(marker), "num_read"]
dt[is.na(dt)] <- 0
dt <- data.frame(dt)
fwrite(dt, output_file, col.names = F, row.names = T, quote = F, sep = "\t")


