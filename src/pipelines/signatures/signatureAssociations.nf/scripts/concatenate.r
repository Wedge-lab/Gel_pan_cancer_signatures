suppressMessages(library(tidyverse))


args = commandArgs(trailingOnly=TRUE)
path = args[1]
pattern = args[2]
output_file = args[3]


input_files <- list.files(".", pattern)

for (file in input_files) {

    results_df <- read_delim(file, delim=",")

    full_results_df <- tryCatch( bind_rows(full_results_df, results_df),
                                error=function (e) {results_df} )

}

write_delim(full_results_df, output_file, delim=",")
