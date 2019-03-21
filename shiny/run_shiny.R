
library(optparse)

# get Rdata file from option
opt <- parse_args(
  OptionParser(
    usage = "%prog results.Rdata",
    option_list=list()
  ), positional_arguments = c(0,1)
)

results_filename <- opt$args
if(length(results_filename) == 0){
  print("No results found! Using example dataset")
  results_filename <- "sQTL_results.Rdata"
}

print(paste0("Loading results from ",results_filename))
load(results_filename)


shiny::runApp(launch.browser=TRUE)
