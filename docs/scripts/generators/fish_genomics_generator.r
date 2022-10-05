############################################################
# For turtles web, 07.21
# This generates the file "fish_genomics.html"
############################################################


cat("Rendering fish_genomics.rmd/html\n")
Sys.setenv(RSTUDIO_PANDOC="C:/Program Files/RStudio/bin/pandoc/")
library(rmarkdown)
library(here)
setwd(here("docs", "scripts", "generators"))
print(getwd())
output_dir = "../.."
render("../markdown/fish_genomics.rmd", output_dir = output_dir, params = list(output_dir = output_dir), quiet = TRUE)