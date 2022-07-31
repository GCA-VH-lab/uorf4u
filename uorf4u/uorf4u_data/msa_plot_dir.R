# Usage: Rscript msa_plot.R [--aa_msa path, --nt_msa path, --sd_msa path].
# Use -h option to show help message.

if (!requireNamespace("optparse", quietly = TRUE)) {
  install.packages("optparse")
}
library("optparse", quietly = T)

option_list = list(
  make_option(
    c("--aa_msa"),
    type = "character",
    default = NA,
    help = "Path to a dir with amino acid MSA fasta files.",
    metavar = "character"
  ),
  make_option(
    c("--nt_msa"),
    type = "character",
    default = NA,
    help = "Path to a dir with nucleotide MSA fasta files.",
    metavar = "character"
  ),
  make_option(
    c("--sd_msa"),
    type = "character",
    default = NA,
    help = "Path to a dir with SD (nt) MSA fasta files.",
    metavar = "character"
  )
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

if (!requireNamespace("ggmsa", quietly = TRUE)) {
  install.packages("ggmsa")
}
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}

suppressPackageStartupMessages(library("ggmsa"))
library("ggplot2")

current_dir = getwd()
types = c('nt', 'aa', 'sd')
for (type in types) {
  input_dir <- ifelse(type == "nt", opt$nt_msa,
                      ifelse(type == "aa", opt$aa_msa,
                             opt$sd_msa))
  if (!is.na(input_dir)) {
    input_dir <- file.path(current_dir, input_dir)
    output_dir = paste(input_dir, "figs", sep = '_')
    dir.create(output_dir)
    files <-
      list.files(path = input_dir,
                 pattern = "*fa",
                 full.names = F)
    color <- ifelse(type == "nt" | type == "sd", "Chemistry_NT",
                    "Chemistry_AA")
    for (file in files) {
      file_path <-  file.path(input_dir, file)
      fig <- ggmsa(file_path,
                   seq_name = TRUE,
                   char_width = 0.65,
                   color = color,
                   border = "#FFFFFF") +
        geom_seqlogo(color = color)
      fig_output <- file.path(output_dir, gsub(".fa$", '.pdf', file))
      ggsave(fig_output, plot = fig, units = 'cm', limitsize = FALSE )
    }
  }
}

?ggsave()
