# Usage: Rscript msa_plot.R --msa_fasta file.fa --output output_dir --seq_type nt|aa [--height x --width y].
# Use -h option to show help message.

if (!requireNamespace("optparse", quietly = TRUE)) {
  install.packages("optparse")
}
library("optparse", quietly = T)

option_list = list(
  make_option(
    c("--msa_fasta"),
    type = "character",
    default = NA,
    help = "Path to a MSA fasta file. (Relatively to the script!).",
    metavar = "character"
  ),
  make_option(
    c("--output"),
    type = "character",
    default = NA,
    help = "Path to the output dir.",
    metavar = "character"
  ),
  make_option(
    c("--height"),
    type = "double",
    default = NA,
    help = "Height of output figure (mm).",
    metavar = "character"
  ),
  make_option(
    c("--width"),
    type = "double",
    default = NA,
    help = "Width of output figure (mm).",
    metavar = "character"
  ),
  make_option(
    c("--seq_type"),
    type = "character",
    default = "nt",
    help = "Sequence type. Can be nt (nucleotides) or aa (aminoacids)",
    metavar = "character"
  )
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

if (! suppressMessages(requireNamespace("ggplot2", quietly = TRUE))) {
  install.packages("ggplot2")
}
if (!suppressMessages(requireNamespace("devtools", quietly=TRUE)))
{
  install.packages("devtools")
}
if (! suppressMessages(requireNamespace("ggmsa", quietly = TRUE))) {
  devtools::install_github("YuLab-SMU/ggmsa")
}

suppressPackageStartupMessages(library("ggmsa"))
suppressPackageStartupMessages(library("ggplot2"))


if (opt$seq_type != "nt" & opt$seq_type != "aa") {
  print("--seq_type argument should be 'nt' or 'aa'")
  exit()
}

file_path <- opt$msa_fasta
output = paste(opt$output, gsub(".fa$", '.pdf', basename(file_path)))
color <- ifelse(opt$seq_type == "nt", "Chemistry_NT",
                "Chemistry_AA")
fig <- ggmsa(
  file_path,
  seq_name = TRUE,
  char_width = 0.65,
  color = color,
  border = "#FFFFFF"
) + geom_seqlogo(color = color)
suppressMessages(
  ggsave(
    output,
    plot = fig,
    units = 'mm',
    limitsize = FALSE,
    width = opt$width,
    height = opt$height
  )
)
