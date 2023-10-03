library(vegan)
library(ggplot2)
library(ggfortify)
library(optparse)

# Set up command line options
option_list <- list(
  make_option(c("-f", "--file"), type="character", default=NULL,
              help="Path to the k-mer count data file", metavar="character"),
  make_option(c("-o", "--output"), type="character", default="pcoa_plot.pdf",
              help="Path to the output PDF file", metavar="character"),
  make_option(c("-c", "--color_by"), type="character", default="class",
              help="Column name to color points by in the PCoA plot", metavar="character")
)

# Parse command line options
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Validate input arguments
if (is.null(opt$file)) {
  stop("Error: Must provide a file path with -f or --file\n", call.=FALSE)
}

# Read in the data
kmer_data <- read.csv(opt$file, header=TRUE, row.names = 1)

# Extract the file_name column and remove non-numeric columns from the data
file_names <- kmer_data$file_name
class_labels <- kmer_data$class
kmer_data$file_name <- NULL
kmer_data$class <- NULL
kmer_data$sample_id <- NULL
kmer_data$unique_id <- NULL

# Compute distance matrix (you can choose a different method if desired)
dist_matrix <- vegdist(kmer_data, method = "euclidean")

# Perform PCoA
pcoa_result <- cmdscale(dist_matrix, eig = TRUE, k = 2)

# Prepare data for plotting
pcoa_df <- as.data.frame(pcoa_result$points)
colnames(pcoa_df) <- c("PCoA1", "PCoA2")
pcoa_df$file_name <- file_names
pcoa_df$class <- class_labels

# Select the column to color by based on the command line argument
color_column <- opt$color_by
if (!(color_column %in% names(pcoa_df))) {
  stop(paste("Error: Invalid color_by argument. Available options are:", paste(names(pcoa_df), collapse=", ")), call.=FALSE)
}

# Plot using ggplot2
pdf(opt$output)
p <- ggplot(pcoa_df, aes_string(x = "PCoA1", y = "PCoA2", color = color_column)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "PCoA of k-mer Counts", x = "PCoA1", y = "PCoA2")
print(p)
dev.off()

