library(vegan)
library(ggplot2)
library(ggfortify)
library(optparse)
library(dplyr)

# Set up command line options
option_list <- list(
  make_option(c("-f", "--file"), type="character", default=NULL,
              help="Path to the k-mer count data file", metavar="character"),
  make_option(c("-o", "--output"), type="character", default="pcoa_plot.pdf",
              help="Path to the output PDF file", metavar="character"))

# Parse command line options
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Validate input arguments
if (is.null(opt$file)) {
  stop("Error: Must provide a file path with -f or --file\n", call.=FALSE)
}

# Function to assign letters to unique file names within each class
assign_letters <- function(data) {
  unique_files <- unique(data$file_name)
  file_to_letter <- setNames(LETTERS[seq_along(unique_files)], unique_files)
  data$file <- file_to_letter[data$file_name]
  return(data)
}

# Read in the data
kmer_data <- read.csv(opt$file, header=TRUE, row.names = 1)
k_size <- nchar(colnames(kmer_data)[4])
kmer_data <- kmer_data %>%
  group_by(file_name) %>%
  slice(sample(1:n(), size = min(n(), 50)))

# Extract the file_name column and remove non-numeric columns from the data
file_names <- kmer_data$file_name

class_labels <- kmer_data$class
kmer_data$file_name <- NULL
kmer_data$class <- NULL
kmer_data$sample_id <- NULL
kmer_data$unique_id <- NULL

# Compute distance matrix (you can choose a different method if desired)
dist_matrix <- vegdist(kmer_data, method = "euclidean", na.rm = TRUE)

# Perform PCoA
pcoa_result <- cmdscale(dist_matrix, eig = TRUE, k = 2)

# Calculate the percentage of variance explained
var_explained <- pcoa_result$eig / sum(pcoa_result$eig) * 100

# Prepare data for plotting
pcoa_df <- as.data.frame(pcoa_result$points)
colnames(pcoa_df) <- c("PCoA1", "PCoA2")
pcoa_df$file_name <- file_names
pcoa_df$class <- class_labels

# Apply the function to each class group
pcoa_df <- pcoa_df %>%
  group_by(class) %>%
  group_modify(~ assign_letters(.)) %>%
  ungroup()

pdf(opt$output, width = 6, height = 4)
p <- ggplot(pcoa_df, aes(x = PCoA1, y = PCoA2, color = file)) +
  geom_point(size = .5, alpha = 1) + #, shape = '.') +
  gghighlight() +
  facet_wrap(~ class) +
  theme_bw() +
  labs(
    title = paste0("PCoA of ", k_size, "-mer counts, 2k subsequences drawn randomly from file (color)"),
    x = paste("PCoA1 (", format(var_explained[1], digits = 2), "%)", sep = ""),
    y = paste("PCoA2 (", format(var_explained[2], digits = 2), "%)", sep = "")
  )
print(p)
dev.off()

