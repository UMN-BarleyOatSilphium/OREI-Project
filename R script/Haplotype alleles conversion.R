# Section 1: Create a dataframe with the list of the markers in each haplotype block

# extract the marker number from the .txt file from haploview
data <- read.table("haplo_chr1.txt", sep="\t", header=FALSE, stringsAsFactors=FALSE, fill = TRUE)


# Extract marker numbers from rows containing "BLOCK"
marker_numbers <- lapply(data$V1, function(line) {
  if (grepl("BLOCK", line)) {
    markers <- gsub("[^0-9 ]", "", line)  # Remove non-numeric characters
    markers <- trimws(markers)            # Remove leading/trailing whitespace
    markers <- strsplit(markers, " ")[[1]]# Split by space
    markers <- as.numeric(markers[-1])    # Remove the first element and convert to numeric
    markers <- markers[!is.na(markers)]  # Remove NA values
    return(markers)
  }
})

# Remove list with NULL values
marker_numbers <- marker_numbers[!sapply(marker_numbers, is.null)]


# Create a dataframe of the markers in each block and transpose the dataframe
marker_df <- data.frame(t(sapply(marker_numbers, function(x) {
  c(x, rep(NA, max(lengths(marker_numbers)) - length(x)))
})))

marker_df <- t(marker_df)

# Rename the columns
colnames(marker_df) <- paste0("chr1_hap", seq_len(ncol(marker_df)))

hap_ranges = marker_df
# count the number of markers per block
non_na_counts <- data.frame(no_of_snps=colSums(!is.na(hap_ranges)))
non_na_counts





# Section 2: Replace the marker numbers with their actual names
# import geno data for chromosome 1
dt <- read.delim("geno_chr1.hmp.txt", head = TRUE, check.names = FALSE)

geno <- dt %>%
  dplyr::select(-c(2:11))

# Replace "NN" with NA in the dataframe
cols <- 2:245

# Replace "NN" with NA in specified columns
geno[, cols][geno[, cols] == "NN"] <- NA


# replace the marker numbers with names
dt <- read.delim("geno_chr1.hmp.txt", head = TRUE, check.names = FALSE)
geno <- dt %>% dplyr::select(1)


# Map the indices to marker names
map_indices_to_snps <- function(index) {
  if (is.na(index) || index > nrow(geno) || index <= 0) {
    return(NA)
  }
  return(geno$"rs#"[index])
}

# Apply the mapping function to the hap_ranges data frame
hap_ranges_named <- hap_ranges %>% mutate_all(~ sapply(., map_indices_to_snps))





# Section 3: Construct haplotypes using ranges from the CSV file

# import geno data for chromosome 1
dt <- read.delim("geno_chr1.hmp.txt", head = TRUE, check.names = FALSE)
geno <- dt %>%
        dplyr::select(-c(2:11))

# Replace "NN" with NA in the dataframe
cols <- 2:245

# Replace "NN" with NA in specified columns
geno[, cols][geno[, cols] == "NN"] <- NA


# function to construct the haplotypes
construct_haplotypes <- function(geno, hap_ranges) {
  # Initialize an empty list to store haplotypes
  haplotypes_list <- list()
  
  # Iterate over each line (column except the first one)
  for (line in colnames(geno)[-1]) {
    # Extract the SNPs for the current line
    snps <- geno[[line]]
    
    # Initialize a temporary list to store haplotypes for the current line
    temp_haps <- list()
    
    # Iterate over each haplotype block
    for (block in colnames(hap_ranges)) {
      # Get the non-NA values in the block
      snp_indices <- hap_ranges[[block]][!is.na(hap_ranges[[block]])]
      
      # Construct the haplotype by concatenating SNPs within the defined range
      hap <- paste(snps[snp_indices], collapse = "")
      
      # Add the haplotype to the temporary list
      temp_haps[[block]] <- hap
    }
    
    # Add the haplotypes for the current line to the main list
    haplotypes_list[[line]] <- temp_haps
  }
  
  # Convert the list to a dataframe
  haplotypes_df <- as.data.frame(do.call(rbind, haplotypes_list))
  
  return(haplotypes_df)
}

# Call the function with the marker data and haplotype ranges
haplotypes_df <- construct_haplotypes(geno, hap_ranges) %>% rownames_to_column("line_name")


# Convert all columns to character and make column line_name as rownames
df_chr1 <- data.frame(lapply(haplotypes_df, as.character), stringsAsFactors=FALSE) %>% column_to_rownames("line_name")




# Section 4: Replace haplotype with count of <=2 with NA

# function to replace haplotypes count <=12 with NA
process_column <- function(column, threshold = 12) {
  # Count the occurrences
  haplo_count <- table(column)
  # Identify values occurring less than the threshold
  haplessthan_threshold <- names(haplo_count[haplo_count <= threshold])
  
  # Replace those values with NA
  column <- ifelse(column %in% haplessthan_threshold, NA, column)
  
  return(column)
}


# Apply the function to all the columns
df_chr1 <- df_chr1 %>%
  mutate(across(1:328, ~process_column(.x, threshold = 12)))

df <- df_chr1





# Section 5: Convert each haplotype that contains NA to NA
# Iterate over each cell
for (col in names(df)) {
  for (i in seq_along(df[[col]])) {
    # Check if the cell contains "NA"
    if (grepl("NA", df[[col]][i], fixed = TRUE)) {
      # Replace the cell with NA
      df[[col]][i] <- NA
    }
  }
}


unique_counts1 <- data.frame(count= sapply(df, function(x) length(unique(na.omit(x)))))
sum(unique_counts1)



# Section 6: Count unique haplotype alleles in each column, calculate frequencies and add

# function to count unique haplotypes and their frequencies in each column
hap_allele_freq <- function(col) {
  # Filter out sequences containing 'NA'
  filtered_col <- col[!grepl("NA", col)]
  # Calculate the frequency of each unique haplotype
  table(filtered_col)
  total_count <- length(filtered_col)
  freq <- table(filtered_col) / total_count
}

# Apply the function to each column
hap_freqs <- apply(df, 2, hap_allele_freq)


# count unique haplotypes and total
unique_counts <- data.frame(count= sapply(df, function(x) length(unique(na.omit(x)))))
sum(unique_counts)

haplotypes_equal_1 <-unique_counts %>% filter(count ==1)

# remove haplotype blocks with only one allele
df <- df %>% select(-rownames(haplotypes_equal_1))


# Section 7: Compute genotype frequencies and replace haplotypes with numbers.
# Assign different number for different alleles with the same frequency

compute_and_replace <- function(col) {
  # Compute genotype frequencies
  freq <- table(col) / sum(!is.na(col))
  
  # Identify haplotypes with frequencies less than or equal to 5%
  low_freq_hap <- names(freq[freq <= 0.05])
  
  # Replace haplotypes with frequencies less than or equal to 5% with 0
  col[col %in% low_freq_hap] <- 0
  
  # Identify haplotypes with frequencies greater than 5%
  high_freq_hap <- names(freq[freq > 0.05])
  
  # Sort haplotypes by frequency
  sorted_hap <- high_freq_hap[order(freq[high_freq_hap])]
  
  # Replace all other haplotypes with numbers in order of increasing frequency
  hap_number <- 1
  for (hap in sorted_hap) {
    col[col == hap] <- hap_number
    hap_number <- hap_number + 1
  }
  
  return(col)
}

# Apply the function to each column
df2 <- lapply(df, compute_and_replace)


# Convert the list to a dataframe and add rownames as in original data
df2_chr1 <- as.data.frame(df2)
rownames(df2_chr1) <- rownames(df)




# do this for all the haplotype block files generated for each chromosome in Haploview
# merge all files with haplotype alleles and numeric values from all the seven chromosomes
