#### Library Load ####
rm(list = ls()) # Start with clean environment
library(rentrez)
library(tidyverse)
library(Biostrings)
library(randomForest)
library(DECIPHER)
library(RSQLite)
library(ggplot2)
library(caret)
library(pROC)

#### Set up initial variables ####

TAXON1 <- "Coleoptera" # 1st Taxonomy of interest
TAXON2 <- "Lepidoptera" # 2nd Taxonomy of interest

# Max proportion of Ns allowed in sequences
prop_Ns_cutoff <- 0.02 

# Gene of interest (name, approx length)
GENE_OI <- c("CytB","800:1200") 

# Random forest downstream hyperparameters
num_trees <- 50
kmer <- 8

#### Helper functions ####

efetch_FASTAs <- function(geneOI, taxon) {
  # Preliminary nucleotide db esearch
  goi_search <- entrez_search(db = "nuccore", 
                              term = paste(taxon, "[ORGN] AND", 
                                           geneOI[1], "[Gene] AND", 
                                           geneOI[2], "[SLEN] AND biomol_genomic[PROP]"))
  
  # Repeat esearch with retmax as count from esearch result above 
  # and use history = T for use with downstream efetch
  goi_search <- entrez_search(db = "nuccore", 
                              term = paste(taxon, "[ORGN] AND", 
                                           geneOI[1], "[Gene] AND", 
                                           geneOI[2], "[SLEN] AND biomol_genomic[PROP]"), 
                              retmax = goi_search$count, 
                              use_history = T)
  
  Sys.sleep(5) # Wait 5s so don't query NCBI too frequently
  
  # Get FASTA sequences from NCBI using efetch
  goi_fastas <- entrez_fetch(db = "nuccore", 
                             rettype = "FASTA", 
                             web_history = goi_search$web_history)
  # Report results
  print(paste(goi_search$count,
              "DNA FASTA sequences fetched matching",
              geneOI[1],
              "from Taxonomic Group",
              taxon))
  
  # Write fetched sequences to fasta file then read into DNAStringSet
  write(goi_fastas, paste(taxon, "_", geneOI[1],".fasta"), sep = "\n")
  taxon1_seq <- readDNAStringSet(paste(taxon, "_", geneOI[1],".fasta"))
  
  return(taxon1_seq)
}

quality_and_length_filt <- function(tax_df, taxon) {
  
  # Report and save initial length distribution stats
  summ_stats <- summary(str_count(tax_df$Sequence))
  print(paste(taxon, "length distribution before QC"))
  print(summ_stats)
  
  print(paste(sum(str_count(tax_df$Sequence, "-")), 
              "gaps and", sum(str_count(tax_df$Sequence, "N")),
              "Ns in all", taxon,
              "sequences were detected before QC"))
  
  # Save quartiles for downstream length filtering
  q1 <- summ_stats[2]
  q3 <- summ_stats[5]
  IQR <- q3 - q1
  
  # Save number of sequences before QC
  before_length <- nrow(tax_df)
  
  
  tax_df <- tax_df %>%
    # Remove starting and ending gaps or Ns
    mutate(Seq2 = str_remove_all(Sequence, "^N+|N+$|-")) %>%
    # Remove sequences with too many Ns
    filter(str_count(Seq2, "N") <= (prop_Ns_cutoff * str_count(Sequence))) %>%
    # Remove sequences too short or too long
    filter(str_count(Seq2) >= (q1 - 1.5 * IQR) &
             str_count(Seq2) <= (q3 + 1.5 * IQR))
  
  # Report length distribution stats after QC
  print(paste(taxon, "length distribution after QC"))
  print(summary(str_count(tax_df$Sequence)))
  
  print(paste(sum(str_count(tax_df$Seq2, "-")), 
              "gaps and", sum(str_count(tax_df$Seq2, "N")),
              "Ns in all", taxon,
              "sequences were detected after QC"))
  
  # Report number of sequences before and after QC
  print(paste(taxon, "sequences before QC:", before_length))
  print(paste(taxon, "sequences after QC:", nrow(tax_df)))
  
  return(tax_df)
}

#### Get, Examine and Clean Sequences ####

# eFetch FASTA sequences that match taxonomic groups
tax1_seq <- efetch_FASTAs(geneOI = GENE_OI, taxon = TAXON1) 
tax2_seq <- efetch_FASTAs(geneOI = GENE_OI, taxon = TAXON2) 

# Prelimimary examination of sequences
BrowseSeqs(tax1_seq)
BrowseSeqs(tax2_seq)

# Build temp dfs
tax1_df <- data.frame(Taxonomy = TAXON1, 
                      Title = names(tax1_seq), 
                      Sequence = paste(tax1_seq))
tax2_df <- data.frame(Taxonomy = TAXON2, 
                      Title = names(tax2_seq), 
                      Sequence = paste(tax2_seq))

# Parse species name from FASTA label
tax1_df$Species <- word(tax1_df$Title, start = 2L, end = 3L)
tax2_df$Species <- word(tax2_df$Title, start = 2L, end = 3L)

# Build quick table to check descriptors for sequences
table(word(tax1_df$Title, start = -1, end = -1))
table(word(tax2_df$Title, start = -1, end = -1))

# Report number of unique species
sprintf("%d initial unique species of %s", length(unique(tax1_df$Species)), TAXON1) 
sprintf("%d initial unique species of %s", length(unique(tax2_df$Species)), TAXON2) 

# QC on both sets of sequences
tax1_df_clean <- quality_and_length_filt(tax1_df, TAXON1)
tax2_df_clean <- quality_and_length_filt(tax2_df, TAXON2)

# Create dfs for QC plotting purposes
len1df_B <- tax1_df %>%
  mutate(Source = paste(TAXON1, "Before")) %>%
  mutate(Seq2 = Sequence)

len2df_B <- tax2_df %>%
  mutate(Source = paste(TAXON2, "Before")) %>%
  mutate(Seq2 = Sequence)

len1df_A <- tax1_df_clean %>%
  mutate(Source = paste(TAXON1, "After"))

len2df_A <- tax2_df_clean %>%
  mutate(Source = paste(TAXON2, "After"))

length_df <- rbind(len1df_B, len2df_B, len1df_A, len2df_A)

# Plot before and after QC sequence length distributions
ggplot(data = length_df, mapping = aes(x = as.vector(str_count(Seq2)), fill = Source)) +
  geom_density(alpha = 0.4, n = 512) +
  labs(x = "Sequence Length (bp)",
       y = "Density",
       title = sprintf("%s and %s Sequence Length Distribution Before and After QC",
                       TAXON1, TAXON2)) +
  theme_bw()

# Examine Sequences after QC
BrowseSeqs(DNAStringSet(tax1_df_clean$Seq2))
BrowseSeqs(DNAStringSet(tax2_df_clean$Seq2))

# Get same number of sequences from both
min_seq <- min(dim(tax1_df_clean)[1],dim(tax2_df_clean)[1])

set.seed(42) # Set seed for reproducible results
tax1_df_clean <- tax1_df_clean %>%
  sample_n(min_seq)

set.seed(42) # Set seed for reproducible results
tax2_df_clean <- tax2_df_clean %>%
  sample_n(min_seq)

# Report number of unique species
sprintf("%d unique species of %s after QC", length(unique(tax1_df$Species)), TAXON1) 
sprintf("%d unique species of %s after QC", length(unique(tax2_df$Species)), TAXON2) 

#### Extracting k-mer Features ####

all_taxa <- rbind(tax1_df_clean, tax2_df_clean)

features_df <- data.frame(Taxonomy = all_taxa$Taxonomy, 
                          Name = all_taxa$Title, 
                          oligonucleotideFrequency(DNAStringSet(all_taxa$Seq2), 
                                                   as.prob = T, width = kmer), 
                          verbose = TRUE)

# Check results of above
table(features_df$Taxonomy)

#### Train/Val Split ####

# Create Validation set
set.seed(42)
valid_df <- features_df %>%
  group_by(Taxonomy) %>%
  sample_n(floor(0.2 * nrow(features_df)/2))

# Check equal representation
table(valid_df$Taxonomy)

# Create Training set
set.seed(42)
train_df <- features_df %>%
  filter(!Name %in% valid_df$Name)

#### Build Initial RF Model ####

set.seed(42)
ranfor_classifier <- randomForest(x = train_df[, 3:ncol(train_df)], 
                                  y = as.factor(train_df$Taxonomy), 
                                  ntree = num_trees, 
                                  importance = T)

# Check the OOB
summary(ranfor_classifier$oob.times)

# Check performance
ranfor_classifier$confusion

# Get predictions from training and validation set
predict_val_initial <- predict(ranfor_classifier, 
                               valid_df[, 3:ncol(train_df)])
predict_train_initial <- predict(ranfor_classifier, 
                                 train_df[, 3:ncol(train_df)])

# Get detailed performance stats to report
confusionMatrix(predict_train_initial, 
                as.factor(train_df$Taxonomy))
confusionMatrix(predict_val_initial, 
                as.factor(valid_df$Taxonomy))

#### Examine Feature Importance ####

feat_importance <- importance(ranfor_classifier, type = 1)

# Create df of features and their importance to classification
feat_importance <- data.frame(Feature = rownames(feat_importance), 
                              Importance = feat_importance)

# Order by descending importance
feat_importance <- feat_importance[order(-feat_importance$MeanDecreaseAccuracy),]

# Examine 30 most important
head(feat_importance, n = 30)

# Examine 30 least important
tail(feat_importance, n = 30)

# Create character vector of features with importance above 0
refined_features <- feat_importance$Feature[feat_importance$MeanDecreaseAccuracy > 0]
length(refined_features)

#### Refine Model using Feature Selection ####

set.seed(42) # So long, and thanks for the fish!

# Create vectors for varying number of features included in the models
num_features <- c(seq(1, 100, 1))
err_rate1 <- vector("numeric", length(num_features))
err_rate2 <- vector("numeric", length(num_features))
val_rate <- vector("numeric", length(num_features))

# Loop through the different number of features, gradually adding more
for (i in 1:length(num_features)) {
  set.seed(42)
  # Build model
  ranfor_classifier_refined <- randomForest(x = train_df[, refined_features[1:(i+1)]], 
                                            y = as.factor(train_df$Taxonomy), 
                                            ntree = num_trees, 
                                            importance = F)
  # Save error rate of last iteration of RF algorithm 
  err_rate1[i] <- ranfor_classifier_refined$err.rate[num_trees,2] # Taxon1
  err_rate2[i] <- ranfor_classifier_refined$err.rate[num_trees,3] # Taxon2
  
  # Predict on validation set
  predict_val <- predict(ranfor_classifier_refined, 
                         valid_df[, refined_features[1:(i+1)]])
  # Save validation error rate
  val_rate[i] <- mean(predict_val == valid_df$Taxonomy)
}

#### Visualizations ####

# Create df for downstream graphing error rate over number of features included in model
graph_err_df <- data.frame(Features = num_features, 
                           error_taxon1 = err_rate1,
                           error_taxon2 = err_rate2, 
                           val_err_rate = (1-val_rate))

# Plot error rate during training (each taxon) and validation (combined)
ggplot(graph_err_df, aes(x = Features)) +
  geom_smooth(aes(y = error_taxon1, color = paste(TAXON1, "Training"))) +
  geom_smooth(aes(y = error_taxon2, color = paste(TAXON2, "Training"))) +
  geom_smooth(aes(y = val_err_rate, color = "Validation (Mean)")) +
  scale_color_manual(values = c("orange", "black", "blue")) +
  labs(
    title = "Error Rates vs. Features Included in RF Model",
    x = "Number of Features",
    y = "Class Error Rate", 
    color = "Dataset Type") +
  theme_bw()

# Plot 20 most important features (motifs) 
ggplot(feat_importance[order(-feat_importance$MeanDecreaseAccuracy)[1:20],], 
       aes(x = MeanDecreaseAccuracy, y = Feature)) +
  geom_point(size = 3, color = "blue") +
  labs(title = "20 Most Impactful DNA Motifs for Classification", 
       x = "Mean Decrease in Model Accuracy", y = "Motifs") +
  theme_bw()
