# cfDNAanalyzer can extract a variety of features. Users can visualize these features to better explore the full landscape of cfDNA characteristics.
# This is a script for visualizing multiple features extracted by cfDNAanalyzer, comparing the similarity of various features, and identifying redundant features.

library(tidyr)
library(dplyr)
library(tidyverse)
library(multidplyr)
library(GenomicRanges)
library(readxl)
library(ggplot2)
library(ggbreak)
library(reshape2)
library(corrplot)

# ------------ Visualize genomic copy number differences ------------
# Users can visualize genomic copy number differences between cancer and healthy samples using the CNA feature in `output_directory/Features/CNA.csv`.
# File `output_directory/Features/CNA.csv` consists of rows representing different samples. 
# The `sample` column holds the sample's file name, followed by a `label` column that indicates the sample's classification. For example, label "1" represents the cancer samples, while label "0" represents the healthy samples.

raw_data = read.csv("output_directory/Features/CNA.csv")
raw_data <- as.data.frame(t(raw_data))
colnames(raw_data) <- raw_data[1, ]
raw_data <- raw_data[-1, ] 

# divide the cancer and healthy samples according to rowname `label`.
label_row <- unlist(raw_data["label", ])
cancer_cols <- which(label_row == 1)
df_cancer <- raw_data[rownames(raw_data) != "label", cancer_cols, drop = FALSE]
healthy_cols <- which(label_row == 0)
df_healthy <- raw_data[rownames(raw_data) != "label", healthy_cols, drop = FALSE]
df_cancer <- type.convert(df_cancer, as.is = FALSE)
df_cancer$average <- rowMeans(df_cancer,na.rm = TRUE)
df_healthy <- type.convert(df_healthy, as.is = FALSE)
df_healthy$average <- rowMeans(df_healthy,na.rm = TRUE)

# create `chr`, `start`, `end`column from rownames.
df_cancer$region <- rownames(df_cancer)          
rownames(df_cancer) <- NULL                       
df_cancer <- separate(df_cancer, region, into = c("chr", "start", "end"), sep = "_")
df_cancer$chr <- sub("^.", "", df_cancer$chr)
df_cancer = df_cancer[,c("chr","start","end","average")]
df_healthy$region <- rownames(df_healthy)          
rownames(df_healthy) <- NULL                       
df_healthy <- separate(df_healthy, region, into = c("chr", "start", "end"), sep = "_")
df_healthy$chr <- sub("^.", "", df_healthy$chr)
df_healthy = df_healthy[,c("chr","start","end","average")]

# eliminate the X chromosome and convert the data in the start and end columns to numeric values while sorting the entire data frame.
df_healthy <- df_healthy[df_healthy$chr != "X", ]
df_healthy$chr <- as.numeric(df_healthy$chr)
df_healthy$start <- as.numeric(df_healthy$start)
df_healthy <- df_healthy[order(df_healthy$chr, df_healthy$start), ]
df_cancer <- df_cancer[df_cancer$chr != "X", ]
df_cancer$chr <- as.numeric(df_cancer$chr)
df_cancer$start <- as.numeric(df_cancer$start)
df_cancer <- df_cancer[order(df_cancer$chr, df_cancer$start), ]
df_healthy$combine <- 1:nrow(df_healthy) 
df_healthy$sample <- "Healthy"
df_healthy$condition <- "Healthy"
df_healthy$color <- "blue"
df_cancer$combine <- 1:nrow(df_cancer) 
df_cancer$sample <- "Cancer"
df_cancer$condition <- "Cancer"
df_cancer$color <- "red"
df_all = merge(df_cancer,df_healthy,by = c("chr", "start", "end"))
df_all <- df_all[order(df_all$chr, df_all$start), ]

# identify the regions that exhibit significant differences between the healthy and cancer samples, altering the color of these regions accordingly.
df_all$diff <- df_all$average.x / df_all$average.y
df_all$color <- ifelse(df_all$diff >= 1.5 | df_all$diff <= 2 / 3 , "green", "blue")
healthy_color = df_all$color
cancer_color <- ifelse(healthy_color == "blue", "red", healthy_color)
color = c(cancer_color,healthy_color)
df = rbind(df_cancer, df_healthy)
df$condition <- factor(df$condition, levels = c("Cancer", "Healthy"))
df$color = color
my_theme <- theme_classic()+theme(axis.text.y = element_text(color = "black", size = 10),
                                  axis.text.x = element_blank(),
                                  axis.ticks.x = element_blank(),
                                  axis.title.x = element_blank(),
                                  strip.background = element_blank(),
                                  strip.placement = "outside")

ggplot(df, aes(x = combine, y = average, group = sample, color = color)) + 
  geom_line(linewidth = 0.5) + 
  facet_grid(cols = vars(chr), rows = vars(condition), switch = "x", space = "free_x", scales = "free") + 
  coord_cartesian(xlim = NULL, ylim = c(1, 3), expand = TRUE) + 
  labs(y = "Copy Number") + 
  my_theme + 
  geom_hline(yintercept = 2, linetype = "dashed", linewidth = 0.3, color = "gray") + 
  scale_color_identity()


# ------------ Visualize the genomic fragment size difference ------------
# Users can visualize genomic fragment size differences between cancer and healthy samples using the FP feature in `output_directory/Features/FP_fragmentation_profile.csv`.
# File `output_directory/Features/FP_fragmentation_profile.csv` consists of rows representing different samples. The `sample` column holds the sample's file name, followed by a `label` column that indicates the sample's classification. For example, label "1" represents the cancer samples, while label "0" represents the healthy samples.

raw_data = read.csv("output_directory/Features/filtered_FP_fragmentation_profile.csv")
raw_data <- as.data.frame(t(raw_data))
colnames(raw_data) <- raw_data[1, ]
raw_data <- raw_data[c(-1,-2), ]

# create `chr`, `start`, `end` column from rownames.
raw_data$region <- rownames(raw_data)          
rownames(raw_data) <- NULL                       
raw_data <- separate(raw_data, region, into = c("chr", "start", "end"), sep = "_")
raw_data$chr <- gsub("chr", "", raw_data$chr)
data <- na.omit(raw_data)
data$arm = data$chr
data <- data[order(data$chr, data$start), ]

# assign the same combine value to every 50 regions for each chromosome.
data <- data %>% 
  group_by(chr) %>%
  mutate(combine = ceiling((1:length(chr)) / 50))
data <- as.data.frame(data)
data$value = as.numeric(data$value)
data <- melt(data, id.vars=c("chr", "start", "end", "arm", "combine"), 
             measure.vars=colnames(data)[-c(1, 2, 3, (ncol(data)-1), ncol(data))], variable.name = "sample")
rst <- data %>% 
  group_by(sample, chr, arm, combine) %>%
  summarize(start = dplyr::first(start),
            end = dplyr::last(end),
            binsize = n(),
            meanRatio = mean(value, na.rm = TRUE),
            varRatio = var(value, na.rm = TRUE))
rst <- rst %>% filter(binsize == 50)

# extract the mean ratio for each sample and the centered ratio for the 5 Mb regions.
meanValue = aggregate(meanRatio ~ sample, rst, FUN = mean)
meanValue <- as.data.frame(apply(meanValue, 2, FUN = rep, each = nrow(rst) / 11), stringsAsFactors = FALSE)  ## Number 11 is total number of healthy and cancer samples
meanValue$sample = as.factor(meanValue$sample)
meanValue$meanRatio = as.numeric(meanValue$meanRatio)
rst$centeredRatio <- rst$meanRatio - meanValue$meanRatio
rst$condition <- c(rep("Healthy", 1042),rep("Cancer", 4689))
rst$color <- c(rep("blue", 1042),rep("red", 4689)) ## The numbers 1042 and 4689 are the number of 5Mb regions of healthy and cancer samples in rst, respectively
rst$condition <- factor(rst$condition, levels = c("Healthy", "Cancer"))
my_theme <- theme_classic() + 
  theme(axis.text.y = element_text(color = "black", size = 10),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        strip.background = element_blank(),
        strip.placement = "outside")

ggplot(rst, aes(x=combine, y=centeredRatio, group=sample, color=color)) + 
  geom_line(linewidth=0.2) + 
  facet_grid(cols = vars(arm), rows = vars(condition), switch = "x", space = "free_x", scales = "free") +
  coord_cartesian(ylim=c(-0.4, 0.5), expand = TRUE) + 
  labs(y="Fragmentation profile") + 
  scale_color_identity() + 
  my_theme


# ------------ Find genomic motif differences ------------
# Users can find the top motifs with the most significant differences between cancer and healthy samples using the EM feature in `output_directory/Features/EM_motifs_frequency.csv` and `output_directory/Features/EM_motifs_mds.csv`.
# File `output_directory/Features/EM_motifs_frequency.csv` and `output_directory/Features/EM_motifs_mds.csv` consist of rows representing different samples. The `sample` column holds the sample's file name, followed by a `label` column that indicates the sample's classification. For example, label "1" represents the cancer samples, while label "0" represents the healthy samples.

raw_data = read.csv("output_directory/Features/EM_motifs_frequency.csv")
raw_data <- as.data.frame(t(raw_data))
colnames(raw_data) <- raw_data[1, ]
raw_data <- raw_data[-1, ] 

# divide the cancer and healthy samples according to rowname `label`.
label_row <- unlist(raw_data["label", ])
cancer_cols <- which(label_row == 1)
df_cancer <- raw_data[rownames(raw_data) != "label", cancer_cols, drop = FALSE]
healthy_cols <- which(label_row == 0)
df_healthy <- raw_data[rownames(raw_data) != "label", healthy_cols, drop = FALSE]
df_cancer <- type.convert(df_cancer, as.is = FALSE)
df_cancer$average_cancer <- rowMeans(df_cancer,na.rm = TRUE)
df_healthy <- type.convert(df_healthy, as.is = FALSE)
df_healthy$average_healthy <- rowMeans(df_healthy,na.rm = TRUE)

# we create `motif` column from rownames and extract average value.
df_cancer$motif <- rownames(df_cancer)          
rownames(df_cancer) <- NULL                       
df_cancer = df_cancer[,c("motif","average_cancer")]
df_healthy$motif <- rownames(df_healthy)          
rownames(df_healthy) <- NULL                       
df_healthy = df_healthy[,c("motif","average_healthy")]
raw_data_mds = read.csv("output_directory/Features/EM_motifs_mds.csv")
cancer_mds = raw_data_mds[raw_data_mds$label == 1, ]
healthy_mds = raw_data_mds[raw_data_mds$label == 0, ]
mean_cancer = mean(cancer_mds$MDS)
sd_cancer = sd(cancer_mds$MDS)
mean_healthy = mean(healthy_mds$MDS)
sd_healthy = sd(healthy_mds$MDS)
FREQ = merge(df_cancer, df_healthy, by = c("motif"))
FREQ$diff <- abs(FREQ$average_healthy - FREQ$average_cancer)
FREQ$rank <- rank(-FREQ$diff)

# identify the indices of the top 10 motifs in terms of frequency and change their colors accordingly.
top_rank_indices <- order(FREQ$rank)[1:10]
FREQ$color <- ifelse(1:nrow(FREQ) %in% top_rank_indices, "blue", "black")

ggplot(FREQ, aes(x = average_cancer, y = average_healthy)) +
  geom_point(aes(color = color)) +
  labs(x = "Average EM for cancer samples (average MDS: 0.954 ± 0.0044)", 
  ## Number 0.954 and 0.00044 are the mean value and standard deviation of cancer samples
       y = "Average EM for healthy samples (average MDS: 0.950 ± 0.0057)") +   
  ## Number 0.950 and 0.00057 are the mean value and standard deviation of healthy samples
  theme_minimal() +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  scale_color_identity() +
  theme(panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),  
        axis.ticks = element_line(colour = "black"),  
        axis.ticks.length = unit(0.25, "cm"))  


# ------------ Compare feature similarities and find redundant information ------------
# Users can compare feature similarities and choose features with no redundant information using features in `output_directory/Features/NOF_meanfuziness.csv` , `output_directory/Features/NOF_occupancy.csv` , `output_directory/Features/WPS_long.csv` , `output_directory/Features/OCF.csv` , `output_directory/Features/EMR_region_mds.csv` , `output_directory/Features/FPR_fragmentation_profile_regions.csv`.
# Files `output_directory/Features/NOF_meanfuziness.csv` , `output_directory/Features/NOF_occupancy.csv` , `output_directory/Features/WPS_long.csv` , `output_directory/Features/OCF.csv` , `output_directory/Features/EMR_region_mds.csv` , `output_directory/Features/FPR_fragmentation_profile_regions.csv` consist of rows representing different samples. The `sample` column holds the sample's file name, followed by a `label` column that indicates the sample's classification. For example, label "1" represents the cancer samples, while label "0" represents the healthy samples.

NO_data = read.csv("output_directory/Features/NOF_occupancy.csv")
NO_data <- as.data.frame(t(NO_data))
colnames(NO_data) <- NO_data[1, ] 
NO_data <- NO_data[c(-1,-2), ] 
NO_data$region <- rownames(NO_data)          
rownames(NO_data) <- NULL                       
NO_data <- separate(NO_data, region, into = c("chr", "start", "end"), sep = "_")

# put the basename (filename deleting ".bam") of every sample into a vector.
basename = c("basename1","basename2",...,"basenameN")
for (i in basename) {
    # region is the input BED file
    region <- read.table("/input.bed", header = F, sep = "\t" )
    names(region) <- c("chr", "start","end")

    NO = NO_data[, c("chr", "start", "end", i)]
    colnames(NO) <- c("chr", "start", "end", "NO")
    region <- merge(region, NO, all.x = TRUE)

    NF = NF_data[, c("chr", "start", "end", i)]
    colnames(NF) <- c("chr", "start", "end", "NF")
    region <- merge(region, NF, all.x = TRUE)
    
    WPS_long = WPS_long_data[, c("chr", "start", "end", i)]
    colnames(WPS_long) <- c("chr", "start", "end", "WPS_long")
    region <- merge(region, WPS_long, all.x = TRUE)
    
    OCF = OCF_data[, c("chr", "start", "end", i)]
    colnames(OCF) <- c("chr", "start", "end", "OCF")
    region <- merge(region, OCF, all.x = TRUE)
    
    EMR = EMR_data[, c("chr", "start", "end", i)]
    colnames(EMR) <- c("chr", "start", "end", "EMR")
    region <- merge(region, EMR, all.x = TRUE)
    
    FPR = FPR_data[, c("chr", "start", "end", i)]
    colnames(FPR) <- c("chr", "start", "end", "FPR")
    region <- merge(region, FPR, all.x = TRUE)
    
    region <- na.omit(region)
    region = region[,c(-1,-2,-3)]
    region <- type.convert(region, as.is = FALSE)
    assign(paste0("correlation_matrix_", i), cor(region))
}

# get all the correlation matrices' names and check is there any matrix with NA.
matrix_names <- ls(pattern = "^correlation_matrix_")
valid_matrices <- list()
for (name in matrix_names) {
  matrix <- get(name)  
  if (!any(is.na(matrix))) {  
    valid_matrices[[name]] <- matrix  
  }
}
n_rows <- nrow(valid_matrices[[1]])
n_cols <- ncol(valid_matrices[[1]])
average_matrix <- matrix(0, n_rows, n_cols)
for (matrix in valid_matrices) {
  average_matrix <- average_matrix + matrix 
}
average_matrix <- average_matrix / length(valid_matrices)

corrplot(average_matrix, 
         method = "square",                 
         type = "lower",                     
         order = "hclust",                   
         addCoef.col = "black",              
         tl.col = "black",                  
         tl.srt = 0,                        
         tl.cex = 0.8,                       
         number.cex = 1)                   


# ------------ Find the differences of PFE for genes differentially expressed between different conditions ------------
# Users can find the differences of PFE for genes differentially expressed between different conditions using features in `output_directory/Features/PFE.csv`.
# File `output_directory/Features/PFE.csv` consists of rows representing different samples. The `sample` column holds the sample's file name, followed by a `label` column that indicates the sample's classification. For example, label "1" represents the cancer samples, while label "0" represents the healthy samples.

control_genes_df = read.table("/cfDNAanalyzer/Epic-seq/code/priordata/control_hg19.txt",header = F)
control_genes = control_genes_df[,3]
under_genes <- read.table("/path/to/under_genes.txt", header = F, sep = "\t")
under_genes <- under_genes[[1]]
over_genes <- read.table("/path/to/over_genes.txt", header = F, sep = "\t")
over_genes <- over_genes[[1]]
raw_data = read.csv("output_directory/Features/PFE.csv")
raw_data <- as.data.frame(t(raw_data))
colnames(raw_data) <- raw_data[1, ]
raw_data <- raw_data[-1, ] 

# divide the cancer and healthy samples according to rowname `label`.
label_row <- unlist(raw_data["label", ])
cancer_cols <- which(label_row == 1)
df_cancer <- raw_data[rownames(raw_data) != "label", cancer_cols, drop = FALSE]
healthy_cols <- which(label_row == 0)
df_healthy <- raw_data[rownames(raw_data) != "label", healthy_cols, drop = FALSE]

# extract the mean PFE value in cancer and healthy group
df_cancer <- type.convert(df_cancer, as.is = FALSE)
df_cancer$average <- rowMeans(df_cancer,na.rm = TRUE)
df_healthy <- type.convert(df_healthy, as.is = FALSE)
df_healthy$average <- rowMeans(df_healthy,na.rm = TRUE)
df_cancer$TSS_ID <- rownames(df_cancer)          
rownames(df_cancer) <- NULL
df_cancer = df_cancer[,c("TSS_ID","average")]                       
df_healthy$TSS_ID <- rownames(df_healthy)          
rownames(df_healthy) <- NULL   
df_healthy = df_healthy[,c("TSS_ID","average")]  

# extract the normalized PFE after removing control genes
pattern <- paste(control_genes, collapse = "|")
regex <- paste0("^((", pattern, ")|(", pattern, ")_[0-9]+)$")
rows_to_remove <- grep(regex, df_cancer[[1]])
df_cancer <- df_cancer[-rows_to_remove, ]
df_cancer$normalized_PFE = scale(df_cancer$average)
normalized_df_cancer = df_cancer[,c("TSS_ID","normalized_PFE")]
df_healthy <- df_healthy[-rows_to_remove, ]
df_healthy$normalized_PFE = scale(df_healthy$average)
normalized_df_healthy = df_healthy[,c("TSS_ID","normalized_PFE")]

# extract PFE of under-expressed and over-expressed genes in cancer and healthy samples.
pattern_1 <- paste(under_genes, collapse = "|")
regex <- paste0("^((", pattern_1, ")|(", pattern_1, ")_[0-9]+)$")
rows_to_keep <- grep(regex, normalized_df_cancer[[1]])
cancer_under <- normalized_df_cancer[rows_to_keep, ]
average_cancer_under = cancer_under$normalized_PFE
pattern_1 <- paste(over_genes, collapse = "|")
regex <- paste0("^((", pattern_1, ")|(", pattern_1, ")_[0-9]+)$")
rows_to_keep <- grep(regex, normalized_df_cancer[[1]])
cancer_over <- normalized_df_cancer[rows_to_keep, ]
average_cancer_over = cancer_over$normalized_PFE
pattern_1 <- paste(under_genes, collapse = "|")
regex <- paste0("^((", pattern_1, ")|(", pattern_1, ")_[0-9]+)$")
rows_to_keep <- grep(regex, normalized_df_healthy[[1]])
healthy_under <- normalized_df_healthy[rows_to_keep, ]
average_healthy_under = healthy_under$normalized_PFE
pattern_1 <- paste(over_genes, collapse = "|")
regex <- paste0("^((", pattern_1, ")|(", pattern_1, ")_[0-9]+)$")
rows_to_keep <- grep(regex, normalized_df_healthy[[1]])
healthy_over <- normalized_df_healthy[rows_to_keep, ]
average_healthy_over = healthy_over$normalized_PFE

# create a data frame in long format that combines the PFE of the under-expressed genes
df_under <- data.frame(
  value = c(average_cancer_under, average_healthy_under),
  group = factor(c(rep("Cancer", length(average_cancer_under)), rep("Healthy", length(average_healthy_under))))
)
df_under$group <- factor(df_under$group, levels = c("Cancer", "Healthy"))
df_over <- data.frame(
  value = c(average_cancer_over, average_healthy_over),
  group = factor(c(rep("Cancer", length(average_cancer_over)), rep("Healthy", length(average_healthy_over))))
)
df_over$group <- factor(df_over$group, levels = c("Cancer", "Healthy"))
wilcox_test_result_under <- wilcox.test(value ~ group, data = df_under)
p_value_under <- wilcox_test_result_under$p.value
wilcox_test_result_over <- wilcox.test(value ~ group, data = df_over)
p_value_over <- wilcox_test_result_over$p.value

ggplot(df_under, aes(x = group, y = value)) +
  geom_boxplot() +
  labs(title = "Boxplot of Cancer PFE vs Healthy PFE in under-expressed genes", x = "Group", y = "PFE Values") +
  theme_minimal()+
  geom_text(
    x = 1.5,   
    y = max(df_under$value) * 0.9,  
    label = paste("p-value:", round(p_value_under, 4)),  
    size = 5
  ) 

ggplot(df_over, aes(x = group, y = value)) +
  geom_boxplot() +
  labs(title = "Boxplot of Cancer PFE vs Healthy PFE in over-expressed genes", x = "Group", y = "PFE Values") +
  theme_minimal() +
  geom_text(
    x = 1.5,    
    y = max(df_over$value) * 0.9,  
    label =  paste("p-value:", round(p_value_over, 4)), 
    size = 5
  )


# ------------ Visualize differences in nucleosome organization ------------
# Users can visualize differences in nucleosome organization around TF in cancer and healthy samples using feature NP in the directory `output_directory/Features/NP_site_list/`. (Take file `output_directory/Features/NP_site_list/site_list1.txt` as an exmaple)
# Files in the directory `output_directory/Features/NP_site_list/` consist of rows representing different samples. The `sample` column holds the sample's file name, followed by a `label` column that indicates the sample's classification. For example, label "1" represents the cancer samples, while label "0" represents the healthy samples.

raw_data = read.csv("/Features/NP_site_list/site_list1.txt")
df_cancer <- raw_data[raw_data$label == 1, ]
df_healthy <- raw_data[raw_data$label == 0, ]
df_cancer <- df_cancer[,c(3,134)]
mean_cancer = rowMeans(df_cancer)
df_healthy <- df_healthy[,c(3,134)]
mean_healthy = rowMeans(df_healthy)
sample <- read.csv("/Features/NP_site_list/site_list1.txt", header = F)
site_number = sample[1,c(3:134)]
site_number <- as.vector(t(site_number))

# put the coverage data of healthy and cancer samples into a data frame and transform it into long data frame form.
df = data.frame(
  Column1  = site_number,
  Cancer = mean_cancer,
  Healthy = mean_healthy
)
df_long <- melt(df, id.vars = "Column1", variable.name = "Group", value.name = "Value")

ggplot(data = df_long, aes(x = Column1, y = Value, color = Group)) +
  geom_line() +  
  labs(x = "Distance from site",
        y = "Coverage") +  
  scale_color_manual(values = c("Cancer" = "red", "Healthy" = "blue")) + 
  theme_minimal() +
  theme(panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "black"),  
        axis.ticks.length = unit(0.25, "cm"))
