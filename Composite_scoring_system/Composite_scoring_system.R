# User can use a composite scoring system that integrates both performance and usability metrics to evaluate model performance. Files `output_directory/Machine_Learning/single_modality/[FeatureName]/single_modality_metrics.csv` will be used.
# Files `output_directory/Machine_Learning/single_modality/[FeatureName]/single_modality_metrics.csv` contains the classifier and feature selection method used, followed by various performance metrics such as accuracy, precision, recall, F1 score, AUC (Area Under the Curve), computation time and memory usage (peak memory).

library(dplyr)
library(tidyr) 
library(ggplot2) 
library(viridis)

all_metric = data.frame()

# Take feature CNA as an example
CNA = read.csv("output_directory/Machine_Learning/single_modality/CNA/single_modality_metrics.csv")
CNA <- CNA[
  CNA$Classifier == "SVM" & 
    CNA$FS_Combination == "filter_FS_0.023618328", 
]
CNA$feature = "CNA"
all_metric = rbind(all_metric,CNA) 

all_metric$accuracy_normalized <- (all_metric$accuracy - min(all_metric$accuracy)) / (max(all_metric$accuracy) - min(all_metric$accuracy))
all_metric$precision_normalized <- (all_metric$precision - min(all_metric$precision)) / (max(all_metric$precision) - min(all_metric$precision))
all_metric$recall_normalized <- (all_metric$recall - min(all_metric$recall)) / (max(all_metric$recall) - min(all_metric$recall))
all_metric$auc_normalized <- (all_metric$auc - min(all_metric$auc)) / (max(all_metric$auc) - min(all_metric$auc))
all_metric$Peak_memory_normalized <- 1 - (all_metric$PeakMemory_MB - min(all_metric$PeakMemory_MB)) / (max(all_metric$PeakMemory_MB) - min(all_metric$PeakMemory_MB))
all_metric$Time_Taken_normalized <- 1 - (all_metric$TotalTime_sec - min(all_metric$TotalTime_sec)) / (max(all_metric$TotalTime_sec) - min(all_metric$TotalTime_sec))
score <- data.frame(
  ID = all_metric$feature,
  Accuracy = all_metric$accuracy_normalized,
  Precision = all_metric$precision_normalized,
  Recall = all_metric$recall_normalized,
  AUC = all_metric$auc_normalized,
  Memory = all_metric$Peak_memory_normalized,
  Time = all_metric$Time_Taken_normalized
)

# compute the composite score
score <- score %>%
  mutate(
    score = (rowMeans(select(., Accuracy, Precision, Recall, AUC)) * 0.8) +
      (rowMeans(select(., Memory, Time)) * 0.2)
  )
sorted_score <- score[order(-score$score), ]
sorted_score <- subset(sorted_score, select = -c(score))
long_score <- sorted_score %>%
  pivot_longer(cols = -ID, names_to = "type", values_to = "value")
long_score = as.data.frame(long_score)
long_score$ID <- factor(long_score$ID, levels = rev(unique(long_score$ID)))
long_score$type <- factor(long_score$type, levels = unique(long_score$type))

ggplot(long_score, aes(x = type, y = ID, fill = value)) +
  geom_tile(color = "white", linewidth = 1/3, linetype = "solid") +
  scale_fill_viridis(option = "mako", direction = 1, na.value = "white", guide = guide_legend(title = "Value")) + 
  theme_minimal() + 
  theme(panel.grid = element_blank(), 
        axis.title = element_blank()) +
  xlab("") + 
  ylab("") + 
  coord_fixed() 
