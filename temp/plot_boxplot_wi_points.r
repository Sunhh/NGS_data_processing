library(ggplot2)
library(dplyr)


a1 <- read.table('geno_pheno', header=T, stringsAsFactors=F, row.names=1, sep="\t", na.strings=c("", "NA"))

svid <- 'v39109636_48'
dataset <- 'FleshBrix_19XXSq'

a1[ a1[,svid] == ".|.", svid] <- NA


k1 <- a1[,svid] == "0|0" & !is.na(a1[,svid]) & !is.na(a1[, dataset])
k2 <- a1[,svid] == "1|1" & !is.na(a1[,svid]) & !is.na(a1[, dataset])
# k3 <- a1[,svid] %in% c("0|1", "1|0") & !is.na(a1[,svid]) & !is.na(a1[, dataset])

#   Value = c(a1[k1, dataset], a1[k2, dataset], a1[k3, dataset]),
data <- data.frame(
  Value = c(a1[k1, dataset], a1[k2, dataset]),
  Group = factor(c(rep('REF', sum(k1)), rep('ALT', sum(k2))), levels= c('REF', 'ALT'))
)

# Calculate the sample sizes for each group
sample_sizes <- data %>%
  dplyr::count(Group) %>%
  dplyr::mutate(label = paste0("n = ", n))

# Perform a t-test
test_result <- t.test(Value ~ Group, data = data)
p_value <- test_result$p.value


# Create a boxplot and add jittered points
pplot <- ggplot(data, aes(x = Group, y = Value)) +
  geom_boxplot(fill= c("tomato", "skyblue")) +
  geom_jitter(width = 0.1, color = "blue", size = 2, alpha = 0.5) +
  labs(x = NULL, y = "Values") +
  geom_text(
    data = sample_sizes,
    aes(x = Group, y = max(data$Value) + 4, label = label),
    vjust = -0.5
  ) +
  geom_text(aes(x = 1.5, y = max(data$Value), label = sprintf("p = %.3f", p_value)), vjust = -0.5) +
  ggtitle(dataset)

print(pplot)


