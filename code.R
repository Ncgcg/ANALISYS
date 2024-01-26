library(tidyverse)
library(rstatix)
library(patchwork)
library(EnvStats)

# UPLOAD -----
data <- list.files(pattern = "*.csv", full.names = F) |> 
  lapply(read_csv) |> array(dimnames = list(
    matrix.names = c("AMPH", "BIN", "CR16", "NWASP", "TBP", 
          "TKS4", "TKS4b", "TKS4L", "TKS5L", "TTP", "WIP")))

# PRE ANALYSIS -----
## 
for (i in 1:length(data)){
  data[[i]]$session <- as.factor(data[[i]]$session)
  data[[i]]$group <- as.factor(data[[i]]$group)
  data[[i]] |> filter(Cq > 0) |> group_by(name, session) |> 
    summarise(group = unique(group), Cq = mean(Cq), 
              session = unique(session), eff = mean(eff), 
              threshold = mean(threshold)) -> data[[i]]
}

## N0 
for (i in 1:length(data)){
  data[[i]]$N0 <- data[[i]]$threshold/(data[[i]]$eff)^data[[i]]$Cq
}

## N0 OUTLIERS
for (i in 1:length(data)){
  quartiles <- quantile(data[[i]]$N0, probs=c(.25, .75), na.rm = FALSE)
  IQR <- IQR(data[[i]]$N0)
  
  Lower <- quartiles[1] - 1.5*IQR
  Upper <- quartiles[2] + 1.5*IQR 
  
  data[[i]] <- subset(data[[i]], data[[i]]$N0 > Lower & data[[i]]$N0 < Upper)
}

for (i in 1:length(data)){
  data[[i]] |> group_by(name) |> 
    summarise(group = unique(group), N0 = mean(N0)) -> data[[i]]
}

## PREVIEW
for (i in 1:length(data)){
  print(ggplot(data[[i]], aes(x = group, y = N0, col = group))+
          geom_boxplot()+
          labs(title = names(data)[i], x = 'group', y = 'N0'))
}

# NORMALIZATION -----
## TBP ADDITION
data[12] <- data[5]

for (i in 1:length(data)){
  list(data[[12]], data[[i]])|> reduce(full_join, by = c('name', 'group')) |> 
  ungroup() |> drop_na(N0.x, N0.y) -> data[[i]]
}

data <- data[-c(5, 12)]

## GOI ~ REF
for (i in 1:length(data)){
  print(ggplot(data[[i]], aes(x = N0.y, y = N0.x, col = group)) + 
          geom_point()+
          labs(title = names(data)[[i]], x = 'TBP', y = names(data)[[i]]))
}

for (i in 1:length(data)){
  names(data)[[i]] |> print()
  summary(aov(data = data[[i]],  N0.y ~ N0.x+group)) |> print()
}
  
## NORMALIZED EXPRESSION
for (i in 1:length(data)){
  data[[i]]$dc <- data[[i]]$N0.y/data[[i]]$N0.x
}

## OUTLIERS
for (i in 1:length(data)){
  quartiles <- quantile(data[[i]]$dc, probs=c(.25, .75), na.rm = FALSE)
  IQR <- IQR(data[[i]]$dc)
  
  Lower <- quartiles[1] - 1.5*IQR
  Upper <- quartiles[2] + 1.5*IQR 
  
  data[[i]] <- subset(data[[i]], data[[i]]$dc > Lower & data[[i]]$dc < Upper)
}

for (i in 1:length(data)){
  names(data)[[i]] |> print()
  data[[i]] |> group_by(group) |> identify_outliers(dc) |> print()
}

for (i in 1:length(data)){
  names(data)[[i]] |> print()
  data[[i]] |> group_by(group) |> summarise(n=n()) |> print()
}

# ANALISYS 1 -----
## PREVIEW
for (i in 1:length(data)){
  print(ggplot(data[[i]], aes(x = group, y = dc))+
          geom_boxplot()+
          labs(title = names(data)[[i]], x = 'group', y = 'NE'))
}

## NORMALITY
shapiro <- function(tib){
  groups <- c('Adjacent', 'Luminal A', 'Luminal B HER-', 
              'Luminal B HER+', 'HER-Enriched', 'Triple Negative')
  for(i in groups){
    print(i)
    tib |> filter(group == i) |> ungroup() |> shapiro_test(dc) |> print()
  }
}

for (i in 1:length(data)){
  names(data)[[i]] |> print()
  shapiro(data[[i]])
}

## HOMOGENEITY
for (i in 1:length(data)){
  names(data)[[i]] |> print()
  data[[i]] |> levene_test(dc ~ group) |> print()
}

## MEDIAN SUBTRACTION
for (i in 1:length(data)){
  data[[i]] |> filter(group == 'Adjacent') |> select (dc) |> drop_na() |> 
    as.list() |> lapply(median) -> m
  data[[i]]$ddc <- data[[i]]$dc - m[[1]]
}

## PLOTTING
for (i in 1:length(data)){
  print(ggplot(data[[i]], aes(x = group, y = ddc))+
          geom_boxplot()+
          labs(title = names(data)[i], x = 'group', y = 'FC'))
}

## TESTS NOVICE LEVEL
for (i in 1:length(data)){
  names(data)[i] |> print()
  data[[i]] |> anova_test(ddc ~ group) |> print()
}

for (i in 1:length(data)){
  names(data)[i] |> print()
  data[[i]] |> tukey_hsd(ddc ~ group) |> print()
}

for (i in 1:length(data)){
  names(data)[i] |> print()
  wdata <- data[[i]]
  wdata$group <- as.character(wdata$group)
  wdata |> wilcox_test(ddc ~ group) |> print()
}

## TESTS APPRANTICE LEVEL
for (i in 1:length(data)){
  names(data)[i] |> print()
  wdata <- data[[i]]
  wdata$group <- as.character(wdata$group)
  #wdata |> WRS2::t1way(formula = ddc ~ group) |> print()
  wdata |> WRS2::lincon(formula = dc ~ group) |> print()
}

## TESTS BOOTSTRAP
for (i in 1:length(data)){
  names(data)[i] |> print()
  wdata <- data[[i]]
  wdata$group <- as.character(wdata$group)
  wdata |> WRS2::t1waybt(formula = ddc ~ group) |> print()
  wdata |> WRS2::mcppb20(formula = ddc ~ group, nboot = 5000) |> print()
  #wdata |> WRS2::Qanova(formula = ddc ~ group) |> print()
}

# ANALISTS 2 -----
## MERGE
dt <- list()
for (i in 1:length(data)){
  dt[[i]] <- data[[i]] |> select(name, group, ddc)
  colnames(dt[[i]]) <- c('name', 'group', names(data)[i])
}

dt|> reduce(full_join, by = c('name', 'group')) |> 
  ungroup() -> dt

## CORRELATIONS
GGally::ggpairs(dt[3:12])
WRS2::pbcor(dt$TKS4L, dt$TKS4b)



## TRASH
lm(data = dt, TKS4 ~ TKS4L + TKS4b)

summary(aov(data = dt, BN ~ group))

## DIMENTIONALITY REDUCTION 
dt[-c(4, 7, 8, 11, 13)] |> drop_na() -> X

pca <- prcomp(X[, 3:length(X)], scale. = T, center = T)
tsne <- Rtsne::Rtsne(X[, 3:length(X)], perplexity = 5)

pcadata <- cbind(X[, 1:2], pca$x)
tsnedata <- data.frame(X[,1:2],
                       x = tsne$Y[,1],
                       y = tsne$Y[,2])

ggplot(pcadata, aes(col = group))+ 
  geom_point(aes(x=PC1, y = PC2))
ggplot(tsnedata, aes(col = group))+ 
  geom_point(aes(x=x,y=y))

# THE END ----
# Thanks!