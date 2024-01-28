library(tidyverse)
library(WRS2)
library(rstatix)
library(patchwork)
library(EnvStats)

# UPLOAD -----
data <- list.files(pattern = "*.csv", full.names = F) |> 
  lapply(read_csv) |> array(dimnames = list(
    matrix.names = c("AMPH", "BIN", "CR16", "NWASP", "TBP", 
          "TKS4", "TKS4b", "TKS4L", "TKS5L", "TTP", "WIP")))

# PREPROCESING -----
## TECH REPEATS
for (i in 1:length(data)){
  data[[i]]$session <- as.factor(data[[i]]$session)
  data[[i]] |> filter(Cq > 0) |> group_by(name, session) |> 
    summarise(group = unique(group), Cq = mean(Cq), 
              session = unique(session), eff = mean(eff), 
              threshold = mean(threshold)) -> data[[i]]
}

## N0 
for (i in 1:length(data)){
  data[[i]]$N0 <- data[[i]]$threshold/(data[[i]]$eff)^data[[i]]$Cq
}

## DUPLICATES
###
for (i in 1:length(data)){
  data[[i]] |> group_by(name) |> 
    summarise(group = group, N0 = N0, n = n(), sd = sd(N0)) -> data[[i]]
}

### VIEW DUPL SD 
for (i in 1:length(data)){
  print(ggplot(data[[i]], aes(x = log10(sd)))+
    geom_density()+
    labs(title = names(data)[i], x = 'SD', y = ''))
}

### DUPL REMOVAL
for (i in 1:length(data)){
  names(data)[[i]] |> print()
  sd <- data[[i]]$sd |> c() |> na.omit()
  mad <- mad(sd)
  if (is.na(mad)){
    print(dim(data[[i]]))} 
  else{
    Lower <- median(sd) - 5*mad
    Upper <- median(sd) + 5*mad
    data[[i]] <- subset(data[[i]], data[[i]]$sd > Lower & data[[i]]$sd < Upper | is.na(data[[i]]$sd))
    print(dim(data[[i]]))
  } 
}

### DUPL AVERAGING
for (i in 1:length(data)){
  names(data)[[i]] |> print()
  data[[i]] |> group_by(name) |> 
    summarise(group = unique(group), N0 = mean(N0)) -> data[[i]]
  print(dim(data[[i]]))
}

for (i in 1:length(data)){
  names(data)[[i]] |> print()
  data[[i]] |> group_by(group) |> summarise(n=n()) |> print()
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

## NORMALIZED EXPRESSION
for (i in 1:length(data)){
  data[[i]]$dc <- data[[i]]$N0.y/data[[i]]$N0.x
}

## MEDIAN SUBTRACTION
for (i in 1:length(data)){
  data[[i]] |> filter(group == 'Adjacent') |> select (dc) |> drop_na() |> 
    as.list() |> lapply(median) -> m
  data[[i]]$ddc <- data[[i]]$dc - m[[1]]
}

## OUTLIERS
for (i in 1:length(data)){
  mad <- mad(data[[i]]$ddc)
  Lower <- median(data[[i]]$ddc) - 5*mad
  Upper <- median(data[[i]]$ddc) + 5*mad
  data[[i]] <- subset(data[[i]], data[[i]]$ddc > Lower & data[[i]]$ddc < Upper)
}

for (i in 1:length(data)){
  names(data)[[i]] |> print()
  data[[i]] |> group_by(group) |> summarise(n=n()) |> print()
}

# PRE ANALYSIS
## GOI ~ REF
for (i in 1:length(data)){
  print(ggplot(data[[i]], aes(x = N0.y, y = N0.x, col = group)) + 
          geom_point()+
          labs(title = names(data)[[i]], x = 'TBP', y = names(data)[[i]]))
}

for (i in 1:length(data)){
  names(data)[[i]] |> print()
  #data[[i]] |> (formula = N0.y ~ N0.x) |> print()
}
  
# ANALISYS 1 -----
## PLOTING
for (i in 1:length(data)){
  print(ggplot(data[[i]], aes(x = group, y = ddc))+
          geom_boxplot()+
          labs(title = names(data)[[i]], x = 'group', y = 'NE'))
}

## NORMALITY
shapiro <- function(tib){
  groups <- c('Adjacent', 'Luminal A', 'Luminal B HER-', 
              'Luminal B HER+', 'HER-Enriched', 'Triple Negative')
  for(i in groups){
    print(i)
    tib |> filter(group == i) |> ungroup() |> shapiro_test(ddc) |> print()
  }
}

for (i in 1:length(data)){
  names(data)[[i]] |> print()
  shapiro(data[[i]])
}

## TESTS NOVICE LEVEL
for (i in 1:length(data)){
  names(data)[i] |> print()
  data[[i]] |> anova_test(ddc ~ group) |> print()
}

for (i in 1:length(data)){
  names(data)[i] |> print()
  data[[i]] |> wilcox_test(ddc ~ group) |> print()
}

## TESTS APPRANTICE LEVEL
for (i in 1:length(data)){
  names(data)[i] |> print()
  data[[i]] |> t1way(formula = ddc ~ group) |> print()
  #data[[i]] |> lincon(formula = ddc ~ group) |> print()
}

## TESTS BOOTSTRAP
for (i in 1:length(data)){
  names(data)[i] |> print()
  #data[[i]] |> t1waybt(formula = ddc ~ group, nboot = 5000) |> print()
  data[[i]] |> mcppb20(formula = ddc ~ group, nboot = 5000) |> print()
  #data[[i]] |> Qanova(formula = ddc ~ group, q = c(0.25, 0.5, 0.75), nboot = 1000) |> print()
  #data[[i]] |> qcomhd(formula = ddc ~ group) |> print()
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
### BEND
dt |> filter(group == 'Adjacent') |> select(3:12) |> pball()
dt |> filter(group == 'Luminal A') |> select(3:12) |> pball()
dt |> filter(group == 'Luminal B HER-') |> select(3:12) |> pball()
dt |> filter(group == 'Luminal B HER+') |> select(3:12) |> pball()
dt |> filter(group == 'Triple Negative') |> select(3:12) |> pball()
dt |> filter(group == 'HER-Enriched') |> select(3:12) |> pball()

### WINS
dt |> filter(group == 'Adjacent') |> select(3:12) |> winall()
dt |> filter(group != 'Adjacent') |> select(3:12) |> winall()

## TRASH
dt$R <- dt$NWASP/dt$CR16
dt |> t1way(formula = R ~ group)
dt |> Qanova(formula = R ~ group)

## DIMENTIONALITY REDUCTION 
dt[-c(4, 7, 8, 11, 12)] |> drop_na() -> X

pca <- prcomp(X[, 3:length(X)], scale. = T, center = T)
tsne <- Rtsne::Rtsne(X[, 3:length(X)], perplexity = 4)

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