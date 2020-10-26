# Plotting and Data viz with R

# Reading in data

```R
library(dplyr)
library(tidyverse)

tblr <- read.csv("data.csv",sep=",")

```

# tidyverse

## Filter data

## Add a new column with mutate

```R
NewTbl <- tblr %>% mutate(BinKb = Start / 1000)

head(NewTbl)
```
## Histogram

## Heatmap

## Treeview
