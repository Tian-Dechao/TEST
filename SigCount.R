rm(list = ls())
library(dplyr)
library(ggplot2)
source('example_data.R')
DF = examp1(n=6)

# function to keep only pval < 0.05
DF2 = FilterDFByPval(X=DF, cn='pval', cutoff=0.05)
# create gene group
DF2 = DF2 %>% 
    mutate(ng = CutClustersByGeneNumber(n))

# confirmed that DF_all agrees with DF_ng
DF_all = DF2 %>% 
    group_by(gender, group) %>%
    summarize(n = n()) 

DF_ng = DF2 %>% 
    group_by(gender, group, ng) %>%
    summarize(n = n()) 
# compute number of clusters 
DF_bg = DF %>% 
    mutate(ng = CutClustersByGeneNumber(n)) 

DF_bg_all = CountGroupElements(X=DF_bg, 
                               grp_cols=c('gender', 'group'),
                               varn = 'n')
DF_bg_ng = CountGroupElements(X=DF_bg, 
                              grp_cols=c('gender', 'group', 'ng'),
                              varn = 'n') 



DF_all_new = CombineDFs(DFfirst = DF_all, DFsecond = DF_bg_all, 
           by_vars = c('gender', 'group'))

DF_ng_new = CombineDFs(DFfirst = DF_ng, DFsecond = DF_bg_ng, 
           by_vars = c('gender', 'group', 'ng'))

# compute proportion of clusters with sig overlap with prior clusters
DF_all_new = DF_all_new %>%
    mutate(prop = n.x / n.y)

DF_ng_new = DF_ng_new %>%
    mutate(prop = n.x / n.y)