rm(list=ls())
library(dplyr)
library(tidyr)
# SE_details are SEs that have contact with bins; probably need to split this into individual rows
# needs to split SE_details into multiple rows if it 

# how to quantify long-range enhancers
# need at least one variables summarize the distance between enhancers and genes 
distance_1d =  function(S1, E1, S2, E2){
    S1 = as.numeric(S1)
    E1 = as.numeric(E1)
    S2 = as.numeric(S2)
    E2 = as.numeric(E2)
    d = c()
    for(i in 1:length(S1)){
        d[i] = .distance_1d(s1=S1[i],
                            e1=E1[i],
                            s2=S2[i],
                            e2=E2[i])
    }
    return(d)
}

.distance_1d = function(s1, e1, s2, e2){
    if(any(is.na(c(s1, e1, s2, e2)))){
        d= NA
    } else {
        
        if(s1 > e2){
            # interval s1 e1 on the right side of s2 e2
            d = s1 - e2
        } else if (e1 < s2){
            # interval s1 e1 on the right side of s2 e2
            d = s2 - e1
        } else {
            d = 0
        }
    }
    return(d)
}

gene_used = read.table('./data/gene_used/k562_used.bed', header=F, sep='\t', stringsAsFactors = F) %>%
    rename(index = V5,
           chrn = V1,
           start_gene = V2,
           end_gene=V3) %>%
    mutate(bin = sprintf('%s:%s', chrn, start_gene)) %>%
    select(index, bin, start_gene, end_gene) 

X = read.table('./data/bin_se_1000K/bin_centered_k562_1000K.txt', 
               header=T, sep='\t', stringsAsFactors = F) %>%
    mutate(SE_interacted_num = ifelse(SE_details != '', 
                     1 + stringr::str_count(SE_details, pattern=';') ,
                     0 ) ) %>%
    separate_rows(SE_details, sep=';') %>%
    separate(SE_details, c('chrn_se', 'start_se', 'end_se'), sep=':|-', remove=F)  %>%
    select(bin, SE_interacted_num, SE_details, start_se, end_se) %>%
    left_join(gene_used, by='bin') %>%
    mutate(dist_se_gene =distance_1d(S1=start_se, E1=end_se, 
                                      S2=start_gene, E2=end_gene))
# summarize by clsuter
# one cluster one row
X.clu = X %>%
    filter(!is.na(X$index)) %>%
    group_by(index) %>%
    summarize(
        ngene = n(), # #genes in a cluster
        nse = length(unique(SE_details[SE_details != ''])), # #SEs have contact with genes in a cluster
        ngene_interacted_se = length(SE_details[SE_details != ''])
    ) %>%
    mutate(pgene_interacted_se = ngene_interacted_se / ngene * 100)
# one cluster-enhancer one row
X.clu.se = X %>%
    filter(!is.na(X$index), # only look at genes assigned to clusters
          SE_details != "" ) %>% # only look at genes interacting with SEs
    group_by(index, SE_details) %>%
    summarize(
        ngene.cluser_se = n() # #genes in a cluster-se pair
    ) %>%
    left_join(X.clu, by='index') %>%
    mutate(pgene.cluster_se = ngene.cluser_se / ngene * 100)
# look at dominate ses in each cluster
X.cluse.dom = X.cluse %>%
    group_by(index) %>%
    summarize(p.gene.se.max = max(p.gene.se))
# clusters that all genes have contact to enhancers, could be one enhancer
# or multiple enhancers
# How to prepare this one 
stop()
stop()
library(ggplot2)
p = ggplot(X.clu, aes(x=ngene, y=pgene_se, color=nse)) + 
    geom_point(alpha=0.5) + 
    theme(legend.position = 'bottom')

p2 = ggplot(X.clu, aes(x=ngene, y=nse, color=pgene_se)) + 
    geom_point(alpha=0.5) + 
    theme(legend.position = 'bottom')

pdf('./figures/test.pdf', width=5, height=5)
print(p)
print(p2)
dev.off()
