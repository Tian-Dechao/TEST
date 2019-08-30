examp1 = function(n){
    x = rep(c('M', 'F'), c(n, n))
    y = rep(c('A', 'B'), c(n/2, n/2))
    y = rep(y, 2)
    set.seed(12345)
    pval = runif(n*2, 0, 0.1)
    n = sample.int(n, size=2*n, replace=T)
    df = data.frame(gender=x, 
                    group=y,
                    n = n,
                    pval = pval,
                    stringsAsFactors = F)
    return(df)
}

FilterDFByPval = function(X, cn='pval', cutoff=0.05){
    X.sub = X %>% 
        filter( ( !!rlang::sym(cn) ) <= cutoff )
    
    return(X.sub)
}

CutClustersByGeneNumber = function(V){
    # create ng group
    # bins and labels are fixed
    nb = c(0, 5, 10, 15, 20, Inf)
    nl = c('<=5', '[6-10]', '[11-15]', '[16-20]', '>20')
    ng = cut(V, breaks = nb, labels = nl, include.lowest = T, ordered_result = T)
    return(ng)    
}


CountGroupElements = function(X, grp_cols=c('gender', 'group'), varn='New'){
    # count number of elements in each group 
    X.grp = X %>%
        group_by_at(vars(one_of(grp_cols))) %>%
        summarize(!!varn :=n())
    return(X.grp)
}

CombineDFs = function(DFfirst, DFsecond, by_vars=c('gender', 'group')){
# combine two tables by merge
    DFnew = DFfirst %>%
        merge(DFsecond, by=by_vars)
    
    return(DFnew)
}

MainComparison = function(DFwp, DFp){
    # DFwp is the DF with pvalues 
    # DFp is prior clusters
    
    DFwp2 = FilterDFwpByPval(X=DFwp, cn='pval', cutoff=0.05)
    # create gene group
    DFwp2 = DFwp2 %>% 
        mutate(ng = CutClustersByGeneNumber(n))
    
    DFwp_all = DFwp2 %>% 
        group_by(gender, group) %>%
        summarize(n = n()) 
    
    DFwp_ng = DFwp2 %>% 
        group_by(gender, group, ng) %>%
        summarize(n = n()) 
}