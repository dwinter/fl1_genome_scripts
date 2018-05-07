 
library(HiCdatR)

gene_frags_raw <- read.table("gene_by_fragment.tsv")
res <- aggregate(V8 ~ V4, FUN=unique, data=gene_frags_raw, simplify=FALSE)
gene_bins <- res$V8
names(gene_bins) <- res$V4

#load our contact matrix and normalize results
f.source.organism.specific.code("../FL1.R")
M <- f.load.one.sample(".", "mapped_contacts_10k_r.txt", 1e4, repetitions=50)

#generate teh gene-by-gene matrix from the contact matrix
#ther _must_ by a faster way to do this...
ngene <- length(gene_bins)
G <- matrix(0, nrow=ngene, ncol=ngene)

for(i in 1:ngene){
  for(j in i:ngene){
    v <- mean(M[ gene_bins[[i]], gene_bins[[j]] ])
    G[i,j] <- v
    G[j,i] <- v
  }
  if( i %% 500 == 0){
      # this takes a while, make sure it is still running...
      message(i)
  }
}

#util functions to work on intervals 
genomic_interval <- function(chrom, start, end){
    structure(list(chrom, as.integer(start),as.integer(end)), .Names=c("chrom", "start", "end"))
}

bed_to_interval <- function(row){
    genomic_interval(row[["chrom"]], row[["start"]], row[["end"]])
}

gdist <- function(a,b){
    if(is.null(a) | is.null(b)){
        return(Inf)
    }
    if(a$chrom != b$chrom){
        return(Inf)
    }
    gap <- if (a$end > b$end) a$start - b$end else b$start - a$end
    if(gap <= 0){
        #overlap
        return(0)
    }
    gap
}

frags_dist <- function(fragA, fragB){
    min(sapply(fragA, function(x) sapply(fragB, function(y) gdist(x,y))))
}


genes_full <- read_bed("../../annotation/M3_full.bed")
gene_intervals <- apply( genes_full, 1, bed_to_interval)
gene_intervals_ordered <- gene_intervals[rownames(G)]
ngene <- length(gene_intervals_ordered)
# 1-D genome distance
D <- matrix(0, nrow=ngene, ngene)
for(i in 1:ngene){
   for(j in 1:ngene){
     D[i,j] <- gdist( gene_intervals_ordered[[i]], gene_intervals_ordered[[j]] )
    }
   if(i %% 100 == 0){
     message(i)
    }
}


cuts <- c(0,1, 5000, 20000, 100000, 1e6, 10e6, Inf)

resample_contacts <- function(nsamp, G, D){
    x <- sample(nrow(G), nsamp)
    tapply(G[x,x], cut(D[x,x], breaks=cuts, include.lowest=TRUE), mean,  na.rm=TRUE)
}

plot_v_null <- function(gene_set, G, D,nrep, return_df=TRUE){
    idx <- rownames(G) %in% gene_set
    n <- sum(idx)
    sim <- replicate(nrep, resample_contacts(n, G, D))
    null <- t(apply(sim, 1, quantile, c(0.025, 0.5, 0.975), na.rm=TRUE))
    bin <- cut(D[idx,idx], breaks=cuts, include.lowest=TRUE)
    obs <- tapply(G[idx,idx], bin, mean, na.rm=TRUE)
    df0 <- data.frame(null, obs, as.numeric(table(bin)))
    names(df0) <- c("lower", "median", "upper", "obs", "nobs")
    df0$bin <- 1:(length(cuts)-1)
    print(ggplot(df0, aes(bin, obs, ymax=upper, ymin=lower)) + geom_ribbon(fill="grey80", alpha=0.6) + geom_line() + scale_y_log10())
    if(return_df){
        return(df0)
    }
}
    




}
```
