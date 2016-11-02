

##########  test_code_5 ###########

data("MouseDeng2014.FitGoM")
omega_mat <- MouseDeng2014.FitGoM$clust_6$omega;

read.data1 = function() {
    x = tempfile()
    download.file('https://cdn.rawgit.com/kkdey/singleCellRNASeqMouseDeng2014/master/data/Deng2014MouseEsc.rda', destfile=x, quiet=TRUE)
    z = get(load((x)))
    return(z)
}

Deng2014MouseESC <- read.data1()


deng.counts <- Biobase::exprs(Deng2014MouseESC)

out <- ExtractHighCorFeatures(omega_mat,
                       deng.counts,
                       num_genes=20)
