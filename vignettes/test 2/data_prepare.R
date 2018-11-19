
MouseJaitinSpleen.topicFit <- readRDS("../inst/extdata/MouseJaitinSpleen-topicFit.rds");
save(out, file="../data/MouseJaitinSpleen.topicFit.rda")
data("MouseJaitinSpleen.topicFit")
devtools::use_data(MouseJaitinSpleen.topicFit, overwrite=TRUE)

MouseDeng2014.FitGoM <- get(load("../data/MouseDeng2014.FitGoM.rda"))
devtools::use_data(MouseDeng2014.FitGoM, overwrite=TRUE)

MouseDeng2014.FitGoM <- list("clust_3"=MouseDeng2014.FitGoM$clust_3,
                             "clust_6"=MouseDeng2014.FitGoM$clust_6)
devtools::use_data(MouseDeng2014.FitGoM, overwrite=TRUE)

MouseJaitinSpleen.FitGoM <- get(load("../data/MouseJaitinSpleen.FitGoM.rda"))
devtools::use_data(MouseJaitinSpleen.FitGoM, overwrite=TRUE)

gtex_omega <- read.table("../inst/extdata/omega_cis_genes_brain_2.txt")
gtex_theta <- read.table("../inst/extdata/theta_cis_genes_brain_2.txt")
GTExV6Brain.FitGoM <- list("omega"=gtex_omega, "theta"=gtex_theta);
devtools::use_data(GTExV6Brain.FitGoM, overwrite=TRUE)

AbundanceGoM <- readRDS("../data/topic_clus_2_maps_abundance.rds")
devtools::use_data(AbundanceGoM, overwrite=TRUE)
