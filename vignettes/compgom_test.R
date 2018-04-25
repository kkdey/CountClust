

#######################  test compGoM function   #######################################


library(CountClust)
data("MouseDeng2014.FitGoM")
length(MouseDeng2014.FitGoM)
topic_fit <- MouseDeng2014.FitGoM$clust_3


read.data <- function() {
       x <- tempfile()
       download.file(paste0("https://cdn.rawgit.com/kkdey/",
                              "singleCellRNASeqMouseDeng2014",
                              "/master/data/Deng2014MouseEsc.rda"),
                     destfile = x, quiet = TRUE)
       z <- get(load((x)))
       return(z)
       }
Deng2014MouseESC <-read.data()

deng.counts <- Biobase::exprs(Deng2014MouseESC)

out <- compGoM(data = t(deng.counts),
        model = MouseDeng2014.FitGoM)

out <- compGoM(data = t(deng.counts),
               model = MouseDeng2014.FitGoM$clust_3)

tt = 10
omega1=matrix(rbind(gtools::rdirichlet(tt*10,c(3,4,2,6)),
                    gtools::rdirichlet(tt*10,c(1,4,6,3)),
                    gtools::rdirichlet(tt*10,c(4,1,2,2))), nrow=3*10*tt);
omega2=matrix(rbind(gtools::rdirichlet(tt*10,c(1,2,4,4)),
                    gtools::rdirichlet(tt*10,c(1,4,6,3)),
                    gtools::rdirichlet(tt*10,c(3,1,5,5))), nrow=3*10*tt);
out <- compare_omega(omega1, omega2)
