#' Polar histogram designed for STRUCTURE plots
#' 
#' This function is adapted from Christophe Ladroue's work in library(phenotypicForest)
#' (https://github.com/chrislad/phenotypicForest/blob/master/R/polarHistogram.R)
#'
#' @param df A data.frame including columns of (family, item, score, value). 
#'            family denotes phenotype, such as tissue type. item denotes 
#'            sample lables, such as brain tissue no. 1. Score denotes cluster 
#'            memebership of each sample. value denotes cluster membership
#'            probabilities of each sample. famiy, item, and score vectors 
#'            are required to be factor.
#' @param guides A vector of values used for visual guides. Typically ranges 
#'                  between 0 to 100.
#' @param outerRadius Length of the radius to the outer edge of the circle.
#' @param circleProportion Percent of the ciricle (in radius) to be used in 
#'                            plotting.
#' @param palette A vector of colors for the clusters.
#' @param Family A factor vector of phenotype for each sample. The level of the famliy vector determins
#'               the order of the samples in the 
#' 
#' @return Returns a polar histogram. 
#'
#' @import plyr
#' @import ggplot2
#' @export
#'
#' @examples
#' polarHistogramStructure()
polarHistogramStructure <-function (df, 
                           family = NULL, 
                           columnNames = NULL, 
                           binSize = 1,
                           spaceItem = 0.2, spaceFamily = 1.2, 
                           familyLabelDistance = 1.2,
                           innerRadius = 0.3, outerRadius = 1,
                           guides = c(10, 20, 40, 80), 
                           alphaStart = -0.3, circleProportion = 0.8,
                           direction = "inwards", familyLabels = FALSE, 
                           normalised = TRUE,
                           palette)
{
    
    if (!is.null(columnNames)) {
        namesColumn <- names(columnNames)
        names(namesColumn) <- columnNames
        df <- plyr::rename(df, namesColumn)
    }
    
    applyLookup <- function(groups, keys, unassigned = "unassigned") {
        lookup <- rep(names(groups), sapply(groups, length, USE.NAMES = FALSE))
        names(lookup) <- unlist(groups, use.names = FALSE)
        p <- lookup[as.character(keys)]
        p[is.na(p)] <- unassigned
        p
    }
    
    if (!is.null(family))
        df$family <- applyLookup(family, df$item)
    df <- plyr::arrange(df, family, item, score)
    if(normalised)
        df <- ddply(df, .(family, item), transform, value = cumsum(value/(sum(value))))
    else {
        maxFamily <- max(plyr::ddply(df,.(family,item), summarise, total = sum(value))$total)
        df <- plyr::ddply(df, .(family, item), transform, value = cumsum(value))
        df$value <- df$value/maxFamily
    }
     
    # indexItem is a vector of numeric values assignd to each levels of the item variable
    # indexFamily is a vector of numeric values assignd to each levels of the family variable
    df <- ddply(df, .(family, item), transform, previous = c(0, head(value, length(value) - 1)))
    
    df2 <- ddply(df, .(family, item), summarise, indexItem = 1)
    df2$indexItem <- cumsum(df2$indexItem)
    df3 <- ddply(df, .(family), summarise, indexFamily = 1)
    df3$indexFamily <- cumsum(df3$indexFamily)
    df <- merge(df, df2, by = c("family", "item"))
    df <- merge(df, df3, by = "family")
    df <- plyr::arrange(df, family, item, score)
    
    affine <- switch(direction,
                     inwards = function(y) (outerRadius - innerRadius) * y + innerRadius,
                     outwards = function(y) (outerRadius - innerRadius) * (1 - y) + innerRadius,
                     stop(paste("Unknown direction")))
    
    
    df <- within(df, {
        xmin <- (indexItem - 1) * binSize + (indexItem - 1) *
            spaceItem + (indexFamily - 1) * (spaceFamily - spaceItem)
        xmax <- xmin + binSize
        ymin <- affine(1 - previous)
        ymax <- affine(1 - value)
    })
    
    if(normalised)
        guidesDF <- data.frame(xmin = rep(df$xmin, length(guides)),
                               y = rep(1 - guides/100, 1, each = nrow(df)))
    else
        guidesDF <- data.frame(xmin = rep(df$xmin, length(guides)),
                               y = rep(1 - guides/maxFamily, 1, each = nrow(df)))
    
    
    guidesDF <- within(guidesDF, {
        xend <- xmin + binSize
        y <- affine(y)
    })
    
    totalLength <- tail(df$xmin + binSize + spaceFamily, 1)/circleProportion - 0
    
    p <- ggplot(df) + geom_rect(aes(xmin = xmin, xmax = xmax,
                                    ymin = ymin, ymax = ymax, fill = score))
    readableAngle <- function(x) {
        angle <- x * (-360/totalLength) - alphaStart * 180/pi + 90
        angle + ifelse(sign(cos(angle * pi/180)) + sign(sin(angle * pi/180)) == -2, 180, 0)
    }
    
    readableJustification <- function(x) {
        angle <- x * (-360/totalLength) - alphaStart * 180/pi + 90
        ifelse(sign(cos(angle * pi/180)) + sign(sin(angle * pi/180)) == -2, 1, 0)
    }
    
    dfItemLabels <- plyr::ddply(df, .(family, item), summarize, xmin = xmin[1])
    dfItemLabels <- within(dfItemLabels, {
        x <- xmin + binSize/2
        angle <- readableAngle(xmin + binSize/2)
        hjust <- readableJustification(xmin + binSize/2)
    })
    
    # add item labels
#     p <- p + geom_text(aes(x = x, label = item, angle = angle,
#                            hjust = hjust), y = 1.02, size = 3, vjust = 0.5, data = dfItemLabels)
    
    p <- p + geom_segment(aes(x = xmin, xend = xend, y = y, yend = y),
                          colour = "white", data = guidesDF)
    
    # add guide labels
    # these are located inside the circle
    if(normalised)
        guideLabels <- data.frame(x = 0, y = affine(1 - guides/100),
                                  #                                  label = paste(guides, "% ", sep = ""))
                                  label = paste(guides/100, sep = " "))
    else
        guideLabels <- data.frame(x = 0, y = affine(1 - guides/maxFamily),
                                  label = paste(guides, " ", sep = ""))
    
    p <- p + geom_text(aes(x = x, y = y, label = label), data = guideLabels,
                       angle = -alphaStart * 180/pi, hjust = 1, size = 4)
    if (familyLabels) {
        familyLabelsDF <- aggregate(xmin ~ family, data = df,
                                    FUN = function(s) mean(s + binSize))
        familyLabelsDF <- within(familyLabelsDF, {
            x <- xmin
            angle <- xmin * (-360/totalLength) - alphaStart * 180/pi
        })
        p <- p + geom_text(aes(x = x, label = family, angle = angle),
                           data = familyLabelsDF, y = familyLabelDistance)
    }
    
    p <- p + theme(panel.background = element_blank(), axis.title.x = element_blank(),
                   axis.title.y = element_blank(), panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(), axis.text.x = element_blank(),
                   axis.text.y = element_blank(), axis.ticks = element_blank())
    
    p <- p + xlim(0, tail(df$xmin + binSize + spaceFamily, 1)/circleProportion)
    p <- p + ylim(0, outerRadius + 0.2)
    p <- p + coord_polar(start = alphaStart)
    #    p <- p + scale_fill_brewer(palette = "Set1", type = "qual")
    p <- p + scale_fill_manual(values = palette)
    p
}
