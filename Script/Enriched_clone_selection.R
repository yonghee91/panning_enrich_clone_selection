install.packages('ggplot2')
install.packages('kernlab')
library(ggplot2)
library(kernlab)

# Set the library name (Changeable)
targetSample <- 'Lib2'

# Set the directories, and file names (panning frequency transition file, positive clone sequence information file) (Changeable)
workDir <- 'Data/test'
writeDir <- 'Output/test'
fileName <- sprintf('%s_panning_frequency_transition.csv', targetSample)
refFileName <- 'positive_clone_seq_info.csv'

# Load data of panning frequency data and positive clones sequence data
refDf <- read.csv(file.path(workDir, refFileName))
refCDR3s <- refDf$HCDR3_AA
refCDR3s <- unique(refCDR3s)

targetDf <- read.csv(file.path(workDir, fileName))

# Run PCA for frequency transition data
estDf <- targetDf[,c(4:ncol(targetDf))]
model <- prcomp(estDf, scale=TRUE)

# PCA analysis data for predicting positive clones by enrichment pattern
writeDf <- data.frame(PC1=model$x[,1], PC2=model$x[,2])
writeDf <- cbind(writeDf, estDf)
if (abs(min(writeDf$PC1)) > abs(max(writeDf$PC1))){
  PC1_sign = F
}else{
  PC1_sign = T
}
if (abs(min(writeDf$PC2)) > abs(max(writeDf$PC2))){
  PC2_sign = F
}else{
  PC2_sign = T
}

if (PC1_sign==T){
  refX1 <- min(abs(writeDf$PC1))
  refX2 <- max(writeDf$PC1)
}else{
  refX1 <- min(abs(writeDf$PC1))
  refX2 <- min(writeDf$PC1)
}
if (PC2_sign==T){
  refY1 <- min(abs(writeDf$PC2))
  refY2 <- max(writeDf$PC2)
}else{
  refY1 <- min(abs(writeDf$PC2))
  refY2 <- min(writeDf$PC2)
}

a = (refY2-refY1) / (refX2-refX1)
b = (refX2*refY1 - refX1*refY2) / (refX2 - refX1)

# Check which PC implies enrichment pattern
tmp <- writeDf[writeDf$PC1==refX2,]
if (apply(tmp[,4:(3+ncol(estDf)-1)], 1, mean)[[1]] > tmp[,3]){
  pos_axis = 'PC1'
}else{
  pos_axis = 'PC2'
}

# Set distance threshold to look up
distThresh <- 0.39

dist <- abs(a*writeDf$PC1 - writeDf$PC2 + b)/(a^2 + 1)^(1/2)
writeDf$dist <- dist
writeDf$positive_predictive <- F

writeDf$dist_from_ref <- ((writeDf$PC1-refX1)^2 + (writeDf$PC2-(a*refX1 + b))^2)^0.5
if(pos_axis=='PC1'){
  writeDf[(writeDf$dist / writeDf$dist_from_ref > distThresh) & (abs(a*writeDf$PC1+b)-abs(writeDf$PC2) > 0) & writeDf$dist_from_ref > 0.0001,]$positive_predictive = T
}else if(pos_axis=='PC2'){
  writeDf[(writeDf$dist / writeDf$dist_from_ref > distThresh) & (abs(a*writeDf$PC1+b)-abs(writeDf$PC2) < 0) & writeDf$dist_from_ref > 0.0001,]$positive_predictive = T
}
writeDf$positive_predictive <- factor(writeDf$positive_predictive, levels = c(T, F))

writeFileName <- paste(targetSample, 'PCA_positive_predictive.png', sep='_')
p1 <- ggplot(data = writeDf, aes(x=PC1, y=PC2, color = positive_predictive)) +
  geom_point(alpha=0.7, size=3) +
  scale_color_manual('positive_predictive', values=c('deeppink', 'slategray3'))
p1 <- p1 + ggtitle('')
p1 <- p1 + labs(x = 'PC1', y = 'PC2')
p1 <- p1 + theme_bw()
p1 <- p1 + theme(axis.title = element_text(size = 15, face='bold'),
                 axis.text.x = element_text(size = 12),
                 axis.text.y = element_text(size = 12),
                 legend.position = 'right')
ggsave(filename = file.path(writeDir, writeFileName), plot = p1, device = 'png', width = 5, height = 4, limitsize = FALSE)


# Save output files
writeFileName <- paste(targetSample, 'clones.csv', sep='_')
cloneDf <- targetDf
cloneDf$positivity <- writeDf$positive_predictive
cloneDf$PC1 <- writeDf$PC1
cloneDf$PC2 <- writeDf$PC2
write.csv(file = file.path(writeDir, writeFileName), x=cloneDf, row.names = F)


### Validation with positive clones which were experimentally confirmed.
# Draw PCA plot for checking which PC (principal component) indicated enriched clones to be positive
# Find positive clones in panning frequency data
matchResult <- apply(targetDf, 1, function(x){x[[3]] %in% refCDR3s})
writeDf <- data.frame(PC1=model$x[,1], PC2=model$x[,2], positive_confirmed=matchResult)
writeDf$positive_confirmed <- factor(writeDf$positive_confirmed, levels = c('TRUE', 'FALSE'))

writeDir <- 'Output/test'
writeFileName <- paste(targetSample, 'PCA_plot_for_validation.png', sep='_')
p1 <- ggplot(data = writeDf, aes(x=PC1, y=PC2, color = positive_confirmed)) +
  geom_point(alpha=0.7, size=3) +
  scale_color_manual('positive_confirmed', values=c('red1', 'slategray3'))
p1 <- p1 + ggtitle('')
p1 <- p1 + labs(x = 'PC1', y = 'PC2')
p1 <- p1 + theme_bw()
p1 <- p1 + theme(axis.title = element_text(size = 15, face='bold'),
                 axis.text.x = element_text(size = 12),
                 axis.text.y = element_text(size = 12),
                 legend.position = 'right')
ggsave(filename = file.path(writeDir, writeFileName), plot = p1, device = 'png', width = 5, height = 4, limitsize = FALSE)
