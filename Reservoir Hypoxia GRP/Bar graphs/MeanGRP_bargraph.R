#### Create bar graph for mean GRP
#### Written by RRB on 5/19/2020

# Load packages
library(ggplot2)
library(patchwork)

# Load data
data <- read.csv(paste("C:/Users/richa/OneDrive/Desktop/Reservoir Hypoxia GRP/Bar graphs/bargraph.csv"))

# Reorder Parameter so displays correctly on x-axis
data$Parameter<-factor(data$Reservoir, levels=c("ACL", "BOL", "PHL"))

data1<-data[ which(data$Species=="LMB"), ]
data2<-data[ which(data$Species=="SAE"), ]
data3<-data[ which(data$Species=="WHC"), ]

# Create jpeg of plot
jpeg(filename = paste("C:/Users/richa/OneDrive/Desktop/MeanGRP_bargraph.jpeg", sep = ""),
     width = 7.3, height = 3, units = "in", res = 600, restoreConsole = TRUE)

# Largemouth Bass Plot
LMBplot <- ggplot(data1, aes(fill=fcDO, y=MeanGRP, x=Reservoir)) + 
              geom_bar(position="dodge", colour="black", stat="identity") +
              geom_errorbar(aes(ymin=MeanGRP-SE_GRP, ymax=MeanGRP+SE_GRP),
                  width=.2,                    # Width of the error bars
                  position=position_dodge(.9)) +
              scale_y_continuous(expand = c(0,0)) +
              scale_fill_manual(labels = c(expression(paste("With ", italic(f)[c], "DO")), 
                                          (expression(paste("Without ", italic(f)[c], "DO")))), 
                                          values=c("white", "gray63")) +
              coord_cartesian(ylim = c(-0.005, 0.025)) +
              theme_bw() + theme(
                                  legend.position = c(0.68, 0.89),
                                 legend.title = element_blank(),
                                 legend.text = element_text(size = 11),
                                 legend.text.align = 0,
                                 legend.key.size = unit(1,"line"),
                                 legend.direction = "vertical",
                                 legend.spacing.x = unit(0.1, 'cm'),
                                 legend.spacing.y = unit(0, "cm"),
                                 legend.margin = margin(0.01, 0.01, 0.01, 0.01, "cm"),
              panel.border = element_rect(fill=NA, colour = "black",size=1), 
              text=element_text(size=14),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(), 
              axis.line = element_line(colour = "black")) +
              geom_hline(yintercept = 0) + 
              annotate("text", x = .65, y = 0.02325, label = "(A)", size = 4.5) +
              labs(y = "Mean GRP (g/g/d)")

             
SAEplot <- ggplot(data2, aes(fill=fcDO, y=MeanGRP, x=Reservoir)) + 
              geom_bar(position="dodge", colour="black", stat="identity", show.legend = FALSE) +
              geom_errorbar(aes(ymin=MeanGRP-SE_GRP, ymax=MeanGRP+SE_GRP),
                  width=.2,                    # Width of the error bars
                  position=position_dodge(.9)) +
              scale_y_continuous(expand = c(0,0)) +
              scale_fill_manual(values=c("white", "gray63")) +
              coord_cartesian(ylim = c(-0.005, 0.025)) +
              theme_bw() + theme(
              panel.border = element_rect(fill=NA, colour = "black",size=1),
              text=element_text(size=14),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(), 
              axis.line = element_line(colour = "black")) +
              geom_hline(yintercept = 0) + 
              annotate("text", x = .65, y = 0.02325, label = "(B)", size = 4.5)

WHCplot <- ggplot(data3, aes(fill=fcDO, y=MeanGRP, x=Reservoir)) + 
              geom_bar(position="dodge", colour="black", stat="identity", show.legend = FALSE) +
              geom_errorbar(aes(ymin=MeanGRP-SE_GRP, ymax=MeanGRP+SE_GRP),
                  width=.2,                    # Width of the error bars
                  position=position_dodge(.9)) +
              scale_y_continuous(expand = c(0,0)) +
              scale_fill_manual(values=c("white", "gray63")) +
              coord_cartesian(ylim = c(-0.005, 0.025)) +
              theme_bw() + theme(
              panel.border = element_rect(fill=NA, colour = "black",size=1),
              text=element_text(size=14),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(), 
              axis.line = element_line(colour = "black")) +
              geom_hline(yintercept = 0) +
              annotate("text", x = .65, y = 0.02325, label = "(C)", size = 4.5)

# Combine plots
patchwork = LMBplot + SAEplot + WHCplot

# Remove title from second subplot
patchwork[[1]] = patchwork[[1]] + theme(axis.title.x = element_blank())

patchwork[[2]] = patchwork[[2]] + theme(axis.text.y = element_blank(),
                                        axis.title.y = element_blank())

# Remove title from third subplot
patchwork[[3]] = patchwork[[3]] + theme(axis.title.x = element_blank(),
                                        axis.text.y = element_blank(),
                                        axis.title.y = element_blank())

patchwork

dev.off()