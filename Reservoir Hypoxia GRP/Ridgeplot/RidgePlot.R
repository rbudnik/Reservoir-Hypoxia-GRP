#### Create ridgeplots with relative density of GRP values
#### Written by RRB on 5/19/2020

# Load packages
library(ggplot2)
library(ggridges)
library(patchwork)

# Load data
data1 = read.csv("C:/Users/richa/OneDrive/Desktop/Reservoir Hypoxia GRP/Ridgeplot/Ridgeplot_LMB.csv") ##read LMB data
data2 = read.csv("C:/Users/richa/OneDrive/Desktop/Reservoir Hypoxia GRP/Ridgeplot/Ridgeplot_SAE.csv") ##read SAE data
data3 = read.csv("C:/Users/richa/OneDrive/Desktop/Reservoir Hypoxia GRP/Ridgeplot/Ridgeplot_WHC.csv") ##read WHC data

# Reorder Parameter so displays correctly on y-axis
data1$Reservoir <- as.character(data1$Reservoir)
data1$Reservoir <- factor(data1$Reservoir, levels=unique(data1$Reservoir))
data1$Reservoir <- factor(data1$Reservoir, levels=c("PHL", "BOL", "ACL"))

data2$Reservoir <- as.character(data2$Reservoir)
data2$Reservoir <- factor(data2$Reservoir, levels=unique(data2$Reservoir))
data2$Reservoir <- factor(data2$Reservoir, levels=c("PHL", "BOL", "ACL"))

data3$Reservoir <- as.character(data3$Reservoir)
data3$Reservoir <- factor(data3$Reservoir, levels=unique(data3$Reservoir))
data3$Reservoir <- factor(data3$Reservoir, levels=c("PHL", "BOL", "ACL"))

# Create jpeg of plot
jpeg(filename = paste("C:/Users/richa/OneDrive/Desktop/Ridgeplot.jpeg", sep = ""),
    width = 3.5, height = 7, units = "in", res = 600, restoreConsole = TRUE)

# Largemouth Bass Plot
LMBplot <- ggplot(data1, aes(x = GRP, y = Reservoir, fill = Scenario)) +
              geom_density_ridges2(alpha=.3, scale=0.95) +
              scale_y_discrete(expand = expand_scale(add = c(0.02, 1))) +
              scale_x_continuous(limit = c(-0.013, 0.036)) +
              xlab("GRP (g/g/d)") +
              scale_fill_cyclical(breaks=c("Hypoxia", "No Hypoxia"), 
                                values=c("white", "grey1"),
                                labels = c(expression("With "*italic("f")[c]*"DO", "Without "*italic("f")[c]*"DO")),
                                guide = "legend") +
              theme_ridges() +
              theme_bw() + 
              theme(legend.position = c(0.53, 0.88),
                    legend.title = element_blank(),
                    legend.text = element_text(size = 10),
                    legend.text.align = 0,
                    legend.key.size = unit(.75,"line"),
                    legend.direction = "vertical",
                    legend.spacing.x = unit(0.1, 'cm'),
                    legend.spacing.y = unit(0, "cm"),
                    legend.margin = margin(0.01, 0.01, 0.01, 0.01, "cm"),
                    panel.border = element_blank(), 
                    panel.grid.major = element_blank(),
                    panel.spacing = unit(100, "lines"),
                    text=element_text(size=14),
                    axis.title.y = element_blank(),
                    axis.text.x = element_text(color = "black"),  
                    axis.text.y = element_text(color = "black"),
                    panel.grid.minor = element_blank(),
                    panel.background = element_rect(colour = "black", size = 1)) +
            geom_vline(xintercept = 0, linetype = "dashed") +
            annotate("text", x = -.012, y = 3.75, label = "(A)", size = 4.5)

# Saugeye Plot
SAEplot <- ggplot(data2, aes(x = GRP, y = Reservoir, fill = Scenario)) +
            geom_density_ridges2(alpha=.3, scale=0.95) +
            scale_y_discrete(expand = expand_scale(add = c(0.02, 1))) +
            scale_x_continuous(limit = c(-0.013, 0.036)) +
            xlab("GRP (g/g/d)") +
            scale_fill_cyclical(breaks=c("Hypoxia", "No Hypoxia"), 
                                values=c("white", "grey1")) +
            theme_ridges() +
            theme_bw() + 
            theme(panel.border = element_blank(), 
                  panel.grid.major = element_blank(),
                  panel.spacing = unit(100, "lines"),
                  text=element_text(size=14),
                  axis.title.y = element_blank(),
                  axis.text.x = element_text(color = "black"),  
                  axis.text.y = element_text(color = "black"),
                  panel.grid.minor = element_blank(),
                  panel.background = element_rect(colour = "black", size = 1)) +
            geom_vline(xintercept = 0, linetype = "dashed") +
            annotate("text", x = -.012, y = 3.75, label = "(B)", size = 4.5)

# White Crappie Plot
WHCplot <- ggplot(data3, aes(x = GRP, y = Reservoir, fill = Scenario)) +
            geom_density_ridges2(alpha=.3, scale=0.95) +
            scale_y_discrete(expand = expand_scale(add = c(0.02, 1))) +
            scale_x_continuous(limit = c(-0.013, 0.036)) +
            xlab("GRP (g/g/d)") +
            scale_fill_cyclical(breaks=c("Hypoxia", "No Hypoxia"), 
                                values=c("white", "grey1")) +
            theme_ridges() +
            theme_bw() + 
            theme(panel.border = element_blank(), 
                  panel.grid.major = element_blank(),
                  panel.spacing = unit(100, "lines"),
                  text=element_text(size=14),
                  axis.title.y = element_blank(),
                  axis.text.x = element_text(color = "black"),  
                  axis.text.y = element_text(color = "black"),
                  panel.grid.minor = element_blank(),
                  panel.background = element_rect(colour = "black", size = 1)) +
              geom_vline(xintercept = 0, linetype = "dashed") +
              annotate("text", x = -.012, y = 3.75, label = "(C)", size = 4.5)

# Combine plots
patchwork = LMBplot / SAEplot / WHCplot

# Remove title from second subplot
patchwork[[1]] = patchwork[[1]] + theme(axis.title.x = element_blank(),
                                        axis.text.x = element_blank())

patchwork[[2]] = patchwork[[2]] + theme(axis.text.x = element_blank(),
                                        axis.title.x = element_blank())

# Remove title from third subplot
patchwork[[3]] = patchwork[[3]] + theme(axis.title.y = element_blank())

patchwork

dev.off()
