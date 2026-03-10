library(ggplot2)
library(dplyr)
library(geomtextpath)

#var1 = "EV"
#var2 = "dCGS"

var1 = "wt"
var2 = "mto1"

for (context in c("CG","CHG","CHH")) {
  
  genome_ann_path = paste0("C:/Users/yonye/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results_040424/results/",var2,"_vs_",var1,"/genome_annotation/")
  
  ###### prepare plot df
  plot_levels_loop = c("Promoters","CDS","Introns","fiveUTRs","threeUTRs","TEG")#,"Transposable_Elements")
  ann_plot_final_df = data.frame(ann = plot_levels_loop, total = NA, gain = NA)
  
  for (ann.loop in plot_levels_loop) {
    DMRsReplicates_loop_path = paste0(genome_ann_path,"/",context,"/",ann.loop,"_",context,"_genom_annotations.csv")
    
    if (file.exists(DMRsReplicates_loop_path)) {
      
      ann_plot_vec = read.csv(DMRsReplicates_loop_path)$regionType
      
    } else {
      ann_plot_vec = data.frame()
    }
    
    ann_plot_final_df[grep(ann.loop,ann_plot_final_df$ann),2] = length(ann_plot_vec)
    ann_plot_final_df[grep(ann.loop,ann_plot_final_df$ann),3] = length(ann_plot_vec[ann_plot_vec == "gain"])
  }
  
  ann_plot_final_df$percent = round((ann_plot_final_df$gain/ann_plot_final_df$total)*100, digits = 1)
  ann_plot_final_df$percent[ann_plot_final_df$percent == 0] = 100
  ann_plot_final_df$percent_loss = round((100-ann_plot_final_df$percent), digits = 1)
  
  plot_levels = c("Promoters","CDS","Introns","5'UTRs","3'UTRs","TEG")#, "TE")
  ann_plot_final_df$ann = plot_levels
  #########################
  y_lim_total = ann_plot_final_df$total
  y_lim_max = max(y_lim_total)
#  y_lim_min = y_lim_max*0.2
#for (YL in 1:length(y_lim_total)) {
#  y_lim_total[YL] = ifelse(y_lim_total[YL] < y_lim_min, y_lim_max*0.9, y_lim_total[YL]*0.85)
#}

  #########################
  # plot
  svg(file = paste0(genome_ann_path,context,"_genom_annotations.svg"), width = 2.45, height = 2, family = "serif")
  
  ann_plot = ann_plot_final_df %>% ggplot() + geom_col(
    aes(x = factor(ann, level=plot_levels), y = total, fill = -percent), # -percent beacouse of the color scale
    colour="black",
    position = "dodge2"
    #alpha = .5
    ) +
    #scale_fill_gradient(low="mediumblue", high="tomato", limits=c(0,100), breaks=seq(0,100,by=10)) +
    scale_fill_gradient2("Regions type\ndistribution", midpoint = -50, low = "#d97777", mid = "#FFFFFF", high = "#7676d6", 
                         limits=c(-100,0), breaks=seq(-100,0,by=50),
                         labels = c("100% Gain","Neutral","100% Loss")) +
    #guides(fill = guide_colourbar()) +
    guides(fill = guide_legend()) +
    #guides(fill = guide_legend(override.aes = list(alpha = .5))) +
    ylim(0,max(ann_plot_final_df$total)*1.15) +
    #coord_polar() + 
    coord_curvedpolar(start = 0) +
    theme(
      panel.background = element_rect(fill = "white", color = "white"),
      panel.grid = element_blank(),
      panel.grid.major.x = element_blank(),
      axis.title.y = element_text(size = 10, face = "bold"),#, margin = margin(t = 0, r = 20, b = 0, l = 0)),
      axis.text.x = element_text(size = 10, face="bold"),
      axis.text.y = element_text(size = 6, face="bold"),#element_blank(),
      title = element_text(size = 8, face = "bold"),
      panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5),
      #legend.title=element_text(size = 8, face = "bold"),
      #legend.text = element_text(size = 8, face = "bold"),
      legend.position = "none") +
    geom_text(aes(x = ann, y = y_lim_max*0.8, label = total), size = 2.65) + 
    labs(title = paste0(var2," vs ",var1," - ",context," context"),
         x = element_blank(),
         y = "Number of DMRs")
  plot(ann_plot)
  
  dev.off()
  
}

###################
### legend
colfunc_ann = list(colfunc_up <- colorRampPalette(c("#d96868", "#FFFFFF")),
                   colfunc_down <- colorRampPalette(c("#6969db", "#FFFFFF")))
legend_ann <- as.raster(matrix(c(colfunc_ann[[1]](20), colfunc_ann[[2]](20)[20:1]), ncol=1))

svg(file = paste0(genome_ann_path,"legend_genom_annotations.svg"), width = 1.34, height = 1.83, family = "serif")
par(mar = c(0,0,2,0))
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = "Regions type\ndistribution")
text(x=0.7, y = seq(0.15,0.85,l=5), labels = c("100% Loss","","50%","","100% Gain"), adj = 0)#c(-2,"",0,"",2))#seq(-2,2,l=5))
rasterImage(legend_ann, 0.2, 0.1, 0.6, 0.9) # xleft, ybottom, xright, ytop
rect(0.2, 0.1, 0.6, 0.9, border="black", lwd = 2)   # xleft, ybottom, xright, ytop
dev.off()
