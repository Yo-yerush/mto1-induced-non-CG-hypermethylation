library(ggplot2)
library(cowplot)

gainORloss = function(DMRsReplicates, context) {
  
  gainORloss_vec = DMRsReplicates$regionType
  
  gain_DMRs = length(grep("gain",gainORloss_vec))
  loss_DMRs = length(grep("loss",gainORloss_vec))
  total = gain_DMRs+loss_DMRs
  pres_gain = round((gain_DMRs/total)*100, 1)
  pres_loss = round((loss_DMRs/total)*100, 1)
  
  #pie_data = data.frame(group = c("", paste0("Loss (",pres_loss,"%)"), paste0("Gain (",pres_gain,"%)")),
  #                      value = c(0,loss_DMRs, gain_DMRs))
  pie_data = data.frame(group = c(paste0(pres_gain,"%"), paste0(pres_loss,"%")),
                        value = c(gain_DMRs, loss_DMRs))
  
  # pie plot
  p.start = ifelse(context == "CG", 1.75,
                   ifelse(context == "CHG", 5.4,
                          6))
  pie_plot = ggplot(pie_data, aes(x = "", y = value, fill = group)) +
    geom_bar(stat = "identity", width = 1, color = "white", lwd = 0.5) +
    coord_polar(theta = "y", start = p.start, direction = -1) +
    scale_fill_manual(values = c("#6c96d9", "#d96c6c")) +
    theme_void() +
    theme(axis.line = element_blank(),
          panel.grid = element_blank(),
          panel.border = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "none",
          text = element_text(family = "serif"))
  
  # text of the plot from the rught side
  text_plot = ggplot() +
    geom_text(aes(x=-0.25, y=1), label=pie_data[2,1], size=4) +
    geom_text(aes(x=-0.25, y=-1), label=pie_data[1,1], size=4) +
    xlim(-1, 1) + ylim(-2, 2) +
    theme_void() + theme(text = element_text(family = "serif"))
  
  comb_plot = cowplot::plot_grid(pie_plot, text_plot,
                                  nrow = 1,
                                  rel_widths = c(2,2))

  final_plot = ggdraw(comb_plot) +
    # hyper line
    draw_line(
      x = c(0.410, 0.575),
      y = c(0.725, 0.725),
      color = "black",
      size = 0.25
    ) +
    # hypo line
    draw_line(
      x = c(0.410, 0.575),
      y = c(0.275, 0.275),
      color = "black",
      size = 0.25
    )
  
  svg(file = paste0("pie_",context,"_gainORloss.svg"), width = 2, height = 1, family = "serif")
  par(mar = rep(0,4)) # bottom, left, top, right
  print(final_plot)
  dev.off()
}

setwd("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/mto1_paper/pie_plots")

CG = read.csv("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results/mto1_vs_wt/DMRs_CG_mto1_vs_wt.csv")
CHG = read.csv("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results/mto1_vs_wt/DMRs_CHG_mto1_vs_wt.csv")
CHH = read.csv("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results/mto1_vs_wt/DMRs_CHH_mto1_vs_wt.csv")

gainORloss(CG, "CG")
gainORloss(CHG, "CHG")
gainORloss(CHH, "CHH")



#### old pie plot
#pie_plot = pie(pie_data$value,
#               #pie_data$group,
#               # change colors
#               col = c("#d96c6c", "#6c96d9"),
#               labels = "",
#               border = "white",
#               #main = paste0(total," DMRs in ",context," context"),
#               radius = 1,
#               lwd = 2,
#               clockwise = T)