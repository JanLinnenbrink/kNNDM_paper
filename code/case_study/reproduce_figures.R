
outtab <- read.csv(file.path(outfolder, "outtab.csv"))
outtab$design <- factor(outtab$design, levels=c("regular", "simpleRandom", "clusterMedium", 
                                                "clusterStrong", "clusterGapped"))
outtab$variate <- as.factor(outtab$variate)

outtab$method <- factor(outtab$method, levels = c("random","spatial","heteroscedastic","intensity",
                                                  "modelbased","knndm_sample","exhaustive"))

collabs <- c("random","regular","clustMed","clustStr","clustGap")

#
library("ggplot2")
library("ggthemes")
library("cowplot")

lw <- 0.6
w <- 0.45
base.size <- 15

(diff_rmse <- ggplot(outtab[outtab$methodID!=0,],
                    aes(x=design, y=rRMSE,colour=method)) +
  geom_boxplot(fill=NA, lwd=lw, width=w,position=position_dodge(0.6)) +
  scale_color_colorblind(name="CV method") +
  geom_hline(aes(yintercept=0)) +
  xlab("") +
  ylab(expression(CV - true~RMSE)) +
  theme_bw(base_size = base.size) +
  theme(legend.position = "None") +
  facet_wrap(~variate) )

(diff_mec <- ggplot(outtab[outtab$methodID!=0,],
                     aes(x=design, y=rMEC,colour=method)) +
    geom_boxplot(fill=NA, lwd=lw, width=w,position=position_dodge(0.6)) +
    scale_color_colorblind(name="CV method") +
    geom_hline(aes(yintercept=0)) +
    xlab("") +
    ylab(expression(CV - true~MEC)) +
    theme_bw(base_size = base.size) +
    theme(legend.position = "None") +
    facet_wrap(~variate) )

legend_bottom <- get_legend(diff_mec +
                             guides(color = guide_legend(nrow = 1)) +
                             theme(legend.position = "bottom"))

pgr <- plot_grid(diff_rmse, diff_mec, legend_bottom,
                 nrow=3, rel_heights = c(1,1,0.5), axis="l", align = "v")
ggsave("figures/4_results_sim_deBruin.pdf", pgr, height=8, width=8)
