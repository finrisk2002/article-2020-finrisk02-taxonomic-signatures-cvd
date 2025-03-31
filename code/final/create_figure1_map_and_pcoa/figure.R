library(cowplot)
fig1A <- fig1A #+ annotate("text", x = .02, y = .94, label = "A", size = 20)
fig1B <- figure.finland# + annotate("text", x = .02, y = 2, label = "B", size = 20)
fig1C <- figure.pcoa #+ annotate("text", x = .02, y = .95, label = "C", size = 20)

low.panel <- plot_grid(fig1B, 
                       fig1C,
		       nrow = 1,
		       rel_widths = c(2,3))

# theme_set(theme_classic(20))
fig1 <- plot_grid(
          fig1A,
	  low.panel,
	  rel_heights = c(7, 12),
	  ncol = 1
	  ) +
	  annotate("text", x = .03, y = .97, label = "A", size = 14) +
	  annotate("text", x = .03, y = .58, label = "B", size = 14) +
	  annotate("text", x = .52,  y = .58, label = "C", size = 14) 

pdf("Fig1.pdf", width=20, height=15); print(fig1); dev.off()
jpeg("Fig1.jpg", width=1400, height=1000); print(fig1); dev.off()

