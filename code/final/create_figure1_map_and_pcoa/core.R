figure.core <- plot_core(core(transform(pseq, "compositional"), detection = .1/100, prevalence = 50/100), plot.type = "heatmap", prevalences = seq(0.1, 1, 0.1), detections = 20)

# print(figure.core)
jpeg("figure_core.jpeg")
print(figure.core)
dev.off()
