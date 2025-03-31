# Load the FINRISK data set

library(microbiome)
pseq <- phfinrisk_genus
df <- meta(pseq)
df$Gender <- df$MEN

# Statistics
tbl <- NULL
tbl <- rbind(tbl, c("Sample size", nrow(df)))
tbl <- rbind(tbl, c("Female/Male", paste(round(100 * prop.table(table(df$MEN)), 1)[["Female"]],
                                         round(100 * prop.table(table(df$MEN)), 1)[["Male"]], sep = "/")))

# Read count
pdf("readcount.pdf"); hist(log10(colSums(abundances(pseq))), 100, main = "Read count"); dev.off()

basesize <- 30
theme_set(theme_classic(basesize))

library(cowplot)
theme_set(theme_cowplot())

source("gender.R")
source("age.R") # BL_AGE
source("bmi.R") # BMI
source("death.R") # DEATH, DEATH_AGE, DEATH_AGEDIFF
source("core.R")


fig <- plot_grid(figure.age,
                 figure.bmi,
		 figure.death,
		 nrow = 3)

theme_set(theme_classic(basesize))
figh <- plot_grid(figure.age,
                 figure.bmi,
		 figure.death,
		 nrow = 1) 
# print(fig)
fig1A <- figh

