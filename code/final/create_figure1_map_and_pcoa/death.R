d <- df
d$Age <- group_age(d$BL_AGE)

# DEATH, DEATH_AGE, DEATH_AGEDIFF

# > table(d$DEATH)
#   0    1 
# 6448  667 
# > prop.table(table(d$DEATH))
#         0          1
# 0.90625439 0.09374561 
# 667 deaths (9.4%) during the followup

dd <- d %>% filter(DEATH == 1) %>%
            mutate(time = floor(DEATH_AGEDIFF)) %>%
            group_by(Gender, time) %>%
      	    summarise(n = n()) %>%
	    mutate(incidences = cumsum(n))

p <- ggplot(dd, aes(x = time, y = incidences, color = Gender)) +
       geom_line() +
       geom_point() +
       scale_color_manual(values = c("red", "blue")) +
       scale_x_continuous(limits = c(0, 15)) +       
       labs(x = "Time (y)", y = "Deaths (n)", title = "") +
       guides(fill = "none", color = "none") +
       theme(
		  axis.text.x = element_text(size = basesize),
		  axis.text.y = element_text(size = basesize),
		  axis.title.x = element_text(size = basesize),
		  axis.title.y = element_text(size = basesize)
		  ) 
    #geom_text(x = 2, y = 400,
    #          label = paste0("Deaths:\n",
    #	                     sum(df$DEATH == 1, na.rm = TRUE)),
    #          size = .4 * basesize,
    #	      color = "black")
    #annotate("text", x = 2, y = 380,
     #         label = paste0("Deaths:\n",
#	                     sum(df$DEATH == 1, na.rm = TRUE)),
 #             size = .4 * basesize,
#	      color = "black") 
       #print(p)       

figure.death <- p



