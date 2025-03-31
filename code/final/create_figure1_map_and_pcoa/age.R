d <- df
d$Gender <- d$MEN
d$Age <- group_age(d$BL_AGE)


p <- ggplot(d, aes(x = Age, fill = Gender)) +
       geom_bar(alpha = 1,
                color = "black",
		position = "dodge") +
       scale_fill_manual(values = c("red", "blue"),
                         labels = c(
			   #paste0("Female (", round(100 * mean(d$Gender == "Female"), 1) , "%)"),
			   #paste0("Male (", round(100 * mean(d$Gender == "Male"), 1) , "%)")
			   "Female",
			   "Male"
			            )) +
       theme(legend.position = c(0.22, 0.95),
                  legend.text = element_text(size = 0.9 * basesize), 
                  legend.title = element_text(size = basesize), 
                  legend.key.size=unit(1, "line"),
          legend.background = element_rect(fill="transparent"), 		  
		  axis.text.x = element_text(size = basesize),
		  axis.text.y = element_text(size = basesize),
		  axis.title.x = element_text(size = basesize),
		  axis.title.y = element_text(size = basesize)
		  ) +
       guides(fill = guide_legend(title = "")) + 
       labs(x = "Age (y)", y = "Frequency (n)", title = "") +
       scale_x_discrete(
                          labels = as.character(seq(20, 70, 10))
       		       ) 
    #annotate("text", x = 5.3, y = 950,
              #label = paste0("Mean Age:\n",#
	                     #round(mean(df$BL_AGE, na.rm = TRUE), 0),
			     #" y"
			     #(",
			     #round(min(df$BMI, na.rm = TRUE), 0),
			     #"-",
			     #round(max(df$BMI, na.rm = TRUE), 0),
			     #")"
		             # 	     ),
              #size = .4 * basesize,
	      #color = "black"
       
# print(p)

figure.age <- p



