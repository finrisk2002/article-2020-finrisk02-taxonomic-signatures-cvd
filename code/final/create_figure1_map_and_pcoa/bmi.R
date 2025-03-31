d <- df
d$Age <- group_age(d$BL_AGE)
d$BMI_group <- group_bmi(d$BMI)

# BMI Missing for N=2 subjects
# sum(is.na(d$BMI))
library(latex2exp)
brs <- c(0,2,4,6,8,10)/100
p <- ggplot(d, aes(x = BMI, fill = Gender)) +
       geom_density(alpha = .5, color = "black") +
       # geom_bar(alpha = 1, color = "black", position = "dodge") +
       scale_fill_manual(values = c("red", "blue")) +
       guides(fill = "none") +
       theme(
		  axis.text.x = element_text(size = basesize),
		  axis.text.y = element_text(size = basesize),
		  axis.title.x = element_text(size = basesize),
		  axis.title.y = element_text(size = basesize)
		  ) +

       labs(x = TeX("BMI ($kg/m^2$)"),
            y = "Frequency (%)", title = "") +
       scale_x_continuous(
			  breaks = seq(20, 50, 10),
                          labels = seq(20, 50, 10)
       		       )  +
       #scale_y_continuous(breaks = c(0,2,4,6,8,10)/100, label = scales::percent)
       #scale_y_continuous(breaks = brs, label = paste0(100 * brs, "%"))
       scale_y_continuous(breaks = brs, label = 100 * brs)       
    #annotate("text", x = 45, y = 0.1,
    #          label = paste0("Mean BMI:\n",#
    #    	                     round(mean(df$BMI, na.rm = TRUE), 0),
#			     " kg/m2"
#			     #(",
#			     #round(min(df$BMI, na.rm = TRUE), 0),
#			     #"-",
#			     #round(max(df$BMI, na.rm = TRUE), 0),
#			     #")"
#			     ),
 #             size = .4 * basesize,
#	      color = "black")        
       

#print(p)
figure.bmi <- p



