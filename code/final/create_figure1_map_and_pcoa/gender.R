# Gender
library(reshape2)
library(dplyr)
df$Gender <- df$MEN
d <- melt(df %>% summarise(Male = mean(Gender == "Male"), Female = mean(Gender == "Female")))

p <- ggplot(d, aes(x = 1, y = value, fill = variable)) +
       geom_bar(stat = "identity", position = "dodge") +
       scale_fill_manual(values = c("red", "blue")) +
       scale_y_continuous(
                          breaks = seq(0, 0.5, 0.1),
			  labels = scales::percent(seq(0, 0.5, 0.1))			  
       ) +
       theme(axis.text.y = element_blank()) +       
       guides(fill = FALSE) +
       coord_flip() +
       labs(x = "", y = "", title = "Gender")

figure.gender <- p

#print(p)
#prop.table(table(df$MEN))


