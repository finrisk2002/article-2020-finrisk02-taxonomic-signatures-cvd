#load("kartta/kartta.RData")
#figure.finland <- figure.finland + theme_bw(10)

#adapted from "MULTIVARIATE DOT-DENSITY MAPS IN R WITH SF & GGPLOT2" https://www.cultureofinsight.com/blog/2018/05/02/2018-04-08-multivariate-dot-density-maps-in-r-with-sf-ggplot2/

library(tidyverse) # dev version of ggplot2 required devtools::install_github('hadley/ggplot2')
library(sf)
library(extrafont)
library(dplyr)
library(readr)
library(purrr)
library(ggplot2)
library(readxl)

#read the number of individuals in each municipality/state
ge_data <- read_excel("numerot3.xlsx")

# shapefile of Finland
uk <- st_read("kuntajako_2017_maa_alueet.shp", stringsAsFactors = FALSE, quiet = TRUE) 

# merge the data
sf_data <- left_join(uk, ge_data) %>% st_as_sf() # I'm losing sf class after join so make sf object again

#create a data frame with the number of dots we want plotted in municipality/state

random_round <- function(x) {
  v=as.integer(x)
  r=x-v
  test=runif(length(r), 0.0, 1.0)
  add=rep(as.integer(0),length(r))
  add[r>test] <- as.integer(1)
  value=v+add
  ifelse(is.na(value) | value<0,0,value)
  return(value)
}

num_dots <- as.data.frame(sf_data) %>% 
  select(Karelia:Lapland) %>% 
  mutate_all(random_round)

#replace NA's in counties without participants as zeros

num_dots[is.na(num_dots)] <- 0

# generates data frame with coordinates for each point + what state it's located in

sf_dots <- map_df(names(num_dots), 
                  ~ st_sample(sf_data, size = num_dots[,.x], type = "random") %>% # generate the points in each polygon
                    st_cast("POINT") %>%                                          # cast the geom set as 'POINT' data
                    st_coordinates() %>%                                          # pull out coordinates into a matrix
                    as_tibble() %>%                                               # convert to tibble
                    setNames(c("lon","lat")) %>%                                  # set column names
                    mutate(Area = .x)                                             # add categorical state
                  ) %>% 
  slice(sample(1:n())) # once map_df binds rows randomise order to avoid bias in plotting order



pal <- c("Karelia" = "#FCB711", "Kuopio" = "#F37021", "Turku" = "#CC004C", "Helsinki" = "#6460AA", "Oulu" = "#0089D0", "Lapland" = "#0DB14B")

# Population counts
counts <- apply(as.matrix(sf_data %>% select(Karelia:Lapland))[, 1:6], 2, function (x) {sum(unlist(x), na.rm = TRUE)})[names(pal)]
sf_dots$Area2 <- sf_dots$Area
for (nam in names(counts)) {
  k <- counts[[nam]]
  newnam <- paste(nam, " (", k, ")", sep = "")
  sf_dots$Area2 <- gsub(nam, newnam, sf_dots$Area2)
  names(pal) <- gsub(nam, newnam, names(pal))  
}
sf_dots$Area2 <- factor(sf_dots$Area2)

# plot it and save as png big enough to avoid over-plotting of the points
p <- ggplot() +
  geom_sf(data = sf_data, fill = "transparent", colour = "black", lwd = 1) +
  geom_point(data = sf_dots, aes(lon, lat, colour = Area2), size=1) +
  scale_colour_manual(values = pal) +
  coord_sf(crs = 3067, datum = NA) +
  labs(x = NULL, y = NULL,
       title = "FINRISK 2002 Participants",
       subtitle = "",
       caption = "") +
  guides(colour = guide_legend(override.aes = list(size = 2), keyheight=0.35, default.unit="inch", title = "")) +
  theme_classic(base_family = "", base_size = 36) +  
  theme(  
        legend.position = c(0.1, 0.6),
        legend.direction = "vertical",
        plot.background = element_rect(fill = "#FFFFFF", color = NA), 
        panel.background = element_rect(fill = "#FFFFFF", color = NA),
        legend.background = element_rect(fill="transparent"), 	
        legend.key = element_rect(fill = "#FFFFFF", colour = NA),
        plot.margin = margin(-2, 0, 0, 0, "cm"),	
        text =  element_text(color = "black")
  ) +
  labs(title = "")

figure.finland <- p

save(figure.finland, file = "kartta.RData")
