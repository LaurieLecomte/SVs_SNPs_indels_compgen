library(rgdal)
library(ggplot2)


canada_poly <- "maps/CAN_adm0.shp"

water_lines <- "maps/watrcrsl.shp"

area<-readOGR(canada_poly)
print('imported can poly')

az_fortify <- fortify(area)

#ggplot(az_fortify)+geom_polygon(aes(long,lat,group=group),colour="grey90")+
#  scale_x_continuous(limits = c(36.928662, 31.347814))+
#  scale_y_continuous(limits = c(-114.293736, -109.057855))
QC_sub <- subset(az_fortify, lat < 65 & lat > 40 & long < -50 & long > -85)
  
HSP <- c(50.233333, -63.6)

base_plot <- 
ggplot(QC_sub) + 
  geom_polygon(aes(long, lat, group = group), col = "black", fill = "grey90") + 
  geom_point(aes(y = HSP[1], x = HSP[2]), col = "red", size = 10) +
  theme_void() +
  geom_text(x = c(-72), y = 54, label = "Quebec", size = 10)

ggsave(plot = base_plot,
       filename = "maps/RO_bassin_on_QC.png", 
       width = 180, 
       height = 200,
       units = "mm",
       dpi = 600)

# Romaine
lower_left <- c(y = 50.238172, x = -63.926478)
upper_left <- c(y = 50.238172, x = -63.926478)
lower_right <- c(y = 50.238172, x = -63.237298)
upper_right <- c(y = 50.425819, x = -63.237298)

RO_PU <- subset(az_fortify, lat < 50.425819 & lat > 50.238172 & long < -63.237298 & long > -63.926478)

plot_RO_bassin <- 
ggplot(RO_PU) + 
  geom_path(aes(long,lat, group=group), colour="black", fill = NA)+ 
  #coord_cartesian(ylim = c(40, 65), xlim = c(-85, -50)) +
 # geom_point(aes(y = HSP[1], x = HSP[2]), col = "red", size = 3) +
  theme_void()

ggsave(plot = plot_RO_bassin,
       filename = "maps/RO_bassin.png", 
       width = 180, 
       height = 200,
       units = "mm",
       dpi = 600)



water <- readOGR(water_lines)
water_fortify <- fortify(water)

RO_PU_water <- subset(water_fortify, lat < 50.425819 & lat > 50.238172 & long < -63.237298 & long > -63.926478)

plot_water_bassin <-
ggplot(RO_PU_water) + 
  geom_polygon(aes(long,lat, group=group), colour="blue") 
  #coord_cartesian(ylim = c(50.238, 50.425), xlim = c(-63.926, -63.237))
  
ggsave(plot = plot_water_bassin,
       filename = "maps/water_bassin.png", 
       width = 180, 
       height = 200,
       units = "mm",
       dpi = 600)

#plot_water_
#ggplot(RO_PU_water) + 
#  geom_polygon(aes(long,lat, group = group), fill = "blue")


#canada_lines <- "C:/Users/15818/Documents/M.Sc/Bioinfo/carto/sde-columbia-iscgm_canada_2007_polbndl-shapefile/columbia_iscgm_canada_2007_polbndl.shp"
#can <-readOGR(canada_lines)
#can_fortify <- fortify(can)

#ggplot(can_fortify) + 
#  geom_polygon(aes(long,lat, group=group), colour="grey70", fill = NA)

#  coord_cartesian(ylim = c(50.238, 50.425), xlim = c(-63.926, -63.237))
