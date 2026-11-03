library(ggplot2)
icos <- read.csv("outputs/icos_metadata.csv")
fluxn <- read.csv("outputs/fluxnet_metadata.csv")


# ==============================================================================
# one row per station, with the archive(s) it belongs to
sites_all <- rbind(icos[, c("site", "igbp", "lat", "lon")], fluxn[, c("site", "igbp", "lat", "lon")])
sites_all <- sites_all[!duplicated(sites_all$site), ]          # keep the first row (ICOS comes first)
in_icos <- sites_all$site %in% icos$site
in_fx   <- sites_all$site %in% fluxn$site
sites_all$archive <- ifelse(in_icos & in_fx, "Both", ifelse(in_icos, "ICOS", "FLUXNET"))

eu    <- subset(sites_all, lon > -25 & lon < 35 & lat > 34 & lat < 72)   # Europe only
world <- map_data("world")


# ==============================================================================
p1 <- ggplot() +
  geom_polygon(data = world, aes(long, lat, group = group), fill = "grey92", colour = "grey60", linewidth = 0.2) +
  geom_point(data = eu, aes(lon, lat, shape = archive, fill = igbp), size = 2.5) +
  scale_shape_manual(values = c(Both = 21, FLUXNET = 24, ICOS = 23)) +   # circle, triangle, diamond
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  coord_quickmap(xlim = c(-12, 35), ylim = c(35, 71)) +
  labs(x = NULL, y = NULL, shape = "Archive", fill = "Forest type") +
  theme_bw()
p1


eu_icos <- subset(icos, lon > -25 & lon < 35 & lat > 34 & lat < 72)


# ==============================================================================
p2 <- ggplot() +
  geom_polygon(data = world, aes(long, lat, group = group), fill = "grey92", colour = "grey60", linewidth = 0.2) +
  geom_point(data = eu_icos, aes(lon, lat, fill = igbp), shape = 23, size = 2.5) +
  geom_text(data = eu_icos, aes(lon, lat, label = site), size = 2.2, hjust = -0.25) +
  coord_quickmap(xlim = c(-12, 35), ylim = c(35, 71)) +
  labs(x = NULL, y = NULL, fill = "Forest type") +
  theme_bw()
p2

# ==============================================================================
ggsave("outputs/map_icos_fluxnet.png", p1, width = 7, height = 7, dpi = 300)
ggsave("outputs/map_icos.png",         p2, width = 7, height = 7, dpi = 300)



