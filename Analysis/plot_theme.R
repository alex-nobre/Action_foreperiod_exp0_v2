#textwidth = 17.5 # cm

# sizes (linewidth, point size)
lsz = 1.2
psz = 1.0
erlsz = 5*lsz/6
erlwth = lsz/12

# jitter spacings
jtwdt = 0.15
jthgt = 0

# alpha
alp =  0.3

# Condition colors
cond_cols <- scale_color_manual(values = c("#48f480", "#F448BD"), labels = c("External", "Action"))

# define a 'theme'
mytheme	 <- theme_classic() + theme(
  # Set font
  text = element_text(family = "Microsoft Sans Serif"),
  # Set axis text size
  axis.text = element_text(size = 6),
  axis.title = element_text(size = 6),
  # Set title elements
  plot.title = element_text(size = 8, hjust = 0.5),
  plot.tag = element_text(size = 8),
  strip.text = element_text(size = 8),
  strip.background = element_blank(),
  # Set legend text size
  legend.title = element_text(size = 6),
  legend.text = element_text(size = 6),
  # Set legend spacing and position
  legend.box.spacing = unit(0, "pt"),
  legend.position = "bottom",
  legend.key = element_blank(),
  # Remove grid and background
  panel.grid = element_blank(),
  panel.background = element_blank(),
  panel.border = element_rect(color = "black", fill = NA),
  axis.line = element_blank())
