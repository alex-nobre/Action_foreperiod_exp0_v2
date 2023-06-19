

# General theme
genTheme <- theme_classic +
  theme(plot.title = element_text(size = 16, hjust = 0.5))


# Colors for factors
condColors <- scale_color_manual(values = c("orange", "blue"))
oneBackFPColors <- scale_color_manual(values = c("blue", "orange", "green", "magenta"))