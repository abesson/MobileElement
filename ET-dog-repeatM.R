# statistics on Mobile Elements from dog genome using script:
#https://bioinfo-fr.net/elements-repetes-genome-humain-apercu-rapide-r-tidyverse
# with repeatMasker.alldata: repeatMaskerCanFam3.1-300418.txt


# libraries
library(rlang)
library(purrr)
library(dplyr)
library(readr)
library(ggplot2)
library(cowplot)
library(tidyverse)
library(forcats)
library(svglite)
library(viridis)
library(xtable)
library(broom)

# export row data 
repeats_dirty = read_tsv("repeatMaskerCanFam3.1-300418.txt")

# focus on chr1 to 38 and sexual to keep only standard chrm
standard_chrm = c(paste0("chr", 1:38), "chrX", "chrY")

# creation of a new table containing only: chr, start, end, strand, name, family and class
repeats = select(repeats_dirty,genoName, genoStart,genoEnd,strand,repName, repFamily,repClass) %>%
  dplyr::rename(chr = genoName, start=genoStart,end=genoEnd,name=repName,family=repFamily, class=repClass) %>%
  filter(chr %in% standard_chrm)

# filter lines with "?" corresponding to uncertain classification
#98% kept with chr and 97,5% kept without ? from original file
repeats=filter(repeats,!(grepl("?",repeats$class, fixed = TRUE) | grepl("?",repeats$family, fixed = TRUE) | grepl("?",repeats$name, fixed = TRUE)))


### TE class ###
# Reorder levels by first frequency
repeats$class = factor(repeats$class %>% fct_infreq %>% fct_rev)

# first bar graph containing the quantity of elements by class
plot_Class=ggplot(repeats,aes(x = class)) +
  geom_bar(stat = "count", fill = "indianred1") +
  geom_text(aes(label = ..count..), y = 10000, hjust = 0, stat = "count") +
  labs(x = "Class", y = "Number of elements") +
  scale_y_continuous(sec.axis = dup_axis()) +
  coord_flip() +
  background_grid(major = "xy", minor = "none")

# second bar graph containing the length of each class into the genome
# length of each chromosome
chr_length = read_tsv("canFam3.sizes", col_names = FALSE) %>%
  dplyr::rename(seqnames = X1, length = X2) %>%
  filter(seqnames %in% standard_chrm) %>%
  mutate(seqnames = factor(seqnames, levels = standard_chrm))
genome_length = sum(as.numeric(chr_length$length))

# new column with width of each element
repeats = mutate(repeats, width = end - start)

# table by TE class (total length for each class)
plot_genomeProp = repeats %>%
  group_by(class) %>%
  summarise(total_length = sum(width)) %>%
  arrange(total_length) %>%
  ggplot(aes(x = class, y = total_length)) +
  geom_bar(stat = "identity", fill = "cornflowerblue") +
  geom_text(aes(label = paste0(round(100 * total_length / genome_length, digits = 1), "%")), y = 10000, hjust = 0) +
  labs(x = "Class", y = "Accumulated length (bp)") +
  scale_y_continuous(sec.axis = sec_axis(~100 * . / genome_length, name = "% of genome")) +
  coord_flip() +
  theme(axis.text.y = element_blank(), axis.title.y = element_blank()) +
  background_grid(major = "xy", minor = "none")

# third boxplot graph containing the distribution size for each TE class
plot_sizeDistrib = ggplot(repeats, aes(x = class, y = width)) +
  geom_boxplot(fill = "mediumorchid", outlier.shape = NA) +
  labs(x = "Class", y = "Element size (bp)") +
  scale_y_continuous(sec.axis = dup_axis()) +
  coord_flip(ylim = c(0, 2500)) +
  theme(axis.text.y = element_blank(), axis.title.y = element_blank()) +
  background_grid(major = "xy", minor = "none")

# graph reorganisation and svg file creation
myoffset = 0.008
firstplotsize = 0.44
svglite("TEclass-dog.svg", width = 10, height = 5)
ggdraw() +
  draw_plot(plot_Class, x = 0.0, y = myoffset, w = firstplotsize, h = 0.96 - 2*myoffset) +
  draw_plot(plot_genomeProp, x = firstplotsize, y = 0.0, w = (1-firstplotsize)/2, h = 0.96) +
  draw_plot(plot_sizeDistrib, x = firstplotsize + (1-firstplotsize)/2, y = 0.0, w = (1-firstplotsize)/2, h = 0.96) +
  draw_plot_label(LETTERS[1:3], x = c(0.14, firstplotsize - 0.01, firstplotsize + (1-firstplotsize)/2 - 0.01), y = 0.92) +
  draw_label("Transposable Element Class", size = 15, x = 0.5, y = 0.97)
dev.off()


### % repeat elements in dog genome about 42% ###
sum(repeats$width) / genome_length * 100


### Mobile element distribution for each dog chromosome ###
# Focus on the majority class (LINE, SINE, LTR)
svglite("TEbyChrom.svg", width = 5, height = 3)
repeats %>%
  mutate(simple_class = if_else(
    !(class %in% c("SINE", "LINE", "LTR")),
    "other",
    as.character(class)
  )) %>%
  group_by(chr, simple_class) %>% # group by chromosome and by element class
  summarize(cum_size = sum(width)) %>%
  left_join(chr_length, by = c("chr" = "seqnames")) %>%
  ungroup %>%
  mutate(
    frac_repeat = cum_size/length,
    simple_class = factor(simple_class, levels = c("other", "LTR", "SINE", "LINE")),
    chr = factor(chr, levels = standard_chrm)
  ) %>%
  ggplot(aes(x = chr, y = frac_repeat, fill = simple_class)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5)) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c("gray", viridis(3, begin = 0.25, end = 0.9))) +
  labs(x = NULL, y = "% of chromosome", title = "Repeated content of each chromosome", fill = "Class") +
  background_grid(major = "xy", minor = "none")
dev.off()


### repeated element families ###
# level of classification (family, subfamily) -type latex or html- 
# to determine the number of families by class element
group_by(repeats, class) %>%
  summarise(n_family = length(unique(family)), n_subfamily = length(unique(name)), n_element = n()) %>%
  arrange(desc(n_family), desc(n_subfamily)) %>%
  xtable() %>%
  print(type = "latex", include.rownames = FALSE)
  
# data preparation
effectif_table = group_by(repeats, class, family) %>%
  summarise(diff_name = length(unique(name)), size = n()) %>%
  ungroup()

# TE containing some families
multiFamily = c("LINE", "SINE", "DNA", "LTR")

# figure panel for element class containing multiple families
myFamilyPlots = map(
  multiFamily,
  ~filter(effectif_table, class == .x) %>%
    arrange(size) %>%
    mutate(family = factor(family, levels = family)) %>%
    mutate(text_dark = if_else(diff_name > 50, TRUE, FALSE)) %>%
    ggplot(aes(x = family, y = size, fill = diff_name)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = size, color = text_dark), y = 1, hjust = 0) +
    scale_fill_viridis(limits = c(0,100), oob = scales::squish) +
    scale_color_manual(values = c("grey80", "black"), guide = FALSE) +
    scale_y_log10() +
    coord_flip(ylim = c(10,1e7)) +
    background_grid(major = "xy", minor = "none") +
    labs(x = NULL, y = NULL, title = .x) +
    theme(legend.position = "none")
)
names(myFamilyPlots) = multiFamily
myFamilyPlots$DNA = myFamilyPlots$DNA + ylab("Number of elements")

# figure panel for TE with only 1 family
plot_other = filter(effectif_table, !(class %in% c("LINE", "SINE", "DNA", "LTR"))) %>%
  arrange(size) %>%
  mutate(family = factor(family, levels = family)) %>%
  ggplot(aes(x = family, y = size, fill = diff_name)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = size), y = 1, hjust = 0, stat = "count") +
  scale_fill_viridis(limits = c(0, 100), oob = scales::squish) +
  scale_y_log10() +
  coord_flip(ylim = c(10, 1e7)) +
  background_grid(major = "xy", minor = "none") +
  labs(x = NULL, y = NULL, title = "Other", fill = "Quantity of subfamilies") +
  theme(legend.position = "bottom")

# legend recuperation
myLegend = get_legend(plot_other)
plot_other = plot_other + theme(legend.position = "none")

# homogeneisation (black magic)
myFamilyPlots = map(myFamilyPlots, ggplotGrob)
plot_other = ggplotGrob(plot_other)
myFamilyPlots.widths = map(myFamilyPlots, ~.x$widths[1:3])
plot_other.widths = plot_other$widths[1:3]

max.widths = grid::unit.pmax(
  plot_other.widths,
  do.call(grid::unit.pmax, myFamilyPlots.widths)
)

plot_other$widths[1:3] = max.widths
myFamilyPlots = map(myFamilyPlots, function(x) {
  x$widths[1:3] = max.widths
  return(x)
})

# panels rearrangement
svglite("familleER.svg", width = 10, height = 7)
ggdraw() +
  draw_text("Families of repeated elements", x = 1/6, y = 0.95, size = 18) +
  draw_plot(myLegend, x = 0.0, y = 0.75, w = 1/3, h = 0.2) +
  draw_plot(myFamilyPlots$DNA, x = 0.0, y = 0.0, w = 1/3, h = 0.8 ) +
  draw_plot(myFamilyPlots$LINE, x = 1/3, y = 0.62, w = 1/3, h = 0.38 ) +
  draw_plot(myFamilyPlots$LTR, x = 1/3, y = 0.31, w = 1/3, h = 0.31 ) +
  draw_plot(myFamilyPlots$SINE, x = 1/3, y = 0.0, w = 1/3, h = 0.31 ) +
  draw_plot(plot_other, x = 2/3, y = 0.0, w = 1/3, h = 0.7 )
dev.off()

# Note: the SINE's family "tRNA-Lys" represents the SINEC family (verification on the repeatMasker file)


### repeats elements families : only SINE and LINE ###

# data preparation
effectif_table = group_by(repeats, class, family) %>%
  summarise(diff_name = length(unique(name)), size = n()) %>%
  ungroup()

# figure panel by TE containing some families
multiFamily = c("LINE", "SINE")

myFamilyPlots = map(
  multiFamily,
  ~filter(effectif_table, class == .x) %>%
    arrange(size) %>%
    mutate(family = factor(family, levels = family)) %>%
    mutate(text_dark = if_else(diff_name > 50, TRUE, FALSE)) %>%
    ggplot(aes(x = family, y = size, fill = diff_name)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = size, color = text_dark), y = 1, hjust = 0) +
    scale_fill_viridis(limits = c(0,100), oob = scales::squish) +
    scale_color_manual(values = c("grey80", "black"), guide = FALSE) +
    scale_y_log10() +
    coord_flip(ylim = c(10,1e7)) +
    background_grid(major = "xy", minor = "none") +
    labs(x = NULL, y = NULL, title = .x) +
    theme(legend.position = "none")
)
names(myFamilyPlots) = multiFamily
myFamilyPlots$SINE = myFamilyPlots$SINE + labs(y = "Nombre d'éléments")

# figure panel for TE with only 1 family
plot_other = filter(effectif_table, !(class %in% c("LINE", "SINE"))) %>%
  arrange(size) %>%
  mutate(family = factor(family, levels = family)) %>%
  ggplot(aes(x = family, y = size, fill = diff_name)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = size), y = 1, hjust = 0, stat = "count") +
  scale_fill_viridis(limits = c(0, 100), oob = scales::squish) +
  scale_y_log10() +
  coord_flip(ylim = c(10, 1e7)) +
  background_grid(major = "xy", minor = "none") +
  labs(x = NULL, y = NULL, title = "Autres", fill = "Quantity of subfamilies") +
  theme(legend.position = "bottom")

# legend recuperation
myLegend = get_legend(plot_other)
plot_other = plot_other + theme(legend.position = "none")

# homogeneisation (black magic)
myFamilyPlots = map(myFamilyPlots, ggplotGrob)
plot_other = ggplotGrob(plot_other)
myFamilyPlots.widths = map(myFamilyPlots, ~.x$widths[1:3])
plot_other.widths = plot_other$widths[1:3]

max.widths = grid::unit.pmax(
  plot_other.widths,
  do.call(grid::unit.pmax, myFamilyPlots.widths)
)

plot_other$widths[1:3] = max.widths
myFamilyPlots = map(myFamilyPlots, function(x) {
  x$widths[1:3] = max.widths
  return(x)
})

# panels rearrangement
svglite("familleER.svg", width = 10, height = 7)
ggdraw() +
  draw_plot(myLegend, x = 2/3, y = 0.5, w = 1/3, h = 0.4) +
  draw_plot(myFamilyPlots$SINE, x = 0.0, y = 0.5, w = 1/3, h = 0.4 ) +
  draw_plot(myFamilyPlots$LINE, x = 1/3, y = 0.5, w = 1/3, h = 0.4 ) 

dev.off()

# Note: the SINE's family "tRNA-Lys" represents the SINEC family (vérification on the repeatMasker file)
