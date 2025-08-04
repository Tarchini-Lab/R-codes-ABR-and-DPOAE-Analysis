##ABR traces pipeline
#Elli Hartig 05/06/22

#load libraries
library(tidyverse)
library(ggpubr)
library(Rmisc)
library(dplyr)


#import data from BioSig exports. This script pulls all csv files from a folder and compiles them
setwd("H:/")
getwd()
traces.raw <- list.files(path ="conditional_16kHz", full.names = TRUE) %>% 
  lapply(read_csv) %>% 
  bind_rows()

#select only useful columns and create a "tibble" dataframe
traces.values <- traces.raw[,c(8,14,49:292)]
traces.df <- as_tibble(traces.values)

#use gather function to reorganize / "tidy" data 
all.sets <- as.factor(c(0:243))

traces.gather <-  traces.values %>%  gather(all.sets, key = "set", value = "uv")

#create a variable for time. Each time point in Biosig export is 0.04 ms)

traces.time <- traces.gather %>%
  mutate(time = as.numeric(set) * 0.04)
traces.time$`Sub. ID`
#add variable for genotype. First define which Sub IDs are your wildtypes, then use an if/else statement and mutate function to create a new column based on genotype info

wildtypes <- c("AGG29", "AGG33", "AGG34", "AGG35", "AGG36", "AGG37", "AGG38", "AGG39", "AGG40", "AGG41", "AGG44", "AGG46", "AGG49", "AGG50", "AGG51", "AGG52", "AGG54", "AGG55", "AGG57", "AGG58", "AGG59", "AGG60", "AGG61", "AGG73", "AGG74", "AGG74", "AGG77", "AGG78")
traces.geno <- traces.time %>% mutate(genotype = if_else(`Sub. ID` %in% wildtypes, "wt", "mut"))

#settings for graph aesthetics
traces.geno$`Level(dB)` <- factor(traces.geno$`Level(dB)`, levels= c("90","85","80","75","70","65","60", "55", "50", "45", "40", "35", "30","25","20"))

traces.cutoff <- traces.geno %>%  filter(time < 7.5)

cbPalette <- c("blue", "red")

spl.names <- c(
  `90` = "90 dB",
  `85` = "85 dB",
  `80` = "80 dB",
  `75` = "75 dB",
  `70` = "70 dB",
  `65` = "65 dB",
  `60` = "60 dB",
  `55` = "55 dB",
  `50` = "50 dB",
  `45` = "45 dB",
  `40` = "40 dB",
  `35` = "35 dB",
  `30` = "30 dB",
  `25` = "25 dB",
  `20` = "20 dB",
  `wt` = "",
  `mut`= ""
)
traces.cutoff$`genotype` <- factor(traces.cutoff$`genotype`, levels= c("wt","mut"))

ggplot(traces.cutoff, aes(x = time, y = uv)) +
  stat_summary(aes(group=genotype, color = genotype), fun=mean, geom="line", size = 1.25) +
  stat_summary(aes(group=genotype, fill = genotype), fun.data=mean_se, geom="ribbon", alpha = 0.4) +
  facet_wrap(. ~ `Level(dB)`, nrow = 8, strip.position = "left", labeller = as_labeller(spl.names)) +
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 11),
        axis.title.x = element_text(size = 22),
        axis.title.y = element_text(size = 22, vjust = 2),
        title = element_text(size = 22),
        plot.title = element_text(hjust = 0.5, vjust = 2),
        legend.text = element_text(size = 18),
        legend.key.height= unit(1, 'cm'),
        legend.key.width= unit(1, 'cm'),
        strip.text.y = element_text(size = 12)) +
  labs(title="ABR traces at 16KHz",
       x ="Time (ms)", y = "Level (uv)") +
  scale_x_continuous(breaks= c(0, 2.5, 5, 7.5)) +
  scale_colour_manual(values=cbPalette, labels = c("+/+","del/Y")) +
  scale_fill_manual(values=cbPalette, labels = c("+/+","del/Y"))



#pull out w1 amp and delay

#may have to tweak time window based on where w1 is

traces.peaks <- traces.geno %>% filter(between(time,1,1.7)) %>% group_by(`Sub. ID`, `Level(dB)`) %>%  slice_max(uv)

traces.peaks <- traces.peaks %>%  dplyr::rename(peak = uv)

#check time window

traces.troughs <- traces.geno %>% filter(between(time,1.7,2.5)) %>% group_by(`Sub. ID`, `Level(dB)`) %>%  slice_min(uv)

traces.troughs <- traces.troughs %>% dplyr::rename(trough = uv)
traces.troughs
traces.bind <- bind_cols(traces.peaks, traces.troughs)

traces.peaks.troughs <- traces.bind[,c(1,2,4,5,6,10)]
traces.peaks.troughs

traces.amp <- traces.peaks.troughs %>% mutate(w1amp = as.numeric(peak) - as.numeric(trough))
traces.amp
traces.measure <- traces.amp %>% filter(`Level(dB)...2` %in% c("90","85","80","75","70","65","60","55","50","45")) %>% 
  dplyr::rename(subjectID =`Sub. ID...1`,
        dB = `Level(dB)...2`,
        w1delay = time...5,
        genotype = genotype...6)
?rename
traces.measure$`genotype` <- factor(traces.measure$`genotype`, levels= c("wt","mut"))

traces.measure


ggplot(traces.measure,aes(y = w1delay, x = genotype, color = genotype)) +
  geom_boxplot(size = 1.5) +
  facet_wrap(~ dB) +
  stat_compare_means(method = "kruskal") +
  stat_compare_means(label = "p.signif", method = "kruskal") +
  scale_colour_manual(values=cbPalette, labels = c("+/+","del/Y"))


traces.summary.amp <- summarySE(traces.measure, measurevar="w1amp", groupvars=c("dB", "genotype"))
ggplot(traces.summary.amp, aes(y= w1amp, x=dB, color = genotype)) +
  geom_errorbar(aes(ymin=w1amp-se, ymax=w1amp+se), width=.1, size = 1)+
  geom_point(aes(fill = genotype), size = 5, shape =23)+
  geom_line(aes(x= as.numeric(dB), y = w1amp, color = genotype), size = 1) +
  theme_classic() +
  scale_colour_manual(values=cbPalette, labels = c("+/+","del/Y")) +
  scale_fill_manual(values=cbPalette, labels = c("+/+","del/Y")) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 5000))


traces.summary.delay <- summarySE(traces.measure, measurevar="w1delay", groupvars=c("dB", "genotype"))
ggplot(traces.summary.delay, aes(y= w1delay, x=dB, color = genotype)) +
  geom_errorbar(aes(ymin=w1delay-se, ymax=w1delay+se), width=.1, size = 1)+
  geom_point(aes(fill = genotype), size = 5, shape =23)+
  geom_line(aes(x= as.numeric(dB), y = w1delay, color = genotype), size = 1) +
  theme_classic() +
  scale_colour_manual(values=cbPalette, labels = c("+/+","del/Y")) +
  scale_fill_manual(values=cbPalette, labels = c("+/+","del/Y")) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 2))
  
