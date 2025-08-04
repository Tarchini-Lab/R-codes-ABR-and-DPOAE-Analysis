## DPOAE data pipeline
# Elli Hartig 09/23/22

# Load libraries
library(tidyverse)

# Import data from BioSig exports. This script pulls all csv files from a folder and compiles them
setwd("H:/")
getwd()
data.raw <- list.files(path = "./AGG_dpoae", full.names = TRUE) %>% 
  lapply(read_csv) %>% 
  bind_rows()

# Average measurements to determine background noise
noisefloor <- rowMeans(data.raw[49:2096])

# Select only columns needed
data.trim <- data.raw[, c(8, 13, 14, 29)]

# Add noise floor data
data.addnf <- cbind(data.trim, `noisefloor(dBv)` = noisefloor)

# Create new columns for noise floor and dp values converted to dB instead of dBv
data.dB <- data.addnf %>% 
  mutate(`noisefloor(dB)` = 20 * log10((10^(`noisefloor(dBv)` / 20)) / 0.05) + 93.9) %>% 
  mutate(`dp(dB)` = 20 * log10((10^(`V3(dBv)` / 20)) / 0.05) + 93.9)

# Create new column for condition
wildtypes <- c("AGG36", "AGG37", "AGG39", "AGG41", "AGG50", "AGG52", "AGG54", "AGG56", "AGG57", "AGG59", "AGG60", "AGG74", "AGG77", "AGG78", "AGG80","AGG84", "AGG85", "AGG33", "AGG34", "AGG35", "AGG38", "AGG40", "AGG44", "AGG46", "AGG49", "AGG73", "AGG76", "AGG79", "AGG82", "AGG83", "AGG86")
hets <- c("AGG43", "AGG45", "AGG48", "AGG53", "AGG55", "AGG58", "AGG75", "AGG81")

data.dB.less <- data.dB[-c(1:175), ]

data.geno <- data.dB.less %>% 
  mutate(genotype = case_when(`Sub. ID` %in% wildtypes ~ "control",
                              `Sub. ID` %in% hets ~ "het",
                              .default = "mut"))
data.geno <- data.geno %>% na.omit()

# Order dataset for aesthetics
data.geno$`AudFreq()` <- factor(data.geno$`AudFreq()`, levels = c("8000", "12000", "16000", "24000", "32000"))
data.geno$genotype <- factor(data.geno$genotype, levels = c("control", "het", "mut"))

# Calculate the noise floor for each frequency and level
noise_floor_by_freq_level <- data.geno %>%
  group_by(`AudFreq()`, `Level()`) %>%
  summarise(noisefloor_dB = mean(`noisefloor(dB)`))

# Merge the noise floor data back into the main dataset
data.geno <- data.geno %>%
  left_join(noise_floor_by_freq_level, by = c("AudFreq()", "Level()"))

# Now, make sure to calculate the noise floor for each level of `Level()`
ggplot(data.geno, aes(x = `Level()`, y = `dp(dB)`)) +
  stat_summary(aes(group = genotype, color = genotype), fun = mean, geom = "point", size = 4) +
  stat_summary(aes(group = genotype, color = genotype), fun = mean, geom = "line", size = 1) +
  stat_summary(aes(group = genotype, color = genotype), fun.data = "mean_se", geom = "errorbar", width = 0.2, size = 1) +
  geom_line(aes(y = noisefloor_dB, group = `AudFreq()`), color = "gray", alpha = 0.6, size = 1) +
  facet_wrap(. ~ `AudFreq()`)

# Create limited dataset for selected frequencies
dpoae.limited <- data.geno %>% filter(`AudFreq()` %in% c("8000", "12000", "16000", "24000", "32000"))

dpoae.limited$`AudFreq()` <- factor(dpoae.limited$`AudFreq()`, 
                                    levels = c("8000", "12000", "16000", "24000", "32000"), 
                                    labels = c("8KHz", "12 KHz", "16 KHz", "24 KHz", "32 KHz"))

dpoae.limited$`Level()` <- factor(dpoae.limited$`Level()`)

# Set color palette
cbPalette <- c("control" = "black", "het" = "red", "mut" = "red")

# New facet label names for dose variable
ggplot(dpoae.limited, aes(x = `Level()`, y = `dp(dB)`)) +
  stat_summary(aes(group = genotype, color = genotype), fun = mean, geom = "point", size = 4) +
  stat_summary(aes(group = genotype, color = genotype), fun = mean, geom = "line", size = 1) +
  stat_summary(aes(group = genotype, color = genotype), fun.data = "mean_se", geom = "errorbar", width = 0.2, size = 1) +
  geom_line(aes(y = noisefloor_dB, group = `AudFreq()`), color = "gray", alpha = 0.6, size = 1) +
  facet_wrap(. ~ `AudFreq()`) +
  scale_colour_manual(values = cbPalette) +
  xlab("Sound Pressure Level (dB)") +
  ylab("Distortion Product (dB)") +
  theme(strip.text.x = element_text(size = 20),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        legend.title = element_blank(),
        legend.text = element_text(size = 18))



###STATS
install.packages("rstatix")
library(rstatix)
# Backticks in the names will mess up the Tukey test, so will having variables listed as numeric
dpoae.stat <- dplyr::rename(dpoae.limited, dp = `dp(dB)`, 
                            freq = `AudFreq()`, 
                            level = `Level()`)
three.way.aov <- aov(data = dpoae.stat, dp ~ genotype * freq * level)
summary(three.way.aov)

tukey <- tukey_hsd(three.way.aov)
print(tukey)

