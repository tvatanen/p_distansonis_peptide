library(tidyverse)
library(RColorBrewer)

load("diabMGS.RData")
load("mgs_ann_top.RData")

# Search DIABIMMUNE assembled gene clusters for peptide sequence RILVELLYLVCSEYL
# (ran separately in high performance computing environment, code and output below)

#$ grep -B 1 RILVELLYLVCSEYL diabimmune.IGC.catalog.nr.TrustedGenesEnriched.faa
#>G78541_k105_8127_9
#LVLNRILVELLYLVCSEYLNHNSRCLGFYSAYKCITNFGMLLHFFVHRRKSCTFVPSKTERFMDRIRLKDKEFELFIPESDIQAAIAKMAVQIKADVEGKNPLFVGVLNGAFMFVAELMRELDVPYELTFARYSSYQGTSSTGILNEIMPVQADIRGRMVILLEDIIDTGFTMSYVMEKLRSEGAADVRLATMLFKPESLKCELTPDYVGLQIPADFIVGHGLDYDELGRSYKDIYKVVE
#$ grep G78541_k105_8127_9 -n diabimmune.IGC.catalog.nr.TrustedGenesEnriched.GeneNames.txt
#3857448:G78541_k105_8127_9

# The index of the gene identified to contain peptide sequence "RILVELLYLVCSEYL" is 3857448 

idx <- 3857448

hits <- c()
for (j in 1:length(MGSobject$i)) { 
  if (idx %in% MGSobject$i[[j]]) {
    hits <- c(hits, j)
  }
}
MGSobject$attr[hits,]
# this gene is on diabMGS0064 found in 533 samples

mgs_ann_top[ hits , ]
# top annotation for this metagenomic speciesis Parabacteroides distansonis

gene_abundances <- read_tsv("peptide_gene_abundance.txt") %>% gather(gid_wgs, abundance)

# LOAD METADATA
# Metadata is not provided in github due to funding and ethical requirements)
# Metadata is available through website https://diabimmune.broadinstitute.org/diabimmune/

metadata <- 
  full_join(subject, samples) %>%
  full_join(samplesWGS) %>%
  full_join(diabetic_table) %>%
  right_join(gene_abundances) %>%
  filter(!(is.na(subjectID)))


ggplot(metadata, aes(y=abundance, x=age_at_collection, color=country)) +
  geom_point() +
  geom_smooth(method = "loess") +
  scale_y_sqrt() +
  theme_bw()

ggplot(metadata, aes(y=abundance, x=country)) +
  geom_boxplot() +
  theme_bw()

metadata %>%
  mutate(has_sequence = abundance > 0) %>%
  group_by(country) %>%
  summarise(prop_sample = sum(has_sequence) / n()) %>%
  ggplot(aes(y=prop_sample, x=country)) +
  geom_bar(stat = "identity") +
  theme_bw()

ggplot(metadata, aes(y=abundance, x=age_at_collection, color=num_aabs > 1)) +
  geom_point() +
  geom_smooth(method = "loess") +
  scale_y_sqrt() +
  theme_bw()

ggplot(metadata %>% filter(cohort == "t1d"), aes(x=age_at_collection, y=abundance, color = t1d_diagnosed)) +
  geom_point() +
  geom_smooth() +
  scale_y_sqrt()

ggplot(metadata %>% filter(cohort == "t1d"), aes(x=age_at_collection, y=abundance, color = num_aabs > 0)) +
  geom_point() +
  geom_smooth(method = "lm") +
  scale_y_sqrt()

ggplot(metadata, aes(x=age_at_collection, y=abundance, color = num_aabs > 0)) +
  geom_point() +
  geom_smooth(method = "lm") +
  scale_y_sqrt()

metadata_aabs <- 
  metadata %>%
  mutate(aab_pos = num_aabs > 1) %>%
  filter(cohort == "t1d")
           
metadata %>%
  group_by(subjectID, country) %>%
  summarise(hrpt = any(abundance > 0)) %>%
  group_by(country) %>%
  summarise(perc_distansonis = sum(hrpt) / n())
# 77% of Finns, 70% of Estonians and 57% of Russians have P distansonis in at least one sample

metadata %>%
  summarise(perc_distansonis = sum(abundance > 0) / n())
# overall 44.7% of metagenomes have this gene

metadata %>%
  group_by(country) %>%
  summarise(perc_distansonis = sum(abundance > 0) / n())
# 55% of Finnish, 47% of Estonian and 33% of Russian samples

metadata %>%
  filter(country %in% c("FIN", "EST")) %>%
  group_by(subjectID, t1d_diagnosed) %>%
  summarise(hrpt = any(abundance > 0)) %>%
  ungroup() %>%
  group_by(t1d_diagnosed) %>%
  summarise(num_distansonis = sum(hrpt),
            n =  n(),
            prop = num_distansonis / n)
# All 5 T1D progressors have P distansonis peptide
# 71.8% of other subjects have this peptide
prop.test(c(140,5), c(195,5))
# Difference is not statistically significant (p = 0.37)

metadata %>%
  filter(country %in% c("FIN", "EST")) %>%
  mutate(seroconverted = num_aabs > 1) %>%
  group_by(subjectID, seroconverted) %>%
  summarise(hrpt = any(abundance > 0)) %>%
  ungroup() %>%
  group_by(seroconverted) %>%
  summarise(num_distansonis = sum(hrpt),
            n = n(),
            prop = num_distansonis / n)
# 92.9% of seroconverters (>1 AABs) have P distansonis peptide
# 71.0% of non-seroconverters have P distansonis peptide
prop.test(c(132,13), c(186,14))
# Difference is not statistically significant (p = 0.14)

# Include Russians
metadata %>%
  mutate(seroconverted = num_aabs > 1) %>%
  group_by(subjectID, seroconverted) %>%
  summarise(hrpt = any(abundance > 0)) %>%
  ungroup() %>%
  group_by(seroconverted) %>%
  summarise(num_distansonis = sum(hrpt),
            n = n(),
            prop = num_distansonis / n)
# 93.3% of seroconverters (>1 AABs) have P distansonis peptide
# 67.3% of non-seroconverters have P distansonis peptide
prop.test(c(171,14), c(254,15))
# Difference is not statistically significant (p = 0.068)

metadata_binomial_model <-
  metadata %>%
  mutate(one_or_more_aab = num_aabs > 0,
         two_or_more_aabs = num_aabs > 1) %>%
  mutate(hrpt = abundance > 0) %>%
  mutate(year = case_when(age_at_collection < 365 ~ "year1",
                          age_at_collection < 365*2 ~"year2",
                          TRUE ~ "year3")) %>%
  mutate(age_norm = (age_at_collection - mean(age_at_collection)) / sd(age_at_collection))

# share metadata
library(openxlsx)
metadata_binomial_model %>% 
  dplyr::select(subjectID, country, gender, sampleID, gid_wgs,
                age_at_collection, cohort,
                IAA, GADA, IA2A, ZNT8A, ICA, num_aabs, one_or_more_aab,
                two_or_more_aabs, hrpt, year) %>%
  write.xlsx("diabimmune_hrpt_aab_t1d_data.xlsx")

metadata_binomial_model %>%
  filter(!(is.na(one_or_more_aab))) %>%
  distinct(subjectID, one_or_more_aab) %>%
  count(one_or_more_aab)
# 40 individuals have one ore more AAB

metadata_binomial_model %>%
  filter(!(is.na(two_or_more_aabs))) %>%
  distinct(subjectID, two_or_more_aabs) %>%
  count(two_or_more_aabs)
# 15 individuals have one ore more AAB

metadata_binomial_model %>%
  filter(!(is.na(IAA))) %>%
  distinct(subjectID, IAA) %>%
  count(IAA)
# 21 individuals have insulin AABs

metadata_binomial_model %>%
  filter(!(is.na(t1d_diagnosed))) %>%
  distinct(subjectID, t1d_diagnosed) %>%
  count(t1d_diagnosed)
# 6 individuals have T1D

metadata_binomial_model %>%
  filter(!(is.na(one_or_more_aab)),
         !(is.na(t1d_diagnosed))) %>%
  nrow()
# N = 1164 samples

library(lme4)
m_aab <- glmer(hrpt ~ age_norm + country + one_or_more_aab +
             (1 | subjectID),
           data = metadata_binomial_model,
           family = binomial,
           control = glmerControl(optimizer = "bobyqa"),
           nAGQ = 10)
summary(m_aab)
# One or more AABs is not associated with hrpt presence (p=0.86)

m_aabs <- glmer(hrpt ~ age_norm + country + two_or_more_aabs +
                  (1 | subjectID),
                data = metadata_binomial_model,
                family = binomial,
                control = glmerControl(optimizer = "bobyqa"),
                nAGQ = 10)
summary(m_aabs)
# Two ore more AABs is associated with hrpt presence (p=0.0017)

m_iaa <- glmer(hrpt ~ age_norm + country + IAA +
                  (1 | subjectID),
                data = metadata_binomial_model,
                family = binomial,
                control = glmerControl(optimizer = "bobyqa"),
                nAGQ = 10)
summary(m_iaa)
# IAA is weakly associated with hrpt presence (p=0.047)

m_t1d <- glmer(hrpt ~ age_norm + country + t1d_diagnosed +
                 (1 | subjectID),
               data = metadata_binomial_model %>%
                 filter(!(is.na(t1d_diagnosed))),
               family = binomial,
               control = glmerControl(optimizer = "bobyqa"),
               nAGQ = 10)
summary(m_t1d)
# T1D status is weakly associated with hrpt presence (p=0.020)

metadata_binomial_model %>%
  filter(!(is.na(one_or_more_aab))) %>%
  group_by(year, one_or_more_aab, country) %>%
  summarise(prop_peptide = sum(hrpt) / n()) %>%
  mutate(country = case_when(country == "FIN" ~ "Finland",
                             country == "EST" ~ "Estonia",
                             country == "RUS" ~ "Russia"),
         country = factor(country,
                          levels = c("Finland", "Estonia", "Russia"))) %>%
  ggplot(aes(x=year, fill=one_or_more_aab, y=prop_peptide)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_bw() +
  xlab("") +
  facet_grid(~country) +
  scale_fill_manual(values = brewer.pal(5, "Set2")[c(5,3)]) +
  ylab("Proportion of samples with hprt:4-18 peptide")
ggsave("one_or_more_aab_vs_peptide_040622.pdf", width = 6, height = 4)

metadata_binomial_model %>%
  filter(!(is.na(one_or_more_aab))) %>%
  group_by(year, one_or_more_aab, country) %>%
  summarise(n_hrpt = sum(hrpt),
            n = n()) %>%
  write.xlsx("one_or_more_aab_vs_peptide_040622.xlsx")

metadata_binomial_model %>%
  filter(!(is.na(two_or_more_aabs))) %>%
  group_by(year, two_or_more_aabs, country) %>%
  summarise(prop_peptide = sum(hrpt) / n()) %>%
  mutate(country = case_when(country == "FIN" ~ "Finland",
                             country == "EST" ~ "Estonia",
                             country == "RUS" ~ "Russia"),
         country = factor(country,
                          levels = c("Finland", "Estonia", "Russia"))) %>%
  ggplot(aes(x=year, fill=two_or_more_aabs, y=prop_peptide)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_bw() +
  xlab("") +
  facet_grid(~country) +
  scale_fill_manual(values = brewer.pal(5, "Set2")[c(5,3)]) +
  ylab("Proportion of samples with hprt:4-18 peptide")
ggsave("two_or_more_aabs_vs_peptide_040622.pdf", width = 6, height = 4)

metadata_binomial_model %>%
  filter(!(is.na(two_or_more_aabs))) %>%
  group_by(year, two_or_more_aabs, country) %>%
  summarise(n_hrpt = sum(hrpt),
            n = n()) %>%
  write.xlsx("two_or_more_aabs_vs_peptide_040622.xlsx")


metadata_binomial_model %>%
  filter(!(is.na(IAA))) %>%
  group_by(year, IAA, country) %>%
  summarise(prop_peptide = sum(hrpt) / n()) %>%
  mutate(country = case_when(country == "FIN" ~ "Finland",
                             country == "EST" ~ "Estonia",
                             country == "RUS" ~ "Russia"),
         country = factor(country,
                          levels = c("Finland", "Estonia", "Russia"))) %>%
  ggplot(aes(x=year, fill=IAA, y=prop_peptide)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_bw() +
  xlab("") +
  facet_grid(~country) +
  scale_fill_manual(values = brewer.pal(5, "Set2")[c(5,3)]) +
  ylab("Proportion of samples with hprt:4-18 peptide")
ggsave("iaa_vs_peptide_040622.pdf", width = 6, height = 4)

metadata_binomial_model %>%
  filter(!(is.na(IAA))) %>%
  group_by(year, IAA, country) %>%
  summarise(n_hrpt = sum(hrpt),
            n = n()) %>%
  write.xlsx("iaa_vs_peptide_040622.xlsx")

metadata_binomial_model %>%
  filter(!(is.na(t1d_diagnosed))) %>%
  group_by(year, t1d_diagnosed, country) %>%
  summarise(prop_peptide = sum(hrpt) / n()) %>%
  mutate(country = case_when(country == "FIN" ~ "Finland",
                             country == "EST" ~ "Estonia",
                             country == "RUS" ~ "Russia"),
         country = factor(country,
                          levels = c("Finland", "Estonia", "Russia"))) %>%
  ggplot(aes(x=year, fill=t1d_diagnosed, y=prop_peptide)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_bw() +
  xlab("") +
  facet_grid(~country) +
  scale_fill_manual(values = brewer.pal(3, "Set2")) +
  ylab("Proportion of samples with hprt4-18 like peptide")
ggsave("t1d_vs_peptide_040622.pdf", width = 6, height = 4)

metadata_binomial_model %>%
  filter(!(is.na(t1d_diagnosed))) %>%
  group_by(year, t1d_diagnosed, country) %>%
  summarise(n_hrpt = sum(hrpt),
            n = n()) %>%
  write.xlsx("t1d_vs_peptide_040622.xlsx")

# same plots per subject, not sample
metadata_binomial_model %>%
  filter(!(is.na(one_or_more_aab))) %>%
  group_by(year, subjectID, one_or_more_aab, country) %>%
  summarise(has_hrpt = any(hrpt)) %>%
  ungroup() %>%
  group_by(year, one_or_more_aab, country) %>%
  summarise(prop_peptide = sum(has_hrpt) / n()) %>%
  mutate(country = case_when(country == "FIN" ~ "Finland",
                             country == "EST" ~ "Estonia",
                             country == "RUS" ~ "Russia"),
         country = factor(country,
                          levels = c("Finland", "Estonia", "Russia"))) %>%
  ggplot(aes(x=year, fill=one_or_more_aab, y=prop_peptide)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_bw() +
  xlab("") +
  facet_grid(~country) +
  scale_fill_manual(values = brewer.pal(5, "Set2")[c(5,3)]) +
  ylab("Proportion of subjects with hprt:4-18 peptide")
ggsave("one_or_more_aab_vs_peptide_per_subject_040622.pdf", width = 6, height = 4)

metadata_binomial_model %>%
  filter(!(is.na(one_or_more_aab))) %>%
  group_by(year, subjectID, one_or_more_aab, country) %>%
  summarise(has_hrpt = any(hrpt)) %>%
  ungroup() %>%
  group_by(year, one_or_more_aab, country) %>%
  summarise(n_hrpt = sum(has_hrpt),
            n = n()) %>%
  write.xlsx("one_or_more_aab_vs_peptide_per_subject_040622.xlsx")

metadata_binomial_model %>%
  filter(!(is.na(two_or_more_aabs))) %>%
  group_by(year, subjectID, two_or_more_aabs, country) %>%
  summarise(has_hrpt = any(hrpt)) %>%
  ungroup() %>%
  group_by(year, two_or_more_aabs, country) %>%
  summarise(prop_peptide = sum(has_hrpt) / n()) %>%
  mutate(country = case_when(country == "FIN" ~ "Finland",
                             country == "EST" ~ "Estonia",
                             country == "RUS" ~ "Russia"),
         country = factor(country,
                          levels = c("Finland", "Estonia", "Russia"))) %>%
  ggplot(aes(x=year, fill=two_or_more_aabs, y=prop_peptide)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_bw() +
  xlab("") +
  facet_grid(~country) +
  scale_fill_manual(values = brewer.pal(5, "Set2")[c(5,3)]) +
  ylab("Proportion of subjects with hprt:4-18 peptide")
ggsave("two_or_more_aab_vs_peptide_per_subject_040622.pdf", width = 6, height = 4)

metadata_binomial_model %>%
  filter(!(is.na(two_or_more_aabs))) %>%
  group_by(year, subjectID, two_or_more_aabs, country) %>%
  summarise(has_hrpt = any(hrpt)) %>%
  ungroup() %>%
  group_by(year, two_or_more_aabs, country) %>%
  summarise(n_hrpt = sum(has_hrpt),
            n = n()) %>%
  write.xlsx("two_or_more_aab_vs_peptide_per_subject_040622.xlsx")

metadata_binomial_model %>%
  filter(!(is.na(IAA))) %>%
  group_by(year, subjectID, IAA, country) %>%
  summarise(has_hrpt = any(hrpt)) %>%
  ungroup() %>%
  group_by(year, IAA, country) %>%
  summarise(prop_peptide = sum(has_hrpt) / n()) %>%
  mutate(country = case_when(country == "FIN" ~ "Finland",
                             country == "EST" ~ "Estonia",
                             country == "RUS" ~ "Russia"),
         country = factor(country,
                          levels = c("Finland", "Estonia", "Russia"))) %>%
  ggplot(aes(x=year, fill=IAA, y=prop_peptide)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_bw() +
  xlab("") +
  facet_grid(~country) +
  scale_fill_manual(values = brewer.pal(5, "Set2")[c(5,3)]) +
  ylab("Proportion of subjects with hprt:4-18 peptide")
ggsave("iaa_vs_peptide_per_subject_040622.pdf", width = 6, height = 4)

metadata_binomial_model %>%
  filter(!(is.na(IAA))) %>%
  group_by(year, subjectID, IAA, country) %>%
  summarise(has_hrpt = any(hrpt)) %>%
  ungroup() %>%
  group_by(year, IAA, country) %>%
  summarise(n_hrpt = sum(has_hrpt),
            n = n()) %>%
  write.xlsx("iaa_vs_peptide_per_subject_040622.xlsx")


metadata_binomial_model %>%
  filter(!(is.na(t1d_diagnosed))) %>%
  group_by(year, subjectID, t1d_diagnosed, country) %>%
  summarise(has_hrpt = any(hrpt)) %>%
  ungroup() %>%
  group_by(year, t1d_diagnosed, country) %>%
  summarise(prop_peptide = sum(has_hrpt) / n()) %>%
  mutate(country = case_when(country == "FIN" ~ "Finland",
                             country == "EST" ~ "Estonia",
                             country == "RUS" ~ "Russia"),
         country = factor(country,
                          levels = c("Finland", "Estonia", "Russia"))) %>%
  ggplot(aes(x=year, fill=t1d_diagnosed, y=prop_peptide)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_bw() +
  xlab("") +
  facet_grid(~country) +
  scale_fill_manual(values = brewer.pal(3, "Set2")) +
  ylab("Proportion of subjects with hprt:4-18 peptide")
ggsave("t1d_vs_peptide_per_subject_040622.pdf", width = 6, height = 4)

metadata_binomial_model %>%
  filter(!(is.na(t1d_diagnosed))) %>%
  group_by(year, subjectID, t1d_diagnosed, country) %>%
  summarise(has_hrpt = any(hrpt)) %>%
  ungroup() %>%
  group_by(year, t1d_diagnosed, country) %>%
  summarise(n_hrpt = sum(has_hrpt),
            n = n()) %>%
  write.xlsx("t1d_vs_peptide_per_subject_040622.xlsx")

# plot without countries separated
metadata_binomial_model %>%
  filter(!(is.na(one_or_more_aab))) %>%
  group_by(year, one_or_more_aab) %>%
  summarise(prop_peptide = sum(hrpt) / n()) %>%
  ggplot(aes(x=year, fill=one_or_more_aab, y=prop_peptide)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_bw() +
  xlab("") +
  scale_fill_manual(values = brewer.pal(5, "Set2")[c(5,3)]) +
  ylab("Proportion of samples with hprt:4-18 peptide")
ggsave("one_or_more_aab_vs_peptide_pooled_040622.pdf", width = 6, height = 4)

metadata_binomial_model %>%
  filter(!(is.na(one_or_more_aab))) %>%
  group_by(year, one_or_more_aab) %>%
  summarise(n_hrpt = sum(hrpt),
            n = n()) %>%
  write.xlsx("one_or_more_aab_vs_peptide_pooled_040622.xlsx")

metadata_binomial_model %>%
  filter(!(is.na(two_or_more_aabs))) %>%
  group_by(year, two_or_more_aabs) %>%
  summarise(prop_peptide = sum(hrpt) / n()) %>%
  ggplot(aes(x=year, fill=two_or_more_aabs, y=prop_peptide)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_bw() +
  xlab("") +
  scale_fill_manual(values = brewer.pal(5, "Set2")[c(5,3)]) +
  ylab("Proportion of samples with hprt:4-18 peptide")
ggsave("two_or_more_aabs_vs_peptide_pooled_040622.pdf", width = 6, height = 4)

metadata_binomial_model %>%
  filter(!(is.na(two_or_more_aabs))) %>%
  group_by(year, two_or_more_aabs) %>%
  summarise(n_hrpt = sum(hrpt),
            n = n()) %>%
  write.xlsx("two_or_more_aabs_vs_peptide_pooled_040622.xlsx")

metadata_binomial_model %>%
  filter(!(is.na(IAA))) %>%
  group_by(year, IAA) %>%
  summarise(prop_peptide = sum(hrpt) / n()) %>%
  ggplot(aes(x=year, fill=IAA, y=prop_peptide)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_bw() +
  xlab("") +
  scale_fill_manual(values = brewer.pal(5, "Set2")[c(5,3)]) +
  ylab("Proportion of samples with hprt:4-18 peptide")
ggsave("iaa_vs_peptide_pooled_040622.pdf", width = 6, height = 4)

metadata_binomial_model %>%
  filter(!(is.na(IAA))) %>%
  group_by(year, IAA) %>%
  summarise(n_hrpt = sum(hrpt),
            n = n()) %>%
  write.xlsx("iaa_vs_peptide_pooled_040622.xlsx")

metadata_binomial_model %>%
  filter(!(is.na(t1d_diagnosed))) %>%
  group_by(year, t1d_diagnosed) %>%
  summarise(prop_peptide = sum(hrpt) / n()) %>%
  ggplot(aes(x=year, fill=t1d_diagnosed, y=prop_peptide)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_bw() +
  xlab("") +
  scale_fill_manual(values = brewer.pal(3, "Set2")) +
  ylab("Proportion of samples with hprt4-18 like peptide")
ggsave("t1d_vs_peptide_pooled_040622.pdf", width = 6, height = 4)

metadata_binomial_model %>%
  filter(!(is.na(t1d_diagnosed))) %>%
  group_by(year, t1d_diagnosed) %>%
  summarise(n_hrpt = sum(hrpt),
            n = n()) %>%
  write.xlsx("t1d_vs_peptide_pooled_040622.xlsx")

# same plots per subject, not sample
metadata_binomial_model %>%
  filter(!(is.na(one_or_more_aab))) %>%
  group_by(year, subjectID, one_or_more_aab) %>%
  summarise(has_hrpt = any(hrpt)) %>%
  ungroup() %>%
  group_by(year, one_or_more_aab) %>%
  summarise(prop_peptide = sum(has_hrpt) / n()) %>%
  ggplot(aes(x=year, fill=one_or_more_aab, y=prop_peptide)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_bw() +
  xlab("") +
  scale_fill_manual(values = brewer.pal(5, "Set2")[c(5,3)]) +
  ylab("Proportion of subjects with hprt:4-18 peptide")
ggsave("one_or_more_aab_vs_peptide_per_subject_pooled_040622.pdf", width = 6, height = 4)

metadata_binomial_model %>%
  filter(!(is.na(one_or_more_aab))) %>%
  group_by(year, subjectID, one_or_more_aab) %>%
  summarise(has_hrpt = any(hrpt)) %>%
  ungroup() %>%
  group_by(year, one_or_more_aab) %>%
  summarise(n_hrpt = sum(has_hrpt),
            n = n()) %>%
  write.xlsx("one_or_more_aab_vs_peptide_per_subject_pooled_040622.xlsx")

metadata_binomial_model %>%
  filter(!(is.na(two_or_more_aabs))) %>%
  group_by(year, subjectID, two_or_more_aabs) %>%
  summarise(has_hrpt = any(hrpt)) %>%
  ungroup() %>%
  group_by(year, two_or_more_aabs) %>%
  summarise(prop_peptide = sum(has_hrpt) / n()) %>%
  ggplot(aes(x=year, fill=two_or_more_aabs, y=prop_peptide)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_bw() +
  xlab("") +
  scale_fill_manual(values = brewer.pal(5, "Set2")[c(5,3)]) +
  ylab("Proportion of subjects with hprt:4-18 peptide")
ggsave("two_or_more_aab_vs_peptide_per_subject_pooled_040622.pdf", width = 6, height = 4)

metadata_binomial_model %>%
  filter(!(is.na(two_or_more_aabs))) %>%
  group_by(year, subjectID, two_or_more_aabs) %>%
  summarise(has_hrpt = any(hrpt)) %>%
  ungroup() %>%
  group_by(year, two_or_more_aabs) %>%
  summarise(n_hrpt = sum(has_hrpt),
            n = n()) %>%
  write.xlsx("two_or_more_aab_vs_peptide_per_subject_pooled_040622.xlsx")

metadata_binomial_model %>%
  filter(!(is.na(IAA))) %>%
  group_by(year, subjectID, IAA) %>%
  summarise(has_hrpt = any(hrpt)) %>%
  ungroup() %>%
  group_by(year, IAA) %>%
  summarise(prop_peptide = sum(has_hrpt) / n()) %>%
  ggplot(aes(x=year, fill=IAA, y=prop_peptide)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_bw() +
  xlab("") +
  scale_fill_manual(values = brewer.pal(5, "Set2")[c(5,3)]) +
  ylab("Proportion of subjects with hprt:4-18 peptide")
ggsave("iaa_vs_peptide_per_subject_pooled_040622.pdf", width = 6, height = 4)

metadata_binomial_model %>%
  filter(!(is.na(IAA))) %>%
  group_by(year, subjectID, IAA) %>%
  summarise(has_hrpt = any(hrpt)) %>%
  ungroup() %>%
  group_by(year, IAA) %>%
  summarise(n_hrpt = sum(has_hrpt),
            n = n()) %>%
  write.xlsx("iaa_vs_peptide_per_subject_pooled_040622.xlsx")

metadata_binomial_model %>%
  filter(!(is.na(t1d_diagnosed))) %>%
  group_by(year, subjectID, t1d_diagnosed) %>%
  summarise(has_hrpt = any(hrpt)) %>%
  ungroup() %>%
  group_by(year, t1d_diagnosed) %>%
  summarise(prop_peptide = sum(has_hrpt) / n()) %>%
  ggplot(aes(x=year, fill=t1d_diagnosed, y=prop_peptide)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_bw() +
  xlab("") +
  scale_fill_manual(values = brewer.pal(3, "Set2")) +
  ylab("Proportion of subjects with hprt:4-18 peptide")
ggsave("t1d_vs_peptide_per_subject_pooled_040622.pdf", width = 6, height = 4)

metadata_binomial_model %>%
  filter(!(is.na(t1d_diagnosed))) %>%
  group_by(year, subjectID, t1d_diagnosed) %>%
  summarise(has_hrpt = any(hrpt)) %>%
  ungroup() %>%
  group_by(year, t1d_diagnosed) %>%
  summarise(n_hrpt = sum(has_hrpt),
            n = n()) %>%
  write.xlsx("t1d_vs_peptide_per_subject_pooled_040622.xlsx")

# Seroconversion rate calculations

seroconversion_rate_data <-
  metadata_binomial_model %>%
  filter(year == "year1",
         !(is.na(two_or_more_aabs))) %>%
  group_by(subjectID, two_or_more_aabs) %>%
  summarise(has_peptide = any(hrpt),
            window = "year1") 

seroconversion_rate_data <- 
  metadata_binomial_model %>%
  filter(year %in% c("year1" ,"year2"),
         !(is.na(two_or_more_aabs))) %>%
  group_by(subjectID, two_or_more_aabs) %>%
  summarise(has_peptide = any(hrpt),
            window = "year1-year2") %>%
  bind_rows(seroconversion_rate_data)

seroconversion_rate_data <- 
  metadata_binomial_model %>%
  filter(!(is.na(two_or_more_aabs))) %>%
  group_by(subjectID, two_or_more_aabs) %>%
  summarise(has_peptide = any(hrpt),
            window = "year1-year3") %>%
  bind_rows(seroconversion_rate_data)

seroconversion_rate_data %>%
  group_by(window, two_or_more_aabs, has_peptide) %>%
  summarise(n = n()) %>%
  pivot_wider(names_from = "two_or_more_aabs", values_from = "n") %>%
  rename(two_or_more_aabs = `TRUE`,
         no_two_aabs = `FALSE`) %>%
  mutate(n = two_or_more_aabs + no_two_aabs,
         seroconversion_rate = two_or_more_aabs / n) %>%
  arrange(has_peptide)

prop.test(x = c(7, 7), n = c(84, 129))
# year 1: p = 0.58

prop.test(x = c(1, 14), n = c(97, 170))
# year 1-2: p = 0.029

prop.test(x = c(1, 14), n = c(84, 185))
# year 1-3: p = 0.068
