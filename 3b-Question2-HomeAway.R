#### Legume Meta Analysis ####
#Authors: Smriti Pehim Limbu, Kim Komatsu, Kathryn Bloodworth

library(car)
library(lme4)
library(lmerTest)
library(ggsci)
library(ggrepel)
library(tidyverse)
library(dplyr)
library(tidyr)
library(ggplot2)
library(glmmTMB)

###bar graph summary statistics function
#barGraphStats(data=, variable="", byFactorNames=c(""))
barGraphStats <- function(data, variable, byFactorNames) {
  count <- length(byFactorNames)
  N <- aggregate(data[[variable]], data[byFactorNames], FUN=length)
  names(N)[1:count] <- byFactorNames
  names(N) <- sub("^x$", "N", names(N))
  mean <- aggregate(data[[variable]], data[byFactorNames], FUN=mean)
  names(mean)[1:count] <- byFactorNames
  names(mean) <- sub("^x$", "mean", names(mean))
  sd <- aggregate(data[[variable]], data[byFactorNames], FUN=sd)
  names(sd)[1:count] <- byFactorNames
  names(sd) <- sub("^x$", "sd", names(sd))
  preSummaryStats <- merge(N, mean, by=byFactorNames)
  finalSummaryStats <- merge(preSummaryStats, sd, by=byFactorNames)
  finalSummaryStats$se <- finalSummaryStats$sd / sqrt(finalSummaryStats$N)
  return(finalSummaryStats)
}  


##### data import #####

# getting clean legume names for standardizing across papers
cleanNames <- read.csv('legume_clean_names.csv') %>% 
  separate(old_genus_species,c("genus","species","subspecies","extra"),sep=" ") %>% 
  unite(col=genus_species, c(genus,species,subspecies,extra), sep='_',na.rm=TRUE) %>% 
  rename(clean_name=new_name) %>% 
  select(-notes)

# getting native/invasive status for each study and legume species
plantStatus <- read.csv('legume_strain diversity_meta analysis_plant associations_edited names_presence absence.csv') %>% 
  left_join(cleanNames) %>%
  rename(paper_id=paper_id) %>% 
  select(paper_id, clean_name, plant_status, sample_country) %>% 
  unique()

#correlation test when nodules number is 3
plantData_n3plus <- plantData %>%
  filter(!is.na(num_nodules), !is.na(strain_richness), num_nodules >= 2)

cor_res <- cor.test(
  plantData_n3plus$num_nodules,
  plantData_n3plus$strain_richness,
  method = "spearman"
)

cor_res

# removing species where legume species isn't clear, native/invasive status is unknown, and cultivated species
plantData <- read.csv("legume_strain diversity_meta analysis_plant data.csv") %>% 
  #change all 999 to NAs
  mutate(across(all_of(c('sites_sampled', 'num_nodules', 'strain_richness', 'proportion_novel_strains',
                         'proportion_familiar_strains', 'proportion_overlapping_strains', 'num_spp_compared')), ~na_if(., 999))) %>% 
  left_join(cleanNames) %>% 
  select(-notes, -genus_species) %>% 
  full_join(plantStatus) %>% 
  mutate(plant_status=ifelse(plant_status=='invasive', 'introduced', plant_status)) %>% 
  filter(cultivation.status!='row crop',
         plant_status!="999",
         !is.na(plant_status),
         species!="sp",
         species!="spp",
         species!="sp.",
         species!="",
         paper_id!=142,
         plant_status!="extinct") %>% 
  # filter(!(sample_country %in% c("Japan-China","Japan_China","Kenya-Sudan","Czech_Republic-France-Georgia-Hungary-Italy-Romania-Spain","Senegal-Mauritania-Tunisia-Burundi","Malawi-Zambia-Kenya","Brazil-Venezuela"))) %>% 
  filter(!is.na(strain_richness)) %>% 
  mutate(strain_richness=as.numeric(strain_richness)) %>% #NA introduced for blanks
  mutate(genetic_region=(ifelse(genetic_region %in% c('16S rRNA', 'PCR-RFLP', 'RFLP-16S'), '16S',
                         ifelse(genetic_region=='IGS-RFLP', 'IGS',
                         genetic_region))))

plantData[plantData==""]<-NA

#### Supplemental Figure S1 ####
#Relating number of nodules collected to strain richness. Below a cutoff of 3, the number of strains is nearly always equal to the number of nodules. At 3 and up, this is less of a problem. Subsetting data to 3+ nodules.
ggplot(data=plantData, aes(x=num_nodules, y=strain_richness)) +
  geom_jitter() +
  geom_abline(slope=1) +
  coord_cartesian(xlim=c(1,10), ylim=c(0,10)) +
  scale_x_continuous(breaks=seq(0,10,2)) +
  scale_y_continuous(breaks=seq(0,10,2)) +
  xlab('Number of Nodules Sampled') + ylab('Strain Richness')


#### Question 2 data and model ####

paper_plant_status <- plantData %>%
  filter(num_nodules>2) %>%
  filter(!is.na(clean_name),
         compares_homeaway==1) %>% 
  group_by(clean_name) %>% 
  filter(all(c("native", "introduced") %in% plant_status)) %>% 
  select(paper_id, genus, species, plant_status, strain_richness) %>% 
  group_by(paper_id, genus, species, plant_status) %>% 
  summarize(strain_richness_mean=mean(strain_richness),
            strain_richness_n=length(strain_richness),
            strain_richness_sd=sd(strain_richness)) %>% 
  ungroup() %>% 
  pivot_wider(names_from=plant_status, values_from=c('strain_richness_mean', 'strain_richness_n', 'strain_richness_sd')) %>% 
  mutate(strain_richness_lnRR=(log(strain_richness_mean_introduced/strain_richness_mean_native)))

str(paper_plant_status)

hist(paper_plant_status$strain_richness_lnRR)

#overall effect size
t.test(paper_plant_status$strain_richness_lnRR, y = NULL,
       alternative = "two.sided",
       mu = 0, paired = FALSE, var.equal = FALSE,
       conf.level = 0.95)



######################## (HOME VS. AWAY) STRAIN RICHNESS DIFFERENCE###############
allSpp_3ab <- plantData %>% 
  filter(num_nodules>2) %>%
  select(clean_name, plant_status, compares_homeaway, compares_natinv, genetic_region) %>% 
  select(clean_name, plant_status, genetic_region) %>% 
  unique() %>%
  mutate(yes=1) %>% 
  pivot_wider(names_from=plant_status, values_from=yes) %>% #figure out which species have the same genetic region from native and introduced ranges
  mutate(both_ranges=native+introduced) %>% 
  filter(both_ranges==2) %>%  
  select(-native, -introduced)
colnames(plantData)

nativeInvasiveGenetic <- plantData %>% 
  right_join(allSpp_3ab) %>% 
  filter(num_nodules>2) %>% 
  select(paper_id, clean_name, plant_status, annual_perennial, growth_form, habitat_type, sample_country,sample_continent, num_nodules, num_plants, genetic_region, strain_richness, cultivation.status) %>%
  mutate(genetic_region2=ifelse(genetic_region %in% c('16S', '16S-23S', '16S-ARDRA', '16S-IGS','16S-RFLP', '16S rDNA', '16S_23S_RFLP', '16S_ARDRA', '16S_BLAST', '16S_PCR_RFLP', '16S_rDNA', '16S_RFLP', 'ARDRA', 'PCR','PCR_RFLP', '16S_PCR-RFLP', 'RFLP_16S','RFLP', 'RFLP_PCR'), '16S',ifelse(genetic_region %in% c('23S', '23S_IVS', 'RFLP-23S'), '23S',ifelse(genetic_region %in% c('AFLP', 'AFLP_Pst-A', 'AFLP_Pst-G', 'AFLP_Pst-GC'), 'AFLP',ifelse(genetic_region %in% c('BOX', 'BOX-AIR', 'box-PCR', 'BOX-PCR','Box_A1R-PCR', 'BOX_PCR', 'BoxA1R','BOXA1R', 'BOXA1R-PCR', 'BOXAIR'), 'BOX',ifelse(genetic_region %in% c('CLUSTAL_W','Cluster Analysis'), 'cluster',ifelse(genetic_region %in% c('ERIC', 'ERIC-PCR', 'ERIC_PCR','RFLP-ERIC'), 'ERIC',ifelse(genetic_region %in% c('IGS', 'IGS_PCR-RFLP', 'RFLP-IGS','IGSS', 'ITS', 'RFLP_ITS'), 'ITS',ifelse(genetic_region %in% c('nif_KD', 'nifD', 'nifD-K', 'nifh','nifH', 'nifH-nifDK', 'nifHD','RFLP_nifH'), 'nif',ifelse(genetic_region %in% c('nodBC', 'nodC', 'nodC-nodA','nodC-RFLP', 'nodA', 'nodA_PCR_RFLP','nodD', 'nodD1', 'nodD2', 'nodDAB','nodDF', 'nodF', 'nodY/K', 'RFLP_nodb3','RFLP_nodA', 'RFLP_nodb1', 'RFLP_nodb4','RFLP_nodb5', 'RFLP_nodC'), 'nod',ifelse(genetic_region %in% c('rep-PCR', 'REP_PCR', 'REP-PCR', 'rep_PCR', 'REP1R-I_REP2-I'), 'REP PCR',ifelse(genetic_region %in% c('PCR-RAPD', 'RAPD'), 'RAPD',ifelse(genetic_region %in% c('recA', 'recA-glnA-dnaK','recA, glnII', 'recA-glnII-atpD'), 'recA',ifelse(genetic_region %in% c('glnA', 'glnB', 'glnII', 'gltA', 'gryB','gyrA', 'gyrB'), 'gln','other')))))))))))))) %>% 
  # filter(genetic_region2=='16S') %>%
  group_by(plant_status, clean_name) %>% #each species tested across different studies is a replicate (averaged across gene regions for species with multiple)
  summarise(strain_richness=mean(strain_richness)) %>%  
  ungroup()



# make long version
nativeInvasiveGenetic_long <- plantData %>% 
  right_join(allSpp_3ab, by = c("clean_name", "genetic_region")) %>% 
  filter(num_nodules > 2) %>% 
  select(paper_id, clean_name, plant_status, annual_perennial, growth_form,
         habitat_type, sample_country, sample_continent, num_nodules,
         num_plants, genetic_region, strain_richness, cultivation.status) %>%
  mutate(
    plant_status = recode(plant_status, "introduced" = "non-native"),
    genetic_region2 = ifelse(genetic_region %in% c('16S', '16S-23S', '16S-ARDRA', '16S-IGS','16S-RFLP', '16S rDNA', '16S_23S_RFLP', '16S_ARDRA', '16S_BLAST', '16S_PCR_RFLP', '16S_rDNA', '16S_RFLP', 'ARDRA', 'PCR','PCR_RFLP', '16S_PCR-RFLP', 'RFLP_16S','RFLP', 'RFLP_PCR'), '16S',
                             ifelse(genetic_region %in% c('23S', '23S_IVS', 'RFLP-23S'), '23S',
                                    ifelse(genetic_region %in% c('AFLP', 'AFLP_Pst-A', 'AFLP_Pst-G', 'AFLP_Pst-GC'), 'AFLP',
                                           ifelse(genetic_region %in% c('BOX', 'BOX-AIR', 'box-PCR', 'BOX-PCR','Box_A1R-PCR', 'BOX_PCR', 'BoxA1R','BOXA1R', 'BOXA1R-PCR', 'BOXAIR'), 'BOX',
                                                  ifelse(genetic_region %in% c('CLUSTAL_W','Cluster Analysis'), 'cluster',
                                                         ifelse(genetic_region %in% c('ERIC', 'ERIC-PCR', 'ERIC_PCR','RFLP-ERIC'), 'ERIC',
                                                                ifelse(genetic_region %in% c('IGS', 'IGS_PCR-RFLP', 'RFLP-IGS','IGSS', 'ITS', 'RFLP_ITS'), 'ITS',
                                                                       ifelse(genetic_region %in% c('nif_KD', 'nifD', 'nifD-K', 'nifh','nifH', 'nifH-nifDK', 'nifHD','RFLP_nifH'), 'nif',
                                                                              ifelse(genetic_region %in% c('nodBC', 'nodC', 'nodC-nodA','nodC-RFLP', 'nodA', 'nodA_PCR_RFLP','nodD', 'nodD1', 'nodD2', 'nodDAB','nodDF', 'nodF', 'nodY/K', 'RFLP_nodb3','RFLP_nodA', 'RFLP_nodb1', 'RFLP_nodb4','RFLP_nodb5', 'RFLP_nodC'), 'nod',
                                                                                     ifelse(genetic_region %in% c('rep-PCR', 'REP_PCR', 'REP-PCR', 'rep_PCR', 'REP1R-I_REP2-I'), 'REP PCR',
                                                                                            ifelse(genetic_region %in% c('PCR-RAPD', 'RAPD'), 'RAPD',
                                                                                                   ifelse(genetic_region %in% c('recA', 'recA-glnA-dnaK','recA, glnII', 'recA-glnII-atpD'), 'recA',
                                                                                                          ifelse(genetic_region %in% c('glnA', 'glnB', 'glnII', 'gltA', 'gryB','gyrA', 'gyrB'), 'gln', 'other'))))))))))))))


#model
nativeInvasiveGeneticModel <- lmer(log(strain_richness) ~ plant_status + (1|clean_name),
                                   data=nativeInvasiveGenetic_long)
summary(nativeInvasiveGeneticModel)
anova(nativeInvasiveGeneticModel)

# true n values from long data
n_labels <- nativeInvasiveGenetic_long %>%
  group_by(clean_name, plant_status) %>%
  summarise(n = n(), .groups = "drop") %>%
  pivot_wider(names_from = plant_status, values_from = n, values_fill = 0) %>%
  mutate(
    species_label = paste0(
      gsub("_", " ", clean_name),
      " (", native, ", ", `non-native`, ")"
    )
  )
colnames(nativeInvasiveGenetic_long)

# mean richness per species within each status for plotting delta
diff_data <- nativeInvasiveGenetic_long %>%
  group_by(clean_name, plant_status) %>%
  summarise(strain_richness = mean(strain_richness, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = plant_status, values_from = strain_richness) %>%
  mutate(delta = `non-native` - native) %>%
  left_join(n_labels %>% select(clean_name, species_label), by = "clean_name")

diff_plot <- diff_data %>%
  mutate(direction = case_when(
    delta > 0 ~ "Non-native higher",
    delta < 0 ~ "Native higher",
    TRUE ~ "No difference"
  ))

xmax <- max(abs(diff_plot$delta), na.rm = TRUE)
xlim_use <- c(-(xmax + 0.5), (xmax + 0.5))

home_away_fig <- ggplot(diff_plot,
                        aes(x = delta, y = reorder(species_label, delta))) +
  geom_segment(aes(x = 0, xend = delta,
                   yend = reorder(species_label, delta)),
               linewidth = 1, color = "grey70") +
  geom_point(aes(fill = direction),
             shape = 21, color = "black", size = 4, stroke = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.9) +
  annotate("text", x = xlim_use[2], y = Inf, label = "Non-native higher",
           hjust = 1, vjust = 0.7, size = 5) +
  annotate("text", x = xlim_use[1], y = Inf, label = "Native higher",
           hjust = 0, vjust = 0.7, size = 5) +
  scale_fill_manual(values = c(
    "Native higher" = "#1b9e77",
    "Non-native higher" = "#d95f02",
    "No difference" = "grey80"
  )) +
  coord_cartesian(xlim = xlim_use, clip = "off") +
  labs(
    x = expression("Non-native " - " Native strain richness ("*Delta*")"),
    y = NULL
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    axis.text.y = element_text(face = "italic", size = 16),
    axis.title.x = element_text(size = 16),
    plot.margin = margin(8, 25, 8, 8)
  )

############################FIGURE 4 (HOME VS. AWAY - SPECIES RICHNESS DIFFERENCE)######

# standardize species names

standardize_species <- function(x) {
  x <- as.character(x)
  x <- str_trim(x)
  x <- str_replace_all(x, "\\s+", "_")
  x <- str_replace_all(x, "_L\\.?$", "")
  x <- str_replace_all(x, "_\\([^)]+\\)", "")
  x <- str_replace_all(x, "[,;:]", "")
  x <- str_replace_all(x, "__+", "_")
  x <- str_replace_all(x, "_$", "")
  return(x)
}

# strain dataset species

nativeInvasiveGenetic_long <- nativeInvasiveGenetic_long %>%
  mutate(clean_name_std = standardize_species(clean_name))

species_strain <- unique(nativeInvasiveGenetic_long$clean_name_std)


# identify rhizobia columns
meta_cols <- c(
  "paper_id","author","year","genus","species","genus_species",
  "introduced","invasive","cultivated","plant_status","strain_richness",
  "compares_homeaway","sample_country","sample_continent",
  "exo_NA","exo_SA","exo_AU","exo_AS","exo_EU","exo_AF",
  "nat_NA","nat_SA","nat_AU","nat_AS","nat_EU","nat_AF",
  "num_nodules","num_plants"
)
df_legume_rhizobia<- read_csv("legume_strain diversity_meta analysis presence absence edited.csv")
colnames(df_legume_rhizobia)
colnames(homeAwayAll)
unique(df_legume_rhizobia$plant_status)

# make join columns same type
df_legume_rhizobia <- df_legume_rhizobia %>%
  mutate(
    paper_id = as.character(paper_id),
    genus = as.character(genus),
    species = as.character(species)
  )
colnames(df_legume_rhizobia)
rhiz_cols <- setdiff(names(df_legume_rhizobia), meta_cols)


# clean rhizobia dataset

dat_rhiz <- df_legume_rhizobia

dat_rhiz <- dat_rhiz %>%
  filter(num_nodules > 2)

dat_rhiz <- dat_rhiz %>%
  mutate(
    plant_status = case_when(
      plant_status == "native" ~ "native",
      plant_status %in% c("introduced","non-native","invasive") ~ "non-native",
      TRUE ~ NA_character_
    )
  )

dat_rhiz <- dat_rhiz %>%
  mutate(clean_name_raw = genus_species)

dat_rhiz <- dat_rhiz %>%
  mutate(
    clean_name_raw = case_when(
      clean_name_raw == "Robinia_pseudoacacia L." ~ "Robinia_pseudoacacia",
      TRUE ~ clean_name_raw
    )
  )

dat_rhiz <- dat_rhiz %>%
  mutate(clean_name_std = standardize_species(clean_name_raw))

dat_rhiz <- dat_rhiz %>%
  filter(!is.na(plant_status), !is.na(clean_name_std))


# overlap check

species_overlap <- intersect(unique(dat_rhiz$clean_name_std), species_strain)

print(length(species_overlap))
print(head(species_overlap, 20))

dat_rhiz2 <- dat_rhiz %>%
  filter(clean_name_std %in% species_overlap)

# species in both native + non-native

species_check <- dat_rhiz2 %>%
  distinct(clean_name_std, plant_status)

species_check <- species_check %>%
  mutate(present = 1)

species_check <- species_check %>%
  pivot_wider(
    names_from = plant_status,
    values_from = present,
    values_fill = 0
  )

species_both <- species_check %>%
  filter(native == 1, `non-native` == 1)

species_both <- species_both$clean_name_std

print(length(species_both))

dat_rhiz2 <- dat_rhiz2 %>%
  filter(clean_name_std %in% species_both)


# long format rhizobia

rhiz_long <- dat_rhiz2 %>%
  pivot_longer(
    cols = all_of(rhiz_cols),
    names_to = "rhiz_species",
    values_to = "present"
  )

rhiz_long <- rhiz_long %>%
  mutate(present = as.numeric(present))

rhiz_long <- rhiz_long %>%
  filter(!is.na(present), present > 0)

# counts for labels

n_labels_rhiz <- dat_rhiz2 %>%
  distinct(clean_name_std, plant_status, paper_id)

n_labels_rhiz <- n_labels_rhiz %>%
  group_by(clean_name_std, plant_status)

n_labels_rhiz <- n_labels_rhiz %>%
  summarise(n = n(), .groups = "drop")

n_labels_rhiz <- n_labels_rhiz %>%
  pivot_wider(
    names_from = plant_status,
    values_from = n,
    values_fill = 0
  )

n_labels_rhiz <- n_labels_rhiz %>%
  mutate(
    species_label = paste0(
      str_replace_all(clean_name_std, "_", " "),
      " (", native, ", ", `non-native`, ")"
    )
  )


# richness
rhiz_richness <- rhiz_long %>%
  group_by(clean_name_std, plant_status)

rhiz_richness <- rhiz_richness %>%
  summarise(
    rhizobia_richness = n_distinct(rhiz_species),
    .groups = "drop"
  )

# delta

diff_data_rhiz <- rhiz_richness %>%
  pivot_wider(
    names_from = plant_status,
    values_from = rhizobia_richness,
    values_fill = 0
  )

diff_data_rhiz <- diff_data_rhiz %>%
  mutate(delta = `non-native` - native)

diff_data_rhiz <- diff_data_rhiz %>%
  left_join(
    n_labels_rhiz %>% select(clean_name_std, species_label),
    by = "clean_name_std"
  )

diff_data_rhiz <- diff_data_rhiz %>%
  mutate(
    Direction = case_when(
      delta > 0 ~ "Non-native higher",
      delta < 0 ~ "Native higher",
      TRUE ~ "No difference"
    )
  )

diff_data_rhiz <- diff_data_rhiz %>%
  arrange(delta)

diff_data_rhiz <- diff_data_rhiz %>%
  mutate(species_label = factor(species_label, levels = species_label))

print(nrow(diff_data_rhiz))


# plot

xmax <- max(abs(diff_data_rhiz$delta), na.rm = TRUE)
xlim_use <- c(-(xmax + 0.5), xmax + 0.5)

home_away_fig_rhiz <- ggplot(diff_data_rhiz, aes(x = delta, y = species_label)) +
  geom_segment(aes(x = 0, xend = delta, yend = species_label),
               linewidth = 1, color = "grey70") +
  geom_point(aes(fill = Direction),
             shape = 21, color = "black", size = 4) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c(
    "Native higher" = "#1b9e77",
    "Non-native higher" = "#d95f02",
    "No difference" = "grey80"
  )) +
  coord_cartesian(xlim = xlim_use) +
  labs(
    x = "Non-native - Native rhizobial richness",
    y = "Plant species"
  ) +theme_classic(base_size=13) + theme(axis.text.y = element_text(face = "italic")) 

print(home_away_fig_rhiz)

#########HOME VS AWAY (SPECIES OVERLAP) 

library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(grid)
# Rhizobial columns
rhiz_cols <- 29:314
rhiz_cols <- names(dat_rhiz)[rhiz_cols]
rhiz_cols <- rhiz_cols[!is.na(rhiz_cols)]
rhiz_cols <- rhiz_cols[rhiz_cols %in% names(dat_rhiz)]
rhiz_cols[rhiz_cols== "E_kummerowiae"] <- "Ensifer_kummerowiae"
#
# SPECIES ORDER (FROM DELTA FIGURE)

species_order <- c(
  "Robinia_pseudoacacia",
  "Mimosa_pudica",
  "Acacia_saligna",
  "Lupinus_angustifolius",
  "Acacia_mangium",
  "Acacia_longifolia",
  "Lathyrus_pratensis",
  "Galega_officinalis",
  "Astragalus_cicer",
  "Lotus_tenuis",
  "Medicago_polymorpha",
  "Medicago_laciniata",
  "Caragana_arborescens",
  "Sesbania_rostrata",
  "Cytisus_scoparius",
  "Acacia_pycnantha",
  "Lotus_corniculatus"
)

plot_dat <- dat_rhiz %>%
  mutate(clean_name = genus_species) %>%
  mutate(
    plant_status = case_when(
      plant_status == "native" ~ "native",
      plant_status %in% c("introduced","non-native","invasive") ~ "non-native",
      TRUE ~ NA_character_
    )
  ) %>%
  mutate(
    clean_name = case_when(
      clean_name == "Robinia_pseudoacacia L." ~ "Robinia_pseudoacacia",
      TRUE ~ clean_name
    )
  ) %>%
  filter(!is.na(plant_status)) %>%
  filter(!is.na(clean_name)) %>%
  filter(clean_name %in% species_order) %>%
  pivot_longer(
    cols = all_of(rhiz_cols),
    names_to = "rhiz_species",
    values_to = "present"
  )  %>% 
  filter(!is.na(present)) %>%
  filter(present %in% c(1,"1",TRUE,"TRUE","Present","present","X","x"))  %>% 
  mutate(rhiz_species = str_replace_all(rhiz_species,"_"," ")) %>%
  mutate(rhiz_species = str_squish(rhiz_species)) %>%
  mutate(rhiz_genus = word(rhiz_species,1)) %>% 
  filter(!is.na(rhiz_genus)) %>%
  filter(rhiz_genus != "") %>%
  distinct(clean_name,plant_status,rhiz_genus)

pair_class <- plot_dat %>%
  group_by(clean_name,rhiz_genus) %>%
  summarise(
    has_native = any(plant_status=="native"),
    has_nonnative = any(plant_status=="non-native"),
    .groups="drop"
  ) %>%
  mutate(
    assoc_group = case_when(
      has_native & !has_nonnative ~ "Native-only",
      has_native & has_nonnative ~ "Shared",
      !has_native & has_nonnative ~ "Non-native-only"
    )
  )

tax_order_df <- pair_class %>%
  count(assoc_group,rhiz_genus,name="freq") %>%
  mutate(
    assoc_group = factor(assoc_group,
                         levels=c("Native-only","Shared","Non-native-only"))
  ) %>%
  arrange(assoc_group,desc(freq),rhiz_genus) %>%
  mutate(
    taxon_plot = paste0(as.character(assoc_group),"___",rhiz_genus)
  )

# matrix
heat_dat <- pair_class %>%
  left_join(
    tax_order_df %>% select(rhiz_genus,assoc_group,taxon_plot),
    by=c("rhiz_genus","assoc_group")
  ) %>%
  mutate(present="Present")

all_combos <- expand_grid(
  clean_name = species_order,
  taxon_plot = tax_order_df$taxon_plot
)

heat_dat_full <- all_combos %>%
  left_join(
    heat_dat %>% select(clean_name,taxon_plot,present),
    by=c("clean_name","taxon_plot")
  ) %>%
  mutate(present = ifelse(is.na(present),"Absent",present)) %>%
  left_join(
    tax_order_df %>% select(taxon_plot,assoc_group,rhiz_genus),
    by="taxon_plot"
  ) %>%
  mutate(
    assoc_group = factor(assoc_group,
                         levels=c("Native-only","Shared","Non-native-only"))
  ) %>%
  mutate(taxon_plot = factor(taxon_plot,levels=tax_order_df$taxon_plot)) %>%
  mutate(clean_name = factor(clean_name,levels=rev(species_order)))

#Y-labels

species_counts <- dat_rhiz %>%
  mutate(clean_name = genus_species) %>%
  mutate(
    plant_status = case_when(
      plant_status=="native" ~ "native",
      plant_status %in% c("introduced","non-native","invasive") ~ "non-native",
      TRUE ~ NA_character_
    )
  ) %>%
  mutate(
    clean_name = case_when(
      clean_name=="Robinia_pseudoacacia L." ~ "Robinia_pseudoacacia",
      TRUE ~ clean_name
    )
  ) %>%
  filter(clean_name %in% species_order) %>%
  group_by(clean_name,plant_status) %>%
  summarise(n=n(),.groups="drop") %>%
  pivot_wider(
    names_from=plant_status,
    values_from=n,
    values_fill=0
  ) %>%
  mutate(
    species_label = paste0(
      str_replace_all(clean_name,"_"," "),
      " (",native,", ",`non-native`,")"
    )
  )

y_labels <- setNames(species_counts$species_label,species_counts$clean_name)

# x labels
x_labels <- heat_dat_full %>%
  distinct(taxon_plot,rhiz_genus) %>%
  arrange(taxon_plot) %>%
  {setNames(.$rhiz_genus,.$taxon_plot)}

#plot
heat_map <- ggplot(heat_dat_full,
                   aes(x=taxon_plot,y=clean_name,fill=present)) +
  geom_tile(color="white",linewidth=0.2) +
  facet_grid(.~assoc_group,scales="free_x",space="free_x") +
  scale_fill_manual(values=c("Absent"="grey92","Present"="black")) +
  scale_x_discrete(labels=x_labels) +
  scale_y_discrete(labels=y_labels) +
  labs(x="Rhizobial genera",y="Plant species") +
  theme_classic(base_size=13) +
  theme(
    strip.background=element_rect(fill="white",color="black"),
    axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size=9,face="italic"),
    axis.text.y=element_text(size=10,face="italic"),
    legend.position="top"
  )

print(heat_map)

ggsave("rhizobia_heatmap.png",heat_map,width=13,height=7,dpi=300)


# Combine the figures vertically
library(patchwork)
combined_homeaway <- (home_away_fig_rhiz /
                        heat_map) +
  plot_annotation(tag_levels = "A")

combined_homeaway


# Save the combined figure
ggsave(
  filename = "homeaway_strain_vs_species_richness.png",
  plot = combined_homeaway,
  width = 10,
  height = 10,
  units = "in",
  dpi = 300
)
