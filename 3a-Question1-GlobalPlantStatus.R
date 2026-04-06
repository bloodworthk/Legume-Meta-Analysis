#### Legume Meta Analysis ####

library(car)
library(lme4)
library(lmerTest)
library(ggsci)
library(ggrepel)
library(tidyverse)
library(dplyr)
library(indicspecies)
library(tibble)
library(ggplot2)
library(permute)
library(glmmTMB)
library(stringr)

##### data import #####

# getting clean legume names for standardizing across papers
cleanNames <- read.csv('legume_clean_names.csv') %>% 
  separate(old_genus_species,c("genus","species","subspecies","extra"),sep=" ") %>% 
  unite(col=genus_species, c(genus,species,subspecies,extra), sep='_',na.rm=TRUE) %>% 
  dplyr::rename(clean_name=new_name) %>% 
  dplyr::select(-notes)

# getting native/invasive status for each study and legume species
plantStatus <- read.csv('legume_strain diversity_meta analysis_plant associations_edited names_presence absence.csv') %>% 
  left_join(cleanNames) %>%
  # rename(paper_id=ï..paper_id) %>%
  dplyr::select(paper_id, clean_name, plant_status, sample_country) %>% 
  unique()

# removing species where legume species isn't clear, native/invasive status is unknown, and cultivated species
plantData <- read.csv("legume_strain diversity_meta analysis_plant data.csv") %>% 
  #change all 999 to NAs
  mutate(across(all_of(c('sites_sampled', 'num_nodules', 'strain_richness', 'proportion_novel_strains',
                         'proportion_familiar_strains', 'proportion_overlapping_strains', 'num_spp_compared')), ~na_if(., 999))) %>% 
  left_join(cleanNames) %>% 
  dplyr::select(-notes, -genus_species) %>% 
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

globalStatus <- read.csv('legume_strain diversity_meta analysis_plant associations_edited names_presence absence.csv') %>% 
  # rename(paper_id=ï..paper_id) %>%
  left_join(cleanNames) %>% 
  mutate(global_plant_status=ifelse((exo_NA+exo_SA+exo_AU+exo_AS+exo_EU+exo_AF)>0,'introduced','native')) %>% 
  dplyr::select(paper_id, clean_name, global_plant_status) %>% 
  unique()

homeAwayAll <- plantData %>% 
  filter(num_nodules>2) %>%
  filter(!is.na(clean_name)) %>% 
  left_join(globalStatus) %>% 
  # group_by(global_plant_status, clean_name) %>%  #averaging over papers, gene regions, etc -- plant species are our "replicates" for this question
  # summarise(strain_richness=mean(strain_richness)) %>% 
  # ungroup() %>% 
  filter(!is.na(strain_richness)) %>% 
  filter(!(strain_source %in% c('999', 'rhizosphere_samples', 'NA', 'Look up Vincent (1970)')))

#### Determine number of countries, species, etc. ####

unique(sort(homeAwayAll$sample_country))

NonNative_Species<-homeAwayAll %>% 
  filter(global_plant_status=="introduced")

unique(sort(NonNative_Species$clean_name))
unique(sort(NonNative_Species$sample_country))

Native_Species<-homeAwayAll %>% 
  filter(global_plant_status=="native")

unique(sort(Native_Species$clean_name))
unique(sort(Native_Species$sample_country))

#local plant status
NonNative_Species_Local<-homeAwayAll %>% 
  filter(plant_status=="introduced")

unique(sort(NonNative_Species_Local$clean_name))
unique(sort(NonNative_Species_Local$sample_country))

Native_Species_Local<-homeAwayAll %>% 
  filter(plant_status=="native")

unique(sort(Native_Species_Local$clean_name))
unique(sort(Native_Species_Local$sample_country))

##### mixed model for native/invasive globally #####

#normality
hist(homeAwayAll$strain_richness)
qqPlot(homeAwayAll$strain_richness)
shapiro.test(homeAwayAll$strain_richness)
# W = 0.70895, p-value < 2.2e-16

hist(log(homeAwayAll$strain_richness))
qqPlot(log(homeAwayAll$strain_richness))
shapiro.test(log(homeAwayAll$strain_richness))
# W = 0.9773, p-value = 2.277e-05




#################GLOBAL COMPARISION-STRAIN RICHNESS (MODEL + FIGURE 3A)
x=glmmTMB(strain_richness ~ global_plant_status + (1|clean_name) + (1|genetic_region),
                  data = homeAwayAll, family = nbinom2(link = "log"))

summary(x)
Anova(x)
res <- resid(x)
plot(fitted(x), res)
qqnorm(res)
plot(density(res))

homeAwayAll <- homeAwayAll %>%
  mutate(
    global_plant_status = factor(
      global_plant_status,
      levels = c("native", "introduced"),
      labels = c("native", "non-native")
    )
  )

# choose y positions for significance line/text
y_max  <- max(homeAwayAll$strain_richness, na.rm = TRUE)
y_line <- y_max * 1.08
y_text <- y_max * 1.13

fig_global_strain_richness <- ggplot(
  homeAwayAll,
  aes(x = global_plant_status, y = strain_richness, fill = global_plant_status)
) +
  geom_boxplot(
    width = 0.6,
    outlier.shape = NA,
    alpha = 0.9,
    color = "black",
    linewidth = 0.5
  ) +
  geom_jitter(
    width = 0.12,
    size = 1.8,
    alpha = 0.18,
    shape = 21,
    color = "black",
    stroke = 0.2
  ) +
  annotate("segment", x = 1, xend = 2, y = y_line, yend = y_line, linewidth = 0.5) +
  annotate("segment", x = 1, xend = 1, y = y_line, yend = y_line * 0.98, linewidth = 0.5) +
  annotate("segment", x = 2, xend = 2, y = y_line, yend = y_line * 0.98, linewidth = 0.5) +
  annotate("text", x = 1.5, y = y_text, label = "*", size = 6) +
  scale_fill_manual(
    values = c(
      "native" = "#1b9e77",
      "non-native" = "#d95f02"
    )
  ) +
  scale_x_discrete(
    labels = c(
      paste0("Native\n(n=", sum(homeAwayAll$global_plant_status == "native", na.rm = TRUE), ")"),
      paste0("Non-native\n(n=", sum(homeAwayAll$global_plant_status == "non-native", na.rm = TRUE), ")")
    )
  ) +
  labs(
    x = "Global plant status",
    y = "Rhizobial strain richness"
  ) +
  expand_limits(y = y_text * 1.03) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14),
    plot.margin = margin(10, 15, 10, 10)
  )


########################GLOBAL COMPARISION - SPECIES RICHNESS
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
homeAwayAll <- homeAwayAll %>%
  mutate(
    paper_id = as.character(paper_id),
    genus = as.character(genus),
    species = as.character(species)
  )

# merge 
merged_df <- df_legume_rhizobia %>%
  left_join(
    homeAwayAll,
    by = c("paper_id", "genus", "species"),
    suffix = c("_rhiz", "_home"),
    relationship = "many-to-many"
  ) %>%
  mutate(
    author = coalesce(as.character(author_home), as.character(author_rhiz)),
    year = coalesce(as.character(year_home), as.character(year_rhiz)),
    sample_country = coalesce(as.character(sample_country_home), as.character(sample_country_rhiz)),
    sample_continent = coalesce(as.character(sample_continent_home), as.character(sample_continent_rhiz)),
    num_nodules = coalesce(as.character(num_nodules_home), as.character(num_nodules_rhiz)),
    num_plants = coalesce(as.character(num_plants_home), as.character(num_plants_rhiz)),
    strain_richness = coalesce(as.character(strain_richness_home), as.character(strain_richness_rhiz)),
    compares_homeaway = coalesce(as.character(compares_homeaway_home), as.character(compares_homeaway_rhiz)),
    plant_status = coalesce(as.character(plant_status_home), as.character(plant_status_rhiz))
  ) %>%
  select(-ends_with("_rhiz"), -ends_with("_home"))


# rhizobia columns
rhizobia_mat <- merged_df %>%
  select(Achromobacter_sp:Paraburkholderia_caribensis)

# presence-absence
rhizobia_pa <- rhizobia_mat %>%
  mutate(across(everything(), ~ ifelse(is.na(.), 0, ifelse(. > 0, 1, 0))))

# grouping variable
group <- merged_df$global_plant_status

# remove missing group rows
keep <- !is.na(group)
rhizobia_pa <- rhizobia_pa[keep, , drop = FALSE]
group <- droplevels(as.factor(group[keep]))


rhizo_cols <- colnames(merged_df)[20:304]

richness_df <- merged_df %>%
  select(global_plant_status, clean_name, genetic_region, all_of(rhizo_cols)) %>%
  pivot_longer(
    cols = all_of(rhizo_cols),
    names_to = "rhizobia_species",
    values_to = "presence"
  ) %>%
  filter(presence > 0) %>%
  group_by(clean_name, global_plant_status, genetic_region) %>%
  summarise(
    rhizobia_richness = n_distinct(rhizobia_species),
    .groups = "drop"
  ) %>%
  filter(!is.na(global_plant_status))   # remove NA
colnames(richness_df)

y=glmmTMB(rhizobia_richness~ global_plant_status + (1|clean_name) + (1|genetic_region),
          data = richness_df, family = nbinom2(link = "log"))

summary(y)
Anova(y)

ggplot(richness_df, aes(x = global_plant_status, y = rhizobia_richness, fill = global_plant_status)) +
  geom_boxplot(width = 0.6, alpha = 0.7) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.6) +
  scale_fill_manual(values = c(
    "native" = "forestgreen",
    "introduced" = "orange"
  )) +
  labs(
    x = "Plant status",
    y = "Rhizobial strain richness"
  ) +
  theme_classic(base_size = 14)


mean_summary <- richness_df %>%
  group_by(global_plant_status) %>%
  summarise(
    mean_richness = mean(rhizobia_richness, na.rm = TRUE),
    sd = sd(rhizobia_richness, na.rm = TRUE),
    n = n(),
    se = sd / sqrt(n)
  )

mean_summary


# ###########################indicator species analysis for publication new

# check structure
str(rhizobia_pa)

rhizobia_pa2 <- rhizobia_pa %>%
  as.data.frame()

if ("clean_name" %in% colnames(rhizobia_pa2)) {
  rownames(rhizobia_pa2) <- rhizobia_pa2$clean_name
  rhizobia_pa2$clean_name <- NULL
}

# convert every column to numeric
rhizobia_pa2[] <- lapply(rhizobia_pa2, function(x) as.numeric(as.character(x)))

colSums(is.na(rhizobia_pa2))

# replace NAs with 0 
rhizobia_pa2[is.na(rhizobia_pa2)] <- 0

rhizobia_mat <- as.matrix(rhizobia_pa2)

# confirm
str(rhizobia_mat)
mode(rhizobia_mat)


# run indicator species analysis
set.seed(123)
ind <- multipatt(
  rhizobia_mat,
  group,
  func = "IndVal.g",
  control = how(nperm = 999)
)

summary(ind)
# Extract results 
ind_table <- ind$sign %>%
  as.data.frame() %>%
  rownames_to_column("rhizobia_species")

group_cols <- colnames(ind_table)[grepl("^s\\.", colnames(ind_table))]

ind_table <- ind_table %>%
  rowwise() %>%
  mutate(
    associated_with = {
      vals <- c_across(all_of(group_cols))
      vals[is.na(vals)] <- 0
      hit <- group_cols[vals == 1]
      if (length(hit) == 1) gsub("^s\\.", "", hit)
      else if (length(hit) > 1) "multiple"
      else "none"
    }
  ) %>%
  ungroup()

# 10 species for native vs non-native
ind_table_balanced <- ind_table %>%
  filter(p.value < 0.05, associated_with %in% c("native", "introduced")) %>%
  mutate(
    rhizobia_species_clean = rhizobia_species,
    rhizobia_species_clean = str_replace_all(rhizobia_species_clean, "_", " "),
    rhizobia_species_clean = str_replace(rhizobia_species_clean, "\\.\\.\\..*", ""),
    p_label = paste0("p = ", formatC(p.value, format = "f", digits = 3))
  ) %>%
  group_by(associated_with) %>%
  slice_max(order_by = stat, n = 10, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    associated_with = recode(associated_with,
                             "introduced" = "Non-native"),
    associated_with = recode(associated_with,
                             "native" = "Native"),
    associated_with = factor(associated_with,
                             levels = c("Native", "Non-native"))
  )

indicator_species_fig <- ggplot(
  ind_table_balanced,
  aes(x = reorder(rhizobia_species_clean, stat), y = stat, fill = associated_with)
) +
  geom_col(width = 0.75, color = "black", linewidth = 0.3) +
  coord_flip() +
  facet_wrap(~ associated_with, scales = "free_y") +
  scale_fill_manual(values = c(
    "Native" = "#1b9e77",
    "Non-native" = "#d95f02"
  )) +
  expand_limits(y = max(ind_table_balanced$stat) * 1.18) +
  labs(
    x = "Rhizobial species",
    y = "Indicator value",
    fill = "Global Plant Status"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.y = element_text(size = 11, face = "italic"),
    legend.position = "none",
    strip.text = element_text(size = 12, face = "bold")
  )

indicator_species_fig


#####Combine figures

fig_global_strain_richness <- fig_global_strain_richness +
  theme(
    axis.title.y = element_text(margin = margin(r = 6)),
    plot.margin = margin(5.5, 5.5, 5.5, 5.5)
  )
library(patchwork)

combined_fig <- free(fig_global_strain_richness, side = "l") /
  indicator_species_fig +
  plot_layout(heights = c(1, 1.4)) +
  plot_annotation(tag_levels = "A")

combined_fig
ggsave(
  "Fig3_indicator_combined.png",
  combined_fig,
  width = 9,
  height = 10,
  dpi = 600
)
