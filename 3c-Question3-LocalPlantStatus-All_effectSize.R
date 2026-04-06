
#### This script addresses question 3 in manuscript: LOCAL COMPARISON


library(ggplot2)
library(lme4)
library(lmerTest)
library(tidyverse)
library(dplyr)
library(readr)
library(olsrr)
library(tidyr)
library(vegan)
library(tibble)


####Read in Dataframes and Clean Data ####
### Starting a new analysis and this chunk is for loading files, Plant_Associations is the cleaned file
Clean_Species<-read_csv("legume_clean_names.csv") %>% 
  filter(old_genus_species!="") %>% 
  separate(old_genus_species,c("genus","species","subspecies","extra"),sep=" ") %>% 
  unite(col=genus_species, c(genus,species,subspecies,extra), sep='_',na.rm=TRUE)

Plant_Associations<-read.csv("legume_strain diversity_meta analysis_plant associations_edited names_presence absence.csv") %>% pivot_longer(names_to = "rhizobia_sp", values_to = "presence_absence", cols=Bradyrhizobium_nitroreducens:Paraburkholderia_caribensis)

##Now filtering out row crop and 999, 0 and 1 presence and absence and joining clean_species
Plant_Associations_Clean <- Plant_Associations %>%
  left_join(Clean_Species) %>%
  filter(cultivation.status!="row crop", genus_species!="999", presence_absence=="1", plant_status!="999", species!="sp",species!="spp",species!="sp.",paper_id!=142, plant_status!="extinct") %>% 
  #filter(!(sample_country %in% c("Japan-China","Japan_China","Kenya-Sudan","Czech_Republic-France-Georgia-Hungary-Italy-Romania-Spain","Senegal-Mauritania-Tunisia-Burundi","Malawi-Zambia-Kenya","Brazil-Venezuela"))) %>% 
  select(-introduced,-invasive,-cultivated) %>% 
  rename(paper_plant_status=plant_status) %>% 
  mutate(global_plant_status=ifelse((exo_NA+exo_SA+exo_AU+exo_AS+exo_EU+exo_AF)>0,1,0)) %>% 
  filter(num_nodules>2) %>% 
  mutate(global_plant_status=ifelse(global_plant_status==1,"nonnative","native")) %>% 
  mutate(paper_plant_status=ifelse(paper_plant_status=="invasive","nonnative",ifelse(paper_plant_status=="introduced","nonnative", paper_plant_status)))

colnames(Plant_Associations_Clean)


Plant_Data<-read.csv("legume_strain diversity_meta analysis_plant data.csv")

##Now filtering out row crop and 999, 0 and 1 presence and absence and joining clean_species
Plant_Data_Clean <- Plant_Data %>%
  rename(notes_plantdata=notes) %>% 
  left_join(Clean_Species) %>%
  filter(cultivation.status!="row crop", genus_species!="999", species!="sp",species!="spp",species!="sp.",paper_id!=142) %>% 
  filter(num_nodules>2) %>% 
  filter(num_nodules!=999)

colnames(Plant_Data_Clean)

#read in paper information dataframe 
Paper_Information<-read.csv("legume_strain diversity_meta analysis_paper information.csv")


#### Figure 5 Dataframe Creation ####
#r Create dataframe with just papers where native/non natives were studied within a given location

#Create a new dataframe that has only plant associations with species 
Plant_Status<-Plant_Associations_Clean %>% 
  select(paper_id,genus_species,new_name,cultivation.status,sample_country, global_plant_status,paper_plant_status) %>% 
  unique()

Native_NonNative_Fig1 <- Plant_Data_Clean %>%
  #merge plant associations so that we can determine plant status, etc.
  merge(Plant_Status,by = c("paper_id", "genus_species", "new_name", "cultivation.status","sample_country"),all=TRUE) %>% 
  select(paper_id,new_name,paper_plant_status,annual_perennial,growth_form,sample_country,sample_continent,num_nodules,genetic_region,strain_richness,compares_natinv) %>% 
  #only keep papers that looked at native and non native species within the same location
  filter(compares_natinv==1) %>% 
  na.omit(plant_paper_status) %>% 
  mutate(Graph_x=ifelse(paper_plant_status=="native","Native (n=53)","Non-native (n=17)")) %>% 
  #remove papers that no longer have a native or non-native partner after filtering steps
  filter(!(paper_id %in% c(289,325))) %>% 
  filter(new_name!="Prosopis chilensis" | sample_country!="Kenya") %>% 
  #average across plant ID, species, and status
  group_by(paper_id,new_name,sample_country,paper_plant_status,annual_perennial,growth_form,sample_continent,Graph_x) %>% 
  summarise(avg_strain_richness=mean(as.numeric(strain_richness))) %>% 
  ungroup()

#Calculate N number
Native_NonNative_N<-Native_NonNative_Fig1 %>% 
  group_by(Graph_x) %>% 
  summarise(status=length(new_name)) #54 native, 18 non native



###########################################LOCAL COMPARISION MODEL AND FIGURE##########
library(stringr)
library(glmmTMB)
library(car)

# Clean species names
Clean_Species <- read_csv("legume_clean_names.csv") %>% 
  filter(old_genus_species != "") %>% 
  separate(old_genus_species, c("genus", "species", "subspecies", "extra"), sep = " ") %>% 
  unite(col = genus_species, c(genus, species, subspecies, extra), sep = "_", na.rm = TRUE)


# Read and clean plant associations
Plant_Associations <- read.csv("legume_strain diversity_meta analysis_plant associations_edited names_presence absence.csv") %>% 
  pivot_longer(
    names_to = "rhizobia_sp",
    values_to = "presence_absence",
    cols = Bradyrhizobium_nitroreducens:Paraburkholderia_caribensis
  )

Plant_Associations_Clean <- Plant_Associations %>%
  left_join(Clean_Species, by = "genus_species") %>%
  filter(
    cultivation.status != "row crop",
    genus_species != "999",
    presence_absence == "1",
    plant_status != "999",
    species != "sp",
    species != "spp",
    species != "sp.",
    paper_id != 142,
    plant_status != "extinct"
  ) %>% 
  select(-introduced, -invasive, -cultivated) %>% 
  rename(paper_plant_status = plant_status) %>% 
  mutate(
    global_plant_status = ifelse((exo_NA + exo_SA + exo_AU + exo_AS + exo_EU + exo_AF) > 0, 1, 0),
    global_plant_status = ifelse(global_plant_status == 1, "non-native", "native"),
    paper_plant_status = ifelse(
      paper_plant_status %in% c("invasive", "introduced"),
      "non-native",
      paper_plant_status
    )
  ) %>% 
  filter(num_nodules > 2)


# Read and clean plant data
Plant_Data <- read.csv("legume_strain diversity_meta analysis_plant data.csv")

Plant_Data_Clean <- Plant_Data %>%
  rename(notes_plantdata = notes) %>% 
  left_join(Clean_Species, by = "genus_species") %>%
  filter(
    cultivation.status != "row crop",
    genus_species != "999",
    species != "sp",
    species != "spp",
    species != "sp.",
    paper_id != 142,
    num_nodules > 2,
    num_nodules != 999
  )

#Plant status table
Plant_Status <- Plant_Associations_Clean %>% 
  select(
    paper_id, genus_species, new_name, cultivation.status,
    sample_country, global_plant_status, paper_plant_status
  ) %>% 
  distinct()

#dataset for figure
Native_NonNative_Fig1 <- Plant_Data_Clean %>%
  merge(
    Plant_Status,
    by = c("paper_id", "genus_species", "new_name", "cultivation.status", "sample_country"),
    all = TRUE
  ) %>% 
  select(
    paper_id, new_name, paper_plant_status, annual_perennial, growth_form,
    sample_country, sample_continent, num_nodules, genetic_region,
    strain_richness, compares_natinv
  ) %>% 
  filter(compares_natinv == 1) %>% 
  filter(!is.na(paper_plant_status)) %>% 
  filter(!(paper_id %in% c(289, 325))) %>% 
  filter(!(new_name == "Prosopis chilensis" & sample_country == "Kenya")) %>% 
  group_by(
    paper_id, new_name, sample_country, paper_plant_status,
    annual_perennial, growth_form, sample_continent
  ) %>% 
  summarise(
    avg_strain_richness = mean(as.numeric(strain_richness), na.rm = TRUE),
    .groups = "drop"
  ) %>% 
  mutate(
    Graph_x = recode(
      paper_plant_status,
      "native" = "native",
      "non-native" = "non-native",
      "introduced" = "non-native",
      "invasive" = "non-native"
    ),
    Graph_x = factor(Graph_x, levels = c("native", "non-native"))
  ) %>%
  filter(!is.na(Graph_x))

colnames(Native_NonNative_Fig1)

#axis labels
Native_NonNative_N <- Native_NonNative_Fig1 %>% 
  count(Graph_x)

n_native <- Native_NonNative_N %>% filter(Graph_x == "native") %>% pull(n)
n_nonnative <- Native_NonNative_N %>% filter(Graph_x == "non-native") %>% pull(n)

#MODEL 
m_local <- glmmTMB(
  avg_strain_richness ~ Graph_x + (1 | new_name) + (1 | paper_id),
  data = Native_NonNative_Fig1,
  family = nbinom2(link = "log")
)

anova_res <- car::Anova(m_local)
print(summary(m_local))
print(anova_res)
res <- resid(m_local)
plot(fitted(m_local), res)
qqnorm(res)
plot(density(res))

library(DHARMa)

sim_res <- simulateResiduals(m_local)

plot(sim_res)

p_val <- anova_res$`Pr(>Chisq)`[1]
chisq_val <- anova_res$Chisq[1]

sig_lab <- case_when(
  p_val < 0.001 ~ "***",
  p_val < 0.01  ~ "**",
  p_val < 0.05  ~ "*",
  TRUE          ~ "ns"
)

p_lab <- ifelse(p_val < 0.001, "P < 0.001", paste0("P = ", signif(p_val, 2)))

#Figure positions
y_max  <- max(Native_NonNative_Fig1$avg_strain_richness, na.rm = TRUE)
y_line <- y_max * 1.08
y_text <- y_max * 1.15


# make sure factor order is correct
Native_NonNative_Fig1 <- Native_NonNative_Fig1 %>%
  mutate(
    Graph_x = recode(as.character(Graph_x),
                     "introduced" = "non-native",
                     "invasive" = "non-native"),
    Graph_x = factor(Graph_x, levels = c("native", "non-native"))
  )
colnames(Native_NonNative_Fig1)

# sample sizes for labels
Native_NonNative_N <- Native_NonNative_Fig1 %>%
  count(Graph_x)

n_native <- Native_NonNative_N %>% filter(Graph_x == "native") %>% pull(n)
n_nonnative <- Native_NonNative_N %>% filter(Graph_x == "non-native") %>% pull(n)


# y positions
y_max  <- max(Native_NonNative_Fig1$avg_strain_richness, na.rm = TRUE)
y_line <- y_max * 1.08
y_text <- y_max * 1.12

fig5 <- ggplot(
  Native_NonNative_Fig1,
  aes(x = Graph_x, y = avg_strain_richness, fill = Graph_x)
) +
  geom_boxplot(
    width = 0.6,
    outlier.shape = NA,
    alpha = 1,
    color = "black",
    linewidth = 0.8
  ) +
  geom_jitter(
    aes(color = Graph_x),
    width = 0.15,
    size = 2.2,
    alpha = 0.22,
    stroke = 0
  ) +
  annotate("segment", x = 1, xend = 2, y = y_line, yend = y_line, linewidth = 0.7) +
  annotate("segment", x = 1, xend = 1, y = y_line, yend = y_line * 0.98, linewidth = 0.7) +
  annotate("segment", x = 2, xend = 2, y = y_line, yend = y_line * 0.98, linewidth = 0.7) +
  annotate("text", x = 1.5, y = y_text, label = sig_lab, size = 6, fontface = "bold") +
  scale_x_discrete(
    labels = c(
      "native" = paste0("Native\n(n=", n_native, ")"),
      "non-native" = paste0("Non-native\n(n=", n_nonnative, ")")
    )
  ) +
  scale_fill_manual(
    values = c(
      "native" = "#1b9e77",
      "non-native" = "#d95f02"
    )
  ) +
  scale_color_manual(
    values = c(
      "native" = "#1b9e77",
      "non-native" = "#d95f02"
    )
  ) +
  labs(
    x = "Local plant status",
    y = "Rhizobial strain richness"
  ) +
  expand_limits(y = y_text * 1.03) +
  theme_classic(base_size = 16) +
  theme(
    legend.position = "none",
    axis.text = element_text(color = "black", size = 13),
    axis.title = element_text(color = "black", size = 15),
    plot.margin = margin(10, 14, 10, 10)
  )

fig5


#######rhizobial richness#######

#clean species name
Clean_Species <- read_csv("legume_clean_names.csv") %>% 
  filter(old_genus_species != "") %>% 
  separate(old_genus_species, c("genus", "species", "subspecies", "extra"), sep = " ") %>% 
  unite(col = genus_species, c(genus, species, subspecies, extra), sep = "_", na.rm = TRUE)


Plant_Associations <- read.csv("legume_strain diversity_meta analysis_plant associations_edited names_presence absence.csv") %>% 
  pivot_longer(
    names_to = "rhizobia_sp",
    values_to = "presence_absence",
    cols = Bradyrhizobium_nitroreducens:Paraburkholderia_caribensis
  )

Plant_Associations_Clean <- Plant_Associations %>%
  left_join(Clean_Species, by = "genus_species") %>%
  filter(
    cultivation.status != "row crop",
    genus_species != "999",
    presence_absence == "1",
    plant_status != "999",
    species != "sp",
    species != "spp",
    species != "sp.",
    paper_id != 142,
    plant_status != "extinct"
  ) %>% 
  select(-introduced, -invasive, -cultivated) %>% 
  rename(paper_plant_status = plant_status) %>% 
  mutate(
    global_plant_status = ifelse((exo_NA + exo_SA + exo_AU + exo_AS + exo_EU + exo_AF) > 0, 1, 0),
    global_plant_status = ifelse(global_plant_status == 1, "non-native", "native"),
    paper_plant_status = ifelse(
      paper_plant_status %in% c("invasive", "introduced"),
      "non-native",
      paper_plant_status
    )
  ) %>% 
  filter(num_nodules > 2)

Plant_Data <- read.csv("legume_strain diversity_meta analysis_plant data.csv")

Plant_Data_Clean <- Plant_Data %>%
  rename(notes_plantdata = notes) %>% 
  left_join(Clean_Species, by = "genus_species") %>%
  filter(
    cultivation.status != "row crop",
    genus_species != "999",
    species != "sp",
    species != "spp",
    species != "sp.",
    paper_id != 142,
    num_nodules > 2,
    num_nodules != 999
  )

Plant_Status <- Plant_Associations_Clean %>% 
  select(
    paper_id, genus_species, new_name, cultivation.status,
    sample_country, global_plant_status, paper_plant_status
  ) %>% 
  distinct()


Rhizobia_With_Status <- Plant_Associations_Clean %>%
  dplyr::select(
    paper_id, genus_species, new_name, cultivation.status,
    sample_country, sample_continent, rhizobia_sp
  ) %>%
  distinct() %>%
  left_join(
    Plant_Data_Clean %>%
      select(
        paper_id, genus_species, new_name, cultivation.status,
        sample_country, compares_natinv
      ) %>%
      distinct(),
    by = c("paper_id", "genus_species", "new_name",
           "cultivation.status", "sample_country")
  ) %>%
  left_join(
    Plant_Status,
    by = c("paper_id", "genus_species", "new_name",
           "cultivation.status", "sample_country")
  )

Native_NonNative_RhizFig <- Rhizobia_With_Status %>%
  filter(compares_natinv == 1) %>%
  filter(!is.na(paper_plant_status)) %>%
  filter(!(paper_id %in% c(289, 325))) %>%
  filter(!(new_name == "Prosopis chilensis" & sample_country == "Kenya")) %>%
  group_by(
    paper_id, new_name, paper_plant_status, sample_country, sample_continent
  ) %>%
  summarise(
    avg_rhizobial_species_richness = n_distinct(rhizobia_sp),
    .groups = "drop"
  ) %>%
  mutate(
    Graph_x = recode(
      paper_plant_status,
      "native" = "native",
      "non-native" = "non-native",
      "introduced" = "non-native",
      "invasive" = "non-native"
    ),
    Graph_x = factor(Graph_x, levels = c("native", "non-native"))
  ) %>%
  filter(!is.na(Graph_x))

Native_NonNative_Rhiz_N <- Native_NonNative_RhizFig %>% count(Graph_x)

n_native <- Native_NonNative_Rhiz_N %>% filter(Graph_x == "native") %>% pull(n)
n_nonnative <- Native_NonNative_Rhiz_N %>% filter(Graph_x == "non-native") %>% pull(n)

m_local_rhiz <- glmmTMB(
  avg_rhizobial_species_richness ~ Graph_x + (1 | new_name) + (1 | paper_id),
  data = Native_NonNative_RhizFig,
  family = nbinom2(link = "log")
)

anova_res <- car::Anova(m_local_rhiz)
p_val <- anova_res$`Pr(>Chisq)`[1]

sig_lab <- case_when(
  p_val < 0.001 ~ "***",
  p_val < 0.01  ~ "**",
  p_val < 0.05  ~ "*",
  TRUE          ~ "ns"
)

y_max  <- max(Native_NonNative_RhizFig$avg_rhizobial_species_richness, na.rm = TRUE)
y_line <- y_max * 1.08
y_text <- y_max * 1.12

fig_rhiz_local <- ggplot(
  Native_NonNative_RhizFig,
  aes(x = Graph_x, y = avg_rhizobial_species_richness, fill = Graph_x)
) +
  geom_boxplot(
    width = 0.6,
    outlier.shape = NA,
    alpha = 1,
    color = "black",
    linewidth = 0.8
  ) +
  geom_jitter(
    aes(color = Graph_x),
    width = 0.15,
    size = 2.2,
    alpha = 0.22,
    stroke = 0
  ) +
  annotate("segment", x = 1, xend = 2, y = y_line, yend = y_line, linewidth = 0.7) +
  annotate("segment", x = 1, xend = 1, y = y_line, yend = y_line * 0.98, linewidth = 0.7) +
  annotate("segment", x = 2, xend = 2, y = y_line, yend = y_line * 0.98, linewidth = 0.7) +
  annotate("text", x = 1.5, y = y_text, label = sig_lab, size = 6, fontface = "bold") +
  scale_x_discrete(
    labels = c(
      "native" = paste0("Native\n(n=", n_native, ")"),
      "non-native" = paste0("Non-native\n(n=", n_nonnative, ")")
    )
  ) +
  scale_fill_manual(values = c(
    "native" = "#1b9e77",
    "non-native" = "#d95f02"
  )) +
  scale_color_manual(values = c(
    "native" = "#1b9e77",
    "non-native" = "#d95f02"
  )) +
  labs(
    x = "Local plant status",
    y = "Rhizobial species richness"
  ) +
  expand_limits(y = y_text * 1.03) +
  theme_classic(base_size = 16) +
  theme(
    legend.position = "none",
    axis.text = element_text(color = "black", size = 13),
    axis.title = element_text(color = "black", size = 15),
    plot.margin = margin(10, 14, 10, 10)
  )

fig_rhiz_local

#########################################ORDINATION############################

Rhizobia_Ordination_Data <- Plant_Associations_Clean %>%
  mutate(
    Graph_x = recode(
      paper_plant_status,
      "native" = "native",
      "introduced" = "non-native",
      "invasive" = "non-native",
      "non-native" = "non-native"
    ),
    Graph_x = factor(Graph_x, levels = c("native", "non-native"))
  ) %>%
  inner_join(
    Native_NonNative_RhizFig %>%
      select(paper_id, new_name, sample_country, sample_continent, Graph_x) %>%
      distinct(),
    by = c("paper_id", "new_name", "sample_country", "sample_continent", "Graph_x")
  ) %>%
  select(paper_id, new_name, sample_country, sample_continent, Graph_x, rhizobia_sp) %>%
  distinct() %>%
  mutate(sample_id = paste(paper_id, new_name, sample_country, sample_continent, Graph_x, sep = "__"))

# check TRUE sample counts
Rhizobia_Ordination_Data %>%
  distinct(sample_id, Graph_x) %>%
  count(Graph_x)

# Community matrix

comm_pa <- Rhizobia_Ordination_Data %>%
  mutate(presence = 1) %>%
  select(sample_id, rhizobia_sp, presence) %>%
  distinct() %>%
  pivot_wider(
    names_from = rhizobia_sp,
    values_from = presence,
    values_fill = 0
  )

meta_ord <- Rhizobia_Ordination_Data %>%
  distinct(sample_id, paper_id, new_name, sample_country, sample_continent, Graph_x)

comm_mat <- comm_pa %>%
  column_to_rownames("sample_id") %>%
  as.matrix()

# remove empty rows if any
keep_rows <- rowSums(comm_mat) > 0
comm_mat <- comm_mat[keep_rows, , drop = FALSE]
meta_ord <- meta_ord %>% filter(sample_id %in% rownames(comm_mat))
meta_ord <- meta_ord %>% slice(match(rownames(comm_mat), sample_id))

# PCoA

dist_jac <- vegdist(comm_mat, method = "jaccard", binary = TRUE)
pcoa_res <- cmdscale(dist_jac, eig = TRUE, k = 2)

ord_dat <- as.data.frame(pcoa_res$points)
colnames(ord_dat) <- c("PCoA1", "PCoA2")
ord_dat$sample_id <- rownames(ord_dat)

ord_dat <- ord_dat %>%
  left_join(meta_ord, by = "sample_id")

eig_vals <- pcoa_res$eig
var_exp <- round(100 * eig_vals / sum(abs(eig_vals)), 1)
ax1_lab <- paste0("PCoA1 (", var_exp[1], "%)")
ax2_lab <- paste0("PCoA2 (", var_exp[2], "%)")

#PERMANOVA
set.seed(123)
perm <- adonis2(dist_jac ~ Graph_x, data = meta_ord)
print(perm)

perm_lab <- paste0(
  "PERMANOVA: R² = ", round(perm$R2[1], 3),
  ", P = ", signif(perm$`Pr(>F)`[1], 2)
)

# Plot
ord_fig <- ggplot(ord_dat, aes(x = PCoA1, y = PCoA2, color = Graph_x, fill = Graph_x)) +
  stat_ellipse(
    geom = "polygon",
    alpha = 0.12,
    level = 0.68,
    linewidth = 0.6
  ) +
  geom_point(size = 3, alpha = 0.85) +
  scale_color_manual(values = c(
    "native" = "#1b9e77",
    "non-native" = "#d95f02"
  )) +
  scale_fill_manual(values = c(
    "native" = "#1b9e77",
    "non-native" = "#d95f02"
  )) +
  labs(
    x = ax1_lab,
    y = ax2_lab,
    color = "Local plant status",
    fill = "Local plant status",
    subtitle = perm_lab
  ) +
  theme_classic(base_size = 15) +
  theme(
    legend.position = "right",
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    plot.title = element_text(face = "bold")
  )

ord_fig


########PIECHART FOR PROPORTION OF SIMILAR VS DIFFERENT STRAINS)#####
Strain_overlap <- read_csv("Strain_overlap.csv")
colnames(Strain_overlap)

pie_data <- Strain_overlap %>% library(dplyr)
library(tidyr)
library(ggplot2)

# clean and check scale
tmp <- Strain_overlap %>%
  mutate(
    Percent_InCommon = as.numeric(Percent_InCommon),
    Percent_Different = as.numeric(Percent_Different)
  )

mean_in_common <- mean(tmp$Percent_InCommon, na.rm = TRUE)
mean_different <- mean(tmp$Percent_Different, na.rm = TRUE)

# if values are proportions (<=1), convert to %
if (mean_in_common <= 1 && mean_different <= 1) {
  mean_in_common <- mean_in_common * 100
  mean_different <- mean_different * 100
}

pie_data <- data.frame(
  category = c("In common", "Different"),
  percent = c(mean_in_common, mean_different)
) %>%
  mutate(percent = round(percent, 1))

pie_data

pie_data
pie_fig <- ggplot(pie_data, aes(x = "", y = percent, fill = category)) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y") +
  geom_text(
    aes(label = paste0(percent, "%")),
    position = position_stack(vjust = 0.5),
    color = "white",
    size = 6,
    fontface = "bold"
  ) +
  scale_fill_manual(
    values = c(
      "In common" = "#3D4B33",
      "Different" = "#6F7F97"
    )
  ) +
  theme_void(base_size = 14) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

pie_fig


##########combine figures###
library(patchwork)
fig5 <- fig5 + theme(legend.position = "none")
fig_rhiz_local <- fig_rhiz_local + theme(legend.position = "none")

combined_fig <- (
  fig5 + pie_fig
) /
  ( ord_fig
  ) +
  plot_annotation(tag_levels = "A")

combined_fig

ggsave(
  filename = "combined_legume_rhizobia_figure.png",
  plot = combined_fig,
  width = 12,
  height = 10,
  units = "in",
  dpi = 300
)
