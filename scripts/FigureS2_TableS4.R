here::i_am("scripts/analysis/zim_brief_report/FigureS2_TableS4.R")
library(here)
library(tidyverse)

run <- "zim_brief_report"
comb_data_meta_rand_id <- readRDS(here("data/processed",run,"dr_processed/combined_dr_data.rds"))

prev_table_ZAN_travel <- comb_data_meta_rand_id %>%
  filter(region == "ZB" & travel == 1) %>%
  group_by(mutation_name) %>%
  summarise(total_ALT_Z_t = sum(genotype >= 1, na.rm = TRUE),
            total_REF_Z_t = sum(genotype == 0, na.rm = TRUE),
            N_Z_t = total_ALT_Z_t + total_REF_Z_t,
            prevalence_Z_t = total_ALT_Z_t/N_Z_t) %>%
  mutate(
    CI = pmap(list(total_ALT_Z_t, N_Z_t), ~ binom.confint(x = ..1, n = ..2, methods = "exact")),
    ci.lower_Z_t = map_dbl(CI, 'lower'),
    ci.upper_Z_t = map_dbl(CI, 'upper')
  ) %>%
  select(-CI)

prev_table_ZAN_notravel <- comb_data_meta_rand_id %>%
  filter(region == "ZB" & travel == 0) %>%
  group_by(mutation_name) %>%
  summarise(total_ALT_Z_nt = sum(genotype >= 1, na.rm = TRUE),
            total_REF_Z_nt = sum(genotype == 0, na.rm = TRUE),
            N_Z_nt = total_ALT_Z_nt + total_REF_Z_nt,
            prevalence_Z_nt = total_ALT_Z_nt/N_Z_nt) %>%
  mutate(
    CI = pmap(list(total_ALT_Z_nt, N_Z_nt), ~ binom.confint(x = ..1, n = ..2, methods = "exact")),
    ci.lower_Z_nt = map_dbl(CI, 'lower'),
    ci.upper_Z_nt = map_dbl(CI, 'upper')
  ) %>%
  select(-CI)

prev_table_TAN <- comb_data_meta_rand_id %>%
  filter(region == "TZ") %>%
  group_by(mutation_name) %>%
  summarise(total_ALT_M = sum(genotype >= 1, na.rm = TRUE),
            total_REF_M = sum(genotype == 0, na.rm = TRUE),
            N_M = total_ALT_M + total_REF_M,
            prevalence_M = total_ALT_M/N_M) %>%
  mutate(
    CI = pmap(list(total_ALT_M, N_M), ~ binom.confint(x = ..1, n = ..2, methods = "exact")),
    ci.lower_M = map_dbl(CI, 'lower'),
    ci.upper_M = map_dbl(CI, 'upper')
  ) %>%
  select(-CI)

reformat <-
  left_join(prev_table_ZAN_travel,
            prev_table_ZAN_notravel,
            by = "mutation_name") %>%
  left_join(.,prev_table_TAN,by = "mutation_name") %>% 
  mutate_if(is.numeric, round,3) %>%
  unite(CI_Z_t,ci.lower_Z_t,ci.upper_Z_t,remove = T, sep = "-") %>%
  unite(CI_Z_nt,ci.lower_Z_nt,ci.upper_Z_nt,remove = T, sep = "-") %>% 
  unite(CI_M,ci.lower_M,ci.upper_M,remove = T, sep = "-")

# Figure S2 ---------------------------------------------------------------
order <- c("CRT-K76T","DHFR-N51I","DHFR-C59R",
           "DHFR-S108N","DHPS-A437G","DHPS-K540E",
           "DHPS-A581G","MDR1-N86Y","MDR1-Y184F",
           "MDR1-D1246Y")
mutations_for_plt <- c("crt-K76T",
                       "dhfr-ts-S108N",
                       "dhfr-ts-N51I",
                       "dhfr-ts-C59R",
                       "dhps-A437G",
                       "dhfr-ts-S108T",
                       "dhps-A581G", 
                       "dhps-K540E", 
                       "mdr1-N86Y", 
                       "mdr1-D1246Y", 
                       "mdr1-Y184F",
                       "mdr2-I492V")
bar_graph_fig2 <- reformat %>% 
  filter(mutation_name %in% mutations_for_plt) %>% 
  select(c("mutation_name",starts_with("prev"),starts_with("CI"))) %>% 
  pivot_longer(cols = starts_with("prevalence_"),
               names_to = "group",
               names_prefix = "prevalence_",
               values_to = "prevalence") %>%
  pivot_longer(cols = starts_with("CI_"),
               names_to = "group_ci",
               names_prefix = "CI_",
               values_to = "CI") %>%
  filter(group == group_ci) %>%
  select(-group_ci) %>%
  mutate(lower = as.numeric(sub("-.*", "", CI)),
         upper = as.numeric(sub(".*-", "", CI))) %>%
  select(-CI) %>% 
  arrange(stringr::str_extract(mutation_name, ".+?(?=-)"),
          # then order by the numeric in the amino acid sequence
          as.integer(stringr::str_extract(gsub("^.*-", "", mutation_name), "\\d+"))) %>%
  mutate(mutation_name = if_else(
    str_detect(mutation_name,"dhfr") == TRUE,
    paste("dhfr-",str_extract(mutation_name, "[^-]+$"),sep = ""),
    mutation_name)) %>% 
  mutate(mutation_name = toupper(mutation_name),#mutation_name = paste("Pf",mutation_name,sep = ""),
         group=ifelse(group=="Z_t","ZB-Travelers",ifelse(group=="Z_nt","ZB-NonTravelers","mTZ"))) %>% 
  mutate(mutation_name = toupper(mutation_name)) %>% 
  mutate(mutation_name = factor(mutation_name, levels = order)) %>%
  ggplot(aes(x = mutation_name, y = prevalence, fill = group)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, color = "black", position = position_dodge(0.7)) +
  scale_fill_manual(values = c("mTZ" = "#007BA5",
                               "ZB-Travelers" = "#F24000", 
                               "ZB-NonTravelers" = "#FFB203"),
                    name="Category") +
  ylab("Prevalence") +
  xlab("Mutation Name") +
  theme_minimal() +
  theme(text = element_text(size = 40),
        #axis.title.y = element_text(size = 40),
        axis.text.x = element_text(size = 25, angle = 45, hjust = 1),
        legend.position = "right")
ggsave(here("plots",run,"FigureS2.png"),
       width = 16, height = 10, dpi = 300,bg="white")

# Table S3 ----------------------------------------------------------------
border_style=officer::fp_border(color="black", width=1.5)
mutations_tbl <- c("k13-P441L","k13-R561H","k13-A675V")
k13_travel <- reformat %>% 
  filter(mutation_name %in% mutations_tbl) %>% 
  mutate(prev_Z_t = paste0(prevalence_Z_t," (",total_ALT_Z_t,"/",N_Z_t,", 95% CI: ",CI_Z_t,")"),
         prev_Z_nt = paste0(prevalence_Z_nt," (",total_ALT_Z_nt,"/",N_Z_nt,", 95% CI: ",CI_Z_nt,")"),
         prev_M=paste0(prevalence_M," (",total_ALT_M,"/",N_M,", 95% CI: ",CI_M,")")) %>% 
  select(c("mutation_name",contains("prev_"))) %>% 
  # order alphabetically
  arrange(stringr::str_extract(mutation_name, ".+?(?=-)"),
          # then order by the numeric in the amino acid sequence
          as.integer(stringr::str_extract(gsub("^.*-", "", mutation_name), "\\d+"))) %>%
  mutate(mutation_name = toupper(mutation_name)) %>% 
  #mutate(mutation_name = paste("Pf",mutation_name,sep = "")) %>% 
  flextable::flextable() %>%
  flextable::autofit() %>% 
  # Remove all existing borders
  flextable::border_remove() %>%
  # add  box theme
  flextable::theme_box() %>%
  flextable::bold(~ prev_Z_t > 0.001,2) %>% 
  flextable::bold(~ prev_Z_nt > 0.001,3) %>% 
  flextable::bold(~ prev_M > 0.001,4) %>% 
  flextable::set_header_labels(
    mutation_name = "Mutation",
    prev_Z_t = "ZB-Travelers",
    prev_Z_nt = "ZB-NonTravelers",
    prev_M="mTZ")  %>% 
  flextable::border(border = border_style) %>% 
  flextable::border(border = border_style, part = "header") %>%
  flextable::align(align = "left", j = c(2:4), part = "all") %>% 
  flextable::fontsize(size = 20,part="header") %>% 
  flextable::fontsize(size = 16, part = "all") %>% 
  flextable::width(j=c(1),width = 1.5) %>%
  flextable::width(j=c(2,3,4),width = 2) %>%
  flextable::height(i=c(2:3),height = 4) %>% 
  flextable::bg(bg = "#BACDEF", part = "header") %>% 
  flextable:: bg(i=1:3,j=1,bg = "#F1D7DA") 
flextable::save_as_image(k13_travel, path = here("tables",run,"TableS4.png"),webshot = "webshot2")
flextable::save_as_docx(k13_travel, path = here("tables",run,"TableS4.docx"))


# overall region ------------------------------------------------------------------
prev_table_ZAN <- comb_data_meta_rand_id %>%
  filter(region == "ZB") %>%
  group_by(mutation_name) %>%
  summarise(total_ALT_Z = sum(genotype >= 1, na.rm = TRUE),
            total_REF_Z = sum(genotype == 0, na.rm = TRUE),
            N_Z = total_ALT_Z + total_REF_Z,
            prevalence_Z = total_ALT_Z/N_Z) %>%
  mutate(
    CI = pmap(list(total_ALT_Z, N_Z), ~ binom.confint(x = ..1, n = ..2, methods = "exact")),
    ci.lower_Z = map_dbl(CI, 'lower'),
    ci.upper_Z = map_dbl(CI, 'upper')
  ) %>%
  select(-CI)

prev_main_zan_tbl <- left_join(prev_table_ZAN,
                               prev_table_TAN,
                               by = "mutation_name") %>%
  mutate_if(is.numeric, round,3) %>%
  unite(CI_Z,ci.lower_Z,ci.upper_Z,remove = T, sep = "-") %>%
  unite(CI_M,ci.lower_M,ci.upper_M,remove = T, sep = "-")

prev_table_sub <- prev_main_zan_tbl %>% 
  filter(mutation_name == "dhps-A581G")

# chi squared test for the in text statistic
tbl_t <- data.frame(present = c(8,
                                142),
                    absent = c(173,
                               1039))
chi_squared_result <- chisq.test(tbl_t)
