here::i_am("scripts/analysis/zim_brief_report/TableS5.R")
library(here)
library(tidyverse)

run <- "zim_brief_report"
comb_data_meta_rand_id <- readRDS(here("data/processed",run,"dr_processed/combined_dr_data.rds"))

district_prev <- comb_data_meta_rand_id %>%
  # filter to point mutations
  filter(!(str_detect(mutation_name,"dup|del|\\*|fs|TL555|TLISC555|LI556|LISC556|TLI555"))) %>% 
  filter(exonic_func == "missense_variant") %>% 
  mutate(mutation_name = gsub("PF3D7_1115700","FP2a",mutation_name),
         mutation_name = gsub("PF3D7_1224000","GCH1",mutation_name),
         mutation_name = gsub("PF3D7_1251200","coronin",mutation_name)) %>% 
  group_by(clinic_district,mutation_name) %>%
  summarise(total_ALT_t = sum(genotype >= 1, na.rm = TRUE),
            total_REF_t = sum(genotype == 0, na.rm = TRUE),
            N_t = total_ALT_t + total_REF_t,
            prevalence = total_ALT_t/N_t) %>%
  mutate(
    CI = pmap(list(total_ALT_t, N_t), ~ binom.confint(x = ..1, n = ..2, methods = "exact")),
    ci.lower_t = map_dbl(CI, 'lower'),
    ci.upper_t = map_dbl(CI, 'upper')
  ) %>% 
  left_join(comb_data_meta_rand_id %>% 
              select(clinic_district,region) %>%  
              distinct(),by="clinic_district") %>% 
  arrange(region, clinic_district) %>% 
  select(-region)

district_prev <- district_prev %>%
  select(-c("CI","total_REF_t")) %>% 
  mutate_if(is.numeric, round,3) %>%
  mutate(ci.lower_t = paste0(ci.lower_t*100,"%"),
         ci.upper_t = paste0(ci.upper_t*100,"%"),
         prevalence = paste0(prevalence*100,"%")) %>% 
  unite(fraction,total_ALT_t,N_t,remove = T, sep = "/") %>% 
  unite(CI,ci.lower_t,ci.upper_t,remove = T, sep = "-") %>%
  pivot_longer(cols = c(fraction, prevalence, CI), 
               names_to = "metric", 
               values_to = "value") %>%
  pivot_wider(names_from = clinic_district, values_from = value)

writexl::write_xlsx(district_prev,
                    here("tables",run,"TableS5.xlsx"))