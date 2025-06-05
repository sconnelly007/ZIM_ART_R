here::i_am("scripts/analysis/zim_brief_report/TableS1.R")
library(here)
library(tidyverse)

run <- "zim_brief_report"
comb_data_meta_rand_id <- readRDS(here("data/processed",run,"dr_processed/combined_dr_data.rds"))

sample_size_c <- comb_data_meta_rand_id %>% 
  select(random_sample_id,clinic) %>% 
  distinct() %>% 
  group_by(clinic) %>% 
  summarise(n=n()) %>% 
  ungroup() %>% 
  left_join(comb_data_meta_rand_id %>% 
              select(clinic,name,clinic_district) %>%  
              distinct(),by="clinic") %>% 
  select(clinic,name,clinic_district,n) %>% 
  arrange(-n)

# create sample size table
border_style = officer::fp_border(color="black", width=3)

# create sample size by clinic
samp_tbl <- sample_size_c %>%  
  flextable::flextable() %>% 
  flextable::autofit() %>% 
  flextable::width(j=c(1,2,3), width = 2) %>% 
  flextable::set_header_labels(
    clinic = "Clinic",
    name = "Clinic Name",
    clinic_district = "District",
    n  = "Sample Size") %>% 
  flextable::border_remove() %>% 
  # add  box theme
  flextable::theme_box() %>%
  # change borders
  flextable::border(border = border_style) %>% 
  flextable::border(border = border_style, part = "header") %>% 
  flextable::align(align = "center",part = "all")
flextable::save_as_image(samp_tbl, path = here("tables",run,"TableS1.png"), webshot = "webshot2")
flextable::save_as_docx(samp_tbl, path = here("tables",run,"TableS1.docx"))
