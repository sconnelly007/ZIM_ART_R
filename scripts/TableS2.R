library(readr)
library(dplyr)
library(flextable)

run <- "zim_brief_report"

# Read in the TSV file
mutation_df <- read_tsv(here("data/metadata",run,"targets.tsv"), show_col_types = FALSE)

# 
mutation_tbl <- mutation_df %>%
  select(
    Chromosome, Position, Gene_Name, AA_Change_Position,
    Aminoacid_Change, Mutation_Name_Single,
    Reference_Resistant, Gene_ID
  ) %>%
  rename(
    `Chr` = Chromosome,
    `Pos` = Position,
    `Gene` = Gene_Name,
    `AA Pos` = AA_Change_Position,
    `AA Change` = Aminoacid_Change,
    `Mutation` = Mutation_Name_Single,
    `Ref_Resistant?` = Reference_Resistant,
    `Gene ID` = Gene_ID
  )

# Reformat
border_style = officer::fp_border(color="black", width=3)
mutation_ft <- flextable(mutation_tbl) %>%
  set_header_labels(
    `Chr` = "Chromosome",
    `Pos` = "Position",
    `Gene` = "Gene",
    `AA Pos` = "AA Position",
    `AA Change` = "Amino Acid Change",
    `Mutation` = "Mutation",
    `Resistance?` = "Ref. Resistant",
    `Gene ID` = "Gene ID"
  ) %>%
  autofit() %>%
  flextable::border_remove() %>% 
  # add  box theme
  flextable::theme_box() %>%
  # change borders
  flextable::border(border = border_style) %>% 
  flextable::border(border = border_style, part = "header") %>% 
  flextable::align(align = "center",part = "all")

# Save flextable
flextable::save_as_image(mutation_ft, path = here("tables",run,"TableS2.png"), webshot = "webshot2")
flextable::save_as_docx(mutation_ft, path = here("tables",run,"TableS2.docx"))


