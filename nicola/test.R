source('lib.R')
library(dplyr)

data = readr::read_tsv('../msk_met_tropism_keap1_samples_maf.tsv') %>%
  mutate(Purity = Purity / 100) %>%
  filter(!is.na(Purity))

# Based on their analysis - Binomial test
# data$class = sapply(1:(data %>% nrow), function(i){
#
#   classification_data = mut_classifier(
#     coverage = data$t_ref_count[i] + data$t_alt_count[i],
#     p = data$Purity[i]/2,
#     alpha_level = alpha_level,
#     model = 'binomial'
#   )
#
#   if(data$t_alt_count[i] > classification_data$nv[2]) return("Clonal LOH")
#   if(data$t_alt_count[i] > classification_data$nv[1]) return("Clonal")
#   return("Subclonal")
# })
#
# sample = "P-0050568-T01-IM6"
#
# single_sample_plotting(data, sample)
#
# single_sample_plotting(data, data$Tumor_Sample_Barcode[324])

samples_list = data$Tumor_Sample_Barcode %>% unique()

data_classified = lapply(samples_list, function(sample){
  
  sample_data = data %>% filter(Tumor_Sample_Barcode == sample)
  
  cli::cli_h1(sample)
  
  bmix = purity_bmix(
    nv = sample_data$t_alt_count,
    coverage = sample_data$t_ref_count + sample_data$t_alt_count,
    purity = sample_data$Purity %>% unique(),
    model = "BetaBinomial",
    eps = 0.01
  )
  
  sample_data$purity_bmix = bmix$purity
  
  # Class per mutation
  sample_data$class = sapply(1:(sample_data %>% nrow), function(i){
    
    classification_data = mut_classifier(
      coverage = sample_data$t_ref_count[i] + sample_data$t_alt_count[i],
      p = sample_data$Purity[i]/2,
      rho = 0.01,
      alpha_level = 0.01,
      model = "BetaBinomial"
    )
    
    if(sample_data$t_alt_count[i] > classification_data$nv[2]) return("Clonal LOH")
    if(sample_data$t_alt_count[i] > classification_data$nv[1]) return("Clonal")
    return("Subclonal")
  })
  
  sample_data$class_bmix = sapply(1:(sample_data %>% nrow), function(i){
    
    if(is.na(sample_data$purity_bmix %>% unique())) return(NA)
    
    classification_data = mut_classifier(
      coverage = sample_data$t_ref_count[i] + sample_data$t_alt_count[i],
      p = sample_data$purity_bmix[i]/2,
      rho = 0.01,
      alpha_level = 0.01,
      model = "BetaBinomial"
    )
    
    if(sample_data$t_alt_count[i] > classification_data$nv[2]) return("Clonal LOH")
    if(sample_data$t_alt_count[i] > classification_data$nv[1]) return("Clonal")
    return("Subclonal")
  })
  
  return(sample_data)
  
}) %>% 
  do.call(rbind, .)

write_tsv(x = data_classified, file = "./data_classified.tsv")

data_classified %>% filter(Hugo_Symbol=="KEAP1") %>% pull(class) %>% table()

data_classified %>% filter(Hugo_Symbol=="KEAP1") %>% pull(class_bmix) %>% table()

plots_list = lapply(samples_list[1:50],
                    single_sample_classify,
                    data = data,
                    alpha_level = 0.01,
                    model = 'BetaBinomial',
                    rho = 0.01)

pdf("Fits_example.pdf", width = 6, height = 8)
plots_list %>% print()
dev.off()
