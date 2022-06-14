source('lib.R')

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

plots_list = lapply(samples_list[1:50],
                    single_sample_classify,
                    data = data,
                    alpha_level = 0.01,
                    model = 'BetaBinomial',
                    rho = 0.01)

pdf("Fits_example.pdf", width = 6, height = 8)
plots_list %>% print()
dev.off()

