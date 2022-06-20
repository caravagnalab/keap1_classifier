source('library.R')
library(dplyr)

# Read input table
data_file = '../msk_met_tropism_keap1_samples_maf.tsv'
data = readr::read_tsv(data_file)

# Convert purity format and filter NA values
data = data %>% 
  mutate(Purity = Purity / 100) %>%
  filter(!is.na(Purity))

# Set  parameters for BetaBinomial test
rho = 0.02    # overdisperion
alpha = 1e-5  # confidence level

# Get list of samples
samples_list = data$Tumor_Sample_Barcode %>% unique()

# Run test
data_classified = lapply(samples_list, function(sample){
  
  # Get mutations of specific sample 
  sample_data = data %>% filter(Tumor_Sample_Barcode == sample)
  cli::cli_h1(sample)
  
  # Classify each mutation
  sample_data$class = sapply(1:(sample_data %>% nrow), function(i){
    
    classification_data = mut_classifier(
      coverage = sample_data$t_ref_count[i] + sample_data$t_alt_count[i],
      p = sample_data$Purity[i]/2,
      rho = rho,
      alpha_level = alpha,
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

## Run test and make plots

plots_list = lapply(samples_list[1:50],
                    single_sample_classify,
                    data = data,
                    alpha_level = alpha,
                    model = 'BetaBinomial',
                    rho = rho)

pdf("./Fit_examples.pdf", width = 6, height = 8)
plots_list %>% print()
dev.off()




