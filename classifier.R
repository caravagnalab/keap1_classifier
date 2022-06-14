setwd("~/cds/regina_elena/")
library(tidyverse)
library(BMix)

# x = x %>% filter(!is.na(Purity))

### Functions
## function that fits VAF spectra with max 3 Binomials
bmix_fit = function(x){
  data = data.frame(
    successes = x$t_alt_count,
    trials = x$t_alt_count+x$t_ref_count)
  fit = bmixfit(data, K.Binomials = 1:3, K.BetaBinomials = 0)
  fit$data = data
  return(fit)
}
## function that estimates purity from binomial peaks fitted from VAF spectra
purity_estimate = function(x, purity, eps){
  if(is.null(x)) return(NA)
  nbin = x$K["B"]
  binomials = sort(x$B.params)
  p_est = NA
  if(nbin == 3) p_est = binomials[2]*2
  if(nbin == 2 & abs((binomials[2]*2)-purity) <= eps)  p_est = binomials[2]*2
  if(nbin == 2 & abs((binomials[2]*2)-purity) > eps & abs((binomials[1]*2)-purity) <= eps) p_est = binomials[1]*2
  if(nbin == 1) p_est = binomials[1]*2
  return(as.double(p_est))
}

## binomial test for subclonal classe
# test_clonality = function(x, sample, purity, alpha){
#   x = x %>% 
#     filter(Tumor_Sample_Barcode == sample, Hugo_Symbol=="KEAP1")
#   nsucc = x$t_alt_count
#   ntrials = x$t_ref_count+x$t_alt_count
#   prob = purity/2
#   output = binom.test(x = nsucc,
#                            n = ntrials,
#                            p = prob,
#                            alternative = "less",
#                            conf.level = 1-alpha)
#   x$p_subclonal = output$p.value
#   x$subclonal = (x$p_subclonal < alpha)
#   x$loh = FALSE
#   x$p_loh = 0
#   if(x$p_subclonal > alpha){
#     output = binom.test(x = nsucc,
#                         n = ntrials,
#                         p = prob,
#                         alternative = "greater",
#                         conf.level = 1-alpha)
#     x$p_loh = output$p.value
#     x$loh = (x$p_loh < alpha)
#   }
#   return(x %>% as_tibble())
# }
test_clonality = function(nv, dp, purity, alpha){
  
  if(is.na(purity)) return(NA)
  
  # Test support for subclonal mutation
  subclonal_pvalue = binom.test(
    x = nv,
    n = dp,
    p = purity / 2,
    alternative = "less",
    conf.level = 1 - alpha
  )$p.value
  
  if(subclonal_pvalue < alpha) return("SUBCLONAL")
  
  # Classify what type of clonal is
  loh_pvalue = binom.test(
    x = nv,
    n = dp,
    p = purity / 2,
    alternative = "greater",
    conf.level = 1 - alpha
  )$p.value
  
  if(loh_pvalue < alpha) return("CLONAL LOH")
  else return("CLONAL")
}


## Plot settings
colormap = c("red2", "magenta", "black", "grey")
names(colormap) = c("KEAP1 Clonal Heterozygous", "KEAP1 Clonal LOH", "KEAP1 Subclonal", "Other")

## Plot functions

plot_vaf_patient = function(x, patient){
  
  x = x %>% 
    filter(Tumor_Sample_Barcode == patient)
  
  purity_path = x %>% pull(Purity) %>% unique()
  purity_bmix = x %>% pull(purity_estimate) %>% unique() %>% as.double()
  
  colormap = c("red2", "magenta", "black", "grey")
  names(colormap) = c("Clonal Heterozygous", "Clonal LOH", "Subclonal", "Other")
  
  ggplot(x, aes(x = VAF, fill = class))+
    geom_histogram(binwidth = 0.01)+
    scale_fill_manual(values = colormap)+
    geom_vline(aes(xintercept = purity_bmix, linetype = "estimated"))+
    # geom_text(aes(x = purity_bmix, y = 0, vjust = -1, label = "estimated purity"), colour="blue")+
    geom_vline(aes(xintercept = purity_path, linetype = "measured"))+
    scale_linetype_manual(name = "purity", values = c(measured = "solid", estimated = "dotted"))+
    # geom_text(aes(x = purity_path, y = 0, label = "measured purity"), colour="orange")+
    xlim(0,1)+
    # ggtitle(paste0("Sample: ", patient, "\nSample purity: ", purity))+
    CNAqc:::my_ggplot_theme()
}

## Main
x =  read.table(file = "./msk_met_tropism_keap1_samples_maf.tsv", 
                header = TRUE, 
                sep = "\t") %>% 
  as_tibble()

samples = x %>% pull(Tumor_Sample_Barcode) %>%
  unique()

fits = lapply(samples[1:50], function(s){
  w = x %>%
    filter(Tumor_Sample_Barcode == s)
  if(nrow(w)>3){
    fit = bmix_fit(w)
  } else {fit = NULL}
  return(fit)
})


names(fits) = samples[1:50]

purity_in = x %>% 
  group_by(Tumor_Sample_Barcode) %>% 
  filter(row_number() == 1) %>% 
  select(Tumor_Sample_Barcode, Purity) %>% 
  ungroup() %>% 
  mutate(Purity = Purity/100)

purity_in$purity_bmix = sapply(samples, function(s) {
  if (!(s %in% names(fits)))
    return(NA)
  purity_estimate(fits[[s]],
                  purity_in %>%
                    filter(Tumor_Sample_Barcode == s) %>%
                    pull(Purity),
                  0.05)
})

purity_in$plot_bmix = lapply(samples, function(s) {
  
  if (!(s %in% names(fits)) | is.null(fits[[s]]))
    return(NA)
  
  x_1 = purity_in %>%
    filter(Tumor_Sample_Barcode == s) %>%
    pull(Purity)
  
  x_2 = purity_in %>%
    filter(Tumor_Sample_Barcode == s) %>%
    pull(purity_bmix)
  
  BMix::plot_clusters(fits[[s]], fits[[s]]$data) +
    geom_vline(aes(xintercept = x_1, linetype = "input")) +
    geom_vline(aes(xintercept = x_2, linetype = "estimated")) + 
    scale_linetype_manual(name = "purity", values = c(input = "solid", estimated = "dotted"))
})

pdf("bmix.pdf", width = 6, height = 6)
purity_in$plot_bmix %>% print
dev.off()

KEAP1_status =  x %>%
  filter(Hugo_Symbol == "KEAP1") %>% 
  select(-Purity) %>% 
  left_join(purity_in, by = 'Tumor_Sample_Barcode') %>% 
  rowwise() %>% 
  dplyr::mutate(
    classification_Purity = test_clonality(
      nv = t_alt_count,
      dp = t_ref_count + t_alt_count,
      purity = Purity,
      alpha = 0.001
    ),
    classification_purity_bmix = test_clonality(
      nv = t_alt_count,
      dp = t_ref_count + t_alt_count,
      purity = purity_bmix,
      alpha = 0.001
    )
  )

748*0.21*0.79

for(i in 1:nrow(KEAP1_status))
{
  cat(i, '\n')
  if(is.na(KEAP1_status$plot_bmix[[i]])) next
  
  ip = NA
  
  if(!is.na(KEAP1_status$classification_Purity[i]))
  ip = KEAP1_status$classification_Purity[i]
  
  bp = NA
  
  if(!is.na(KEAP1_status$classification_purity_bmix[i]))
    bp = KEAP1_status$classification_purity_bmix[i]
  
  
  KEAP1_status$plot_bmix[[i]] = 
  KEAP1_status$plot_bmix[[i]] +
    labs(
      title = paste("Input: ", ip),
      subtitle = paste("BMix: ", KEAP1_status$classification_purity_bmix[i])
      ) +
    geom_vline(xintercept = KEAP1_status$VAF[i], color = 'purple', size = 2, alpha = .6)
}
  
pdf("clasification.pdf", width = 4, height = 4)
KEAP1_status$plot_bmix %>% print
dev.off()







fits = lapply(samples, function(s){
                  w = x %>%
                    filter(Tumor_Sample_Barcode == s)
                  if(nrow(w)>3){
                    fit = bmix_fit(x %>%
                               filter(Tumor_Sample_Barcode == s))
                    pe = purity_estimate(fit)
                  } else {pe = NA}
                  w$purity_estimate = pe
                  w
                }) %>% 
  do.call(rbind, .)

alpha = 0.05

test = apply(fits %>% filter(!is.na(Purity)), 1, test_clonality) %>% 
  do.call(rbind, .)

test = test %>% 
  mutate(class = case_when(
    Hugo_Symbol == "KEAP1" & !loh & subclonal ~ "Subclonal",
    Hugo_Symbol == "KEAP1" & !subclonal & loh ~ "Clonal LOH",
    Hugo_Symbol == "KEAP1"  & !subclonal & !loh ~ "Clonal Heterozygous",
    TRUE ~ "Other"
  ))

test %>% pull(class) %>% table()

toplot = test %>% 
  group_by(Tumor_Sample_Barcode) %>% 
  filter(any(class %in% c("Clonal Heterozygous", "Clonal LOH", "Subclonal", "Other"))) %>% 
  mutate(Purity = as.double(Purity)/100,
         VAF = as.double(VAF))
  
plist = lapply(toplot %>%
                 pull(Tumor_Sample_Barcode) %>% 
                 unique(), function(s){plot_vaf_patient(toplot, s)+
                     theme(axis.title.x = element_blank(),
                           axis.title.y = element_blank())})

samples_w_purity = toplot %>% pull(Tumor_Sample_Barcode) %>% unique()

lapply(1:length(samples_w_purity), function(x){
  ggsave(
  filename = paste0("./classifier_VAF/",samples_w_purity[x],"_.png"),
  plot = plist[[x]],
  device = "png",
  dpi = "print",
  width = 8,
  height = 8)})


plt = ggpubr::ggarrange(plotlist = plist[1:20], 
                        ncol = 5, 
                        nrow = 4, 
                        common.legend = TRUE, 
                        legend = "bottom")

z = test %>% 
  mutate(Purity = as.double(Purity)/100,
         VAF = as.double(VAF),
         t_ref_count = as.integer(t_ref_count),
         t_alt_count = as.integer(t_alt_count))

plt_n = ggplot(z %>% filter(Hugo_Symbol=="KEAP1"))+
  geom_point(aes(x=VAF,y=(t_alt_count+t_ref_count), color=class))+
  xlim(0,1)+
  ylab("Sequencing depth")+
  scale_color_manual(values = colormap)+
  facet_wrap(~Purity)+
  CNAqc:::my_ggplot_theme()

ggsave(filename = "./dp_vaf.png",
      plot = plt_n,
      device = "png",
      dpi = "print",
      width = 8,
      height = 8)

dat_text <- data.frame(
  Tumor_Sample_Barcode = toplot %>% summarise() %>% pull(Tumor_Sample_Barcode),
  purity = toplot %>% summarise(n = mean(Purity)) %>% pull(n),
  label   = paste0("\u03C0","=" , toplot %>% summarise(n = mean(Purity)) %>% pull(n))
)

plist = lapply(toplot %>% pull(Purity) %>% unique(), function(p){
  ggplot(toplot %>% filter(Purity==p))+geom_histogram(aes(x=VAF, fill=class), binwidth = 0.01)+
  scale_fill_manual(values = colormap)+
  facet_wrap(~Tumor_Sample_Barcode, scales = "free_y")+
  CNAqc:::my_ggplot_theme()+
  theme(legend.text=element_text(size=14), legend.title = element_text(size=14))+
  theme(axis.title.x = element_text(size=14))+
  ggtitle(label = paste0("Purity ", p))
})

purities = toplot %>% pull(Purity) %>% unique()
lapply(purities, function(x){
  plt = plist[[which(purities==x)]]
  ggsave(filename = paste0("./classifier_",x,".png"), 
         plot = plt, 
         device = "png", 
         dpi = "print",
         width = 12, 
         height = 16)
})

plt_final = plt+geom_text(
    data    = dat_text,
    mapping = aes(x = -Inf, y = -Inf, label = label),
    hjust   = -2,
    vjust   = -5.5
  )

ggsave(filename = "./classifier.png", 
       plot = plt_final, 
       device = "png", 
       dpi = "print",
       width = 18, 
       height = 22)
toplot
x %>% pull(Tumor_Sample_Barcode) %>% unique() %>%  length()
