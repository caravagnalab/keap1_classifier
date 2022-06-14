mut_classifier = function(coverage = 748,
                          p = .41 / 2,
                          rho = 0.01,
                          alpha_level = 1e-3,
                          model = 'binomial')
{
  # NVs range from
  nvs = 1:coverage

  log_p = NULL
  if (model == 'binomial')
  {
    # P(X <= nv) for X~Bin(.|coverage, p)
    log_p = sapply(nvs, dbinom, size = coverage, prob = p)
  }
  else
  {
    # P(X <= nv) for X~Bin(.|coverage, p)
    log_p = sapply(
      nvs,
      VGAM::dbetabinom,
      size = coverage,
      prob = p,
      rho = rho
    )
  }

  p_x = cumsum(log_p)


  # Expectation of X
  e_p = coverage * p

  # n: P(X <= n) < alpha
  # N: P(X > N) < 1 - alpha
  l_a = which(p_x < alpha_level, arr.ind = TRUE) %>% max
  r_a = which(p_x > 1 - alpha_level, arr.ind = TRUE) %>% min

  # Cutoffs in VAF space
  vafs = nvs / coverage
  l_v = vafs[l_a]
  r_v = vafs[r_a]

  inputs = data.frame(nv = nvs,
                      p = p_x,
                      VAF = vafs)

  return(list(
    model = model,
    density = inputs,
    nv = c(l_a, r_a),
    vaf = c(l_v, r_v)
  ))
}

test_power_plot = function(coverage = 748,
                           p = .41 / 2,
                           alpha_level = 1e-3,
                           rho = 0.01,
                           model = 'binomial')
{
  data_plot = mut_classifier(
    coverage = coverage,
    p = p,
    alpha_level = alpha_level,
    model = model,
    rho = rho
  )

  inputs = data_plot$density
  nvs = data_plot$density$nv %>% as.vector()
  l_a = data_plot$nv[1]
  r_a = data_plot$nv[2]
  l_v = data_plot$vaf[1]
  r_v = data_plot$vaf[2]

  model_string = "Binomial"
  if(data_plot$model == 'betabinomial')
    model_string = paste0("Beta-Binomial (rho = ", rho, ')')

  col_loh = 'steelblue'
  col_subclonal = 'indianred3'
  col_clonal = 'forestgreen'

  ggplot() +
    CNAqc:::my_ggplot_theme() +
    geom_rect(
      data = data.frame(
        xmin = 0,
        xmax = nvs[l_a],
        ymin = 0,
        ymax = 1
      ),
      aes(
        xmin = xmin,
        xmax = xmax,
        ymin = ymin,
        ymax = ymax
      ),
      fill = col_subclonal,
      alpha = .2
    ) +
    geom_rect(
      data = data.frame(
        xmax = coverage,
        xmin = nvs[r_a],
        ymin = 0,
        ymax = 1
      ),
      aes(
        xmin = xmin,
        xmax = xmax,
        ymin = ymin,
        ymax = ymax
      ),
      fill = col_loh,
      alpha = .2
    ) +
    geom_rect(
      data = data.frame(
        xmax = nvs[l_a],
        xmin = nvs[r_a],
        ymin = 0,
        ymax = 1
      ),
      aes(
        xmin = xmin,
        xmax = xmax,
        ymin = ymin,
        ymax = ymax
      ),
      fill = col_clonal,
      alpha = .2
    ) +
    # geom_vline(xintercept = e_p,
    #            col = 'black',
    #            linetype = 6) +
    geom_point(
      data = inputs %>% filter(nv < nvs[l_a]),
      aes(x = nv, y = p),
      size = .6,
      color = col_subclonal
    ) +
    geom_point(
      data = inputs %>% filter(nv > nvs[r_a]),
      aes(x = nv, y = p),
      size = .6,
      color = col_loh
    ) +
    geom_point(
      data = inputs %>% filter(nv >= nvs[l_a], nv <= nvs[r_a]),
      aes(x = nv, y = p),
      size = 1,
      color = col_clonal
    ) +
    geom_point(
      data = inputs,
      aes(x = nv, y = VAF),
      size = .6,
      shape = 3,
      color = 'gray'
    ) +
    geom_point(
      data = inputs %>% filter(nv == nvs[l_a]),
      aes(x = nv, y = VAF),
      size = 3,
      color = col_subclonal
    ) +
    geom_segment(
      data = inputs %>% filter(nv == nvs[l_a]),
      aes(
        x = nvs[l_a],
        y = VAF,
        xend = coverage,
        yend = VAF
      ),
      linetype = 'dashed',
      color = col_subclonal
    ) +
    geom_segment(
      data = inputs %>% filter(nv == nvs[l_a]),
      aes(
        x = nvs[l_a],
        y = VAF,
        xend = nvs[l_a],
        yend = 0
      ),
      linetype = 'dashed',
      color = col_subclonal
    ) +
    geom_point(
      data = inputs %>% filter(nv == nvs[r_a]),
      aes(x = nv, y = VAF),
      size = 3,
      color = col_loh
    ) +
    geom_segment(
      data = inputs %>% filter(nv == nvs[r_a]),
      aes(
        x = nvs[r_a],
        y = VAF,
        xend = coverage,
        yend = VAF
      ),
      linetype = 'dashed',
      color = col_loh
    ) +
    geom_segment(
      data = inputs %>% filter(nv == nvs[r_a]),
      aes(
        x = nvs[r_a],
        y = VAF,
        xend = nvs[r_a],
        yend = 1
      ),
      linetype = 'dashed',
      color = col_loh
    ) +
    geom_text(
      x = nvs[l_a] - 5,
      y = 1,
      label = 'subclonal',
      hjust = 1,
      size = 3,
      color = col_subclonal
    ) +
    geom_text(
      x = nvs[r_a] + 5,
      y = 0,
      label = 'clonal LOH',
      hjust = 0,
      size = 3,
      color = col_loh
    ) +
    labs(
      x = 'NV',
      y = "1 - P(X > NV)",
      caption = model_string,
      title = paste0("Coverage ", coverage, ' with purity ', p),
      subtitle = paste0('Alpha-level: ', alpha_level)
    ) +
    coord_cartesian(clip = 'off') +
    scale_y_continuous(sec.axis = sec_axis(~ ., name = 'VAF | coverage')) +
    theme(axis.line.y.right = element_line(color = "gray")) +
    annotation_custom(
      grid::textGrob(paste('>', round(r_v, 2)), gp = grid::gpar(col = col_loh, fontsize = 8)),
      xmin = coverage,
      xmax = coverage,
      ymin = r_v + .03,
      ymax = r_v + .03
    ) +
    annotation_custom(
      grid::textGrob(
        paste('<', round(l_v, 2)),
        gp = grid::gpar(col = col_subclonal, fontsize = 8)
      ),
      xmin = coverage,
      xmax = coverage,
      ymin = l_v - .03,
      ymax = l_v - .03
    ) +
    annotate(
      "text",
      x = 0,
      y = .95,
      label = "Subclonal",
      hjust = 0,
      size = 3,
      color = col_subclonal
    ) +
    annotate(
      "text",
      x = coverage,
      y = .95,
      label = "Clonal LOH",
      hjust = 1,
      size = 3,
      color = col_loh
    ) +
    annotate(
      "text",
      x = coverage * p,
      y = .95,
      label = "Clonal",
      hjust = .5,
      size = 3,
      color = col_clonal
    )
}

example_cartoons = function()
{
  ggpubr::ggarrange(
    test_power_plot(140, .3),
    test_power_plot(140, .3, model = 'betabinomial', rho = 0.001),
    test_power_plot(140, .3, model = 'betabinomial', rho = 0.01),
    ncol = 1,
    nrow = 3
  )

  ggpubr::ggarrange(
    test_power_plot(700, .4),
    test_power_plot(700, .4, model = 'betabinomial', rho = 0.001),
    test_power_plot(700, .4, model = 'betabinomial', rho = 0.01),
    ncol = 1,
    nrow = 3
  )
}

# Per-sample classifier
single_sample_classify = function(data, sample, alpha_level, model, rho)
{
  cli::cli_h1(sample)

  sample_data = data %>% filter(Tumor_Sample_Barcode == sample)

  # Class per mutation
  sample_data$class = sapply(1:(sample_data %>% nrow), function(i){

    classification_data = mut_classifier(
      coverage = sample_data$t_ref_count[i] + sample_data$t_alt_count[i],
      p = sample_data$Purity[i]/2,
      alpha_level = alpha_level,
      model = model,
      rho = rho
    )

    if(sample_data$t_alt_count[i] > classification_data$nv[2]) return("Clonal LOH")
    if(sample_data$t_alt_count[i] > classification_data$nv[1]) return("Clonal")
    return("Subclonal")
  })

  # KEAP1 stuff
  keap_1_status = sample_data %>%
    filter(Hugo_Symbol == 'KEAP1') %>%
    mutate(class = paste0(class, ' (', t_alt_count, '/', t_ref_count + t_alt_count, ')')) %>%
    pull(class) %>%
    paste(collapse = ', ')

  keap_1_vafs = sample_data %>%
    filter(Hugo_Symbol == 'KEAP1') %>%
    pull(VAF)

  # Fit plot
  fit_plot = sample_data %>%
    ggplot() +
    CNAqc:::my_ggplot_theme() +
    geom_vline(xintercept = keap_1_vafs, color = 'indianred3', linetype = 6) +
    geom_vline(xintercept = sample_data$Purity[1], color = 'black', linetype = 'dashed', size = .6) +
    geom_vline(xintercept = sample_data$Purity[1]/2, color = 'gray', linetype = 'dashed') +
    geom_histogram(binwidth = 0.01, aes(VAF, fill = class)) +
    ggsci::scale_fill_jama() +
    xlim(0, 1) +
    labs(title = sample,
         subtitle = paste0("Purity: ", sample_data$Purity[1], " - KEAP1: ", keap_1_status)) +
    guides(fill = guide_legend(''), color = guide_legend('', override.aes = aes(fill = NA)))

  # Test power for KEAP 1
  keap1_data = sample_data %>% filter(Hugo_Symbol == 'KEAP1')

  keap_1_plots = lapply(1:(keap1_data %>% nrow), function(i){
    test_power_plot(
      coverage = keap1_data$t_ref_count[i] + keap1_data$t_alt_count[i],
      p = keap1_data$Purity[i]/2,
      alpha_level = alpha_level,
      model = model,
      rho = rho
    ) +
      geom_vline(xintercept = keap1_data$t_alt_count[i], linetype = 'dashed', size = .5)
  })

  # Fig assembly
  lp = append(list(fit_plot), keap_1_plots)

  figure = ggpubr::ggarrange(
    plotlist = lp,
    nrow = lp %>% length,
    ncol = 1
  )

  return(figure)
}
