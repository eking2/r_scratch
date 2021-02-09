library(tidyverse)
library(cowplot)
library(ggrepel)

# specific activities
calc_activity <- function(slope, blank, dil, conc, enz_vol=10, rxn_vol=100, ext=6.22) {
    
    # ext input as mM^-1 cm^-1
    # vols are uL
    # conc is mg/ml
    
    # enz in reaction, mg
    enz_rxn <- conc * 1/dil * enz_vol * 1/1000
    
    # mA/min
    d_slope <- slope - blank
    
    # mA/min to mM/min, use A = eps * path * conc
    path <- rxn_vol/400
    res <- d_slope * 1/ext * 1/path * 1/1000 # mM/min
    
    # molar to mol
    res <- res * rxn_vol/10^6  # mmol/min
    res <- res/enz_rxn * 10^6 # nmol/min-mg
    
    return (round(res, 3))
}


div_prop <- function(f, sigma_a, a, sigma_b, b) {
    # no covariance
    sigma_f <- abs(f) * ( (sigma_a / a)^2 + (sigma_b / b)^2 )^0.5
    return (sigma_f)
}


# read in data, make tidy
concs <- read_csv('conc.csv')
slopes <- read_csv('slopes.csv') %>%
    pivot_longer(cols=c(`1`, `2`), names_to='rep', values_to='slope') %>%
    drop_na() %>%
    inner_join(concs)

blanks <- slopes %>%
    filter(cofa_name == 'no') %>%
    select(c('sample_name', 'slope')) %>%
    rename('blank' = 'slope')    

slopes_clean <- slopes %>%
    inner_join(blanks) %>%
    filter(cofa_name != 'no') %>%
    mutate(activity = calc_activity(slope, blank, DF, conc))

write_csv(slopes_clean, 'slopes_clean.csv')

# get specific activity ratios
# propogate errors
slopes_agg <- slopes_clean %>%
    group_by(sample_name, cofa_name) %>%
    summarize(mean_sa = mean(activity), sd_sa = sd(activity)) %>%
    pivot_wider(id_cols=sample_name, names_from=cofa_name, names_glue='{cofa_name}_{.value}', values_from=c(mean_sa, sd_sa)) %>%
    mutate(`NMN/NAD` = NMN_mean_sa / NAD_mean_sa) %>%
    mutate(`NMN/NADP` = NMN_mean_sa / NADP_mean_sa) %>%
    mutate(`NMN/NAD_sd` = div_prop(`NMN/NAD`, NMN_sd_sa, NMN_mean_sa, NAD_sd_sa, NAD_mean_sa)) %>%
    mutate(`NMN/NADP_sd` = div_prop(`NMN/NADP`, NMN_sd_sa, NMN_mean_sa, NADP_sd_sa, NADP_mean_sa)) 
    
p <- slopes_agg %>%
    ggplot(aes(x=`NMN/NADP`, y=`NMN/NAD`, color=NMN_mean_sa)) +
    geom_point(aes(size = NMN_mean_sa)) +
    geom_errorbarh(aes(xmax = `NMN/NADP` + `NMN/NADP_sd`, xmin=`NMN/NADP` - `NMN/NADP_sd`), color='black', height=0.08) +
    geom_errorbar(aes(ymax= `NMN/NAD` + `NMN/NAD_sd`, ymin = `NMN/NAD` - `NMN/NAD_sd`), color='black', width=0.06) +
    labs(x='NMN/NADP', y='NMN/NAD', color='NMN Specific Activity\n(nmol/min/mg)', title='Cofactor Activity Ratios',
         subtitle='Comparing orthogonality at 37°C, DL-G3P') +
    theme_half_open(12) +
    background_grid() +
    scale_size(guide='none') +
    scale_color_viridis_c() +
    geom_text_repel(aes(label=sample_name), color='black', size=2.5, box.padding=0.45, min.segment.length = 0.6) +
    theme(legend.position = 'bottom', legend.key.width=unit(25, 'pt'))

save_plot('cofa_activity.pdf', p, base_height = 5, base_width = 6)

