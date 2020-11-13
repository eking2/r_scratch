library(tidyverse)
library(broom)
library(cowplot)

df <- read_csv('bradford.csv')

# split samples and standards
standards <- df %>%
    filter(Dil == 'std')

samples <- df %>%
    filter(Dil != 'std') %>%
    mutate(Dil = as.numeric(Dil))

# standard curve with polynomial 2
# check fit r squared
curve <- standards %>%
    group_by(Dil) %>%
    nest() %>%
    mutate(model = map(data, ~lm(Sample ~ I(abs_595^2) + abs_595, data=.x)),
           tidied = map(model, tidy),
           glanced = map(model, glance)) %>%
    mutate(r2 = map_dbl(glanced, 'adj.r.squared')) %>%
    unnest(tidied) %>%
    select(Dil, r2, term, estimate) %>%
    pivot_wider(id_cols=c(Dil, r2), names_from=term, values_from=estimate) %>%
    rename(sample = 1, r2=2, intercept=3, `x^2`=4, x=5)

# fit samples to standard curve
# mult by dilution factor to get final concentration 
samples <- samples %>%
    mutate(conc_ug_ml = (curve$`x^2` * abs_595^2 + curve$x * abs_595 + curve$intercept) * Dil) %>%
    mutate(conc_mg_ml = conc_ug_ml / 1000)

samples_agg <- samples %>%
    group_by(Sample) %>%
    summarize(mean_conc_mg_ml = mean(conc_mg_ml), sd_conc_mg_ml = sd(conc_mg_ml)) %>%
    mutate(perr = 100 * sd_conc_mg_ml / mean_conc_mg_ml) %>%
    mutate(w_glycerol = mean_conc_mg_ml * 0.6) %>%
    mutate_at(2:5, round, 3)

# plot expression
p <- samples_agg %>%
    ggplot(aes(x=Sample, y=mean_conc_mg_ml, fill=Sample)) + 
    geom_bar(stat='identity', color='black', width=0.7) +
    geom_errorbar(aes(ymin=mean_conc_mg_ml - sd_conc_mg_ml, ymax=mean_conc_mg_ml + sd_conc_mg_ml), width=0.2) +
    geom_jitter(data=samples, aes(y=conc_mg_ml), width=0.3, shape=21, color='black', fill='red', size=2.5, alpha=0.7) +
    theme_half_open(14) +
    labs(title='Purified Protein Concentrations', subtitle='Calculated from Bradford Assay', x='', y='Concentration (mg/mL)') +
    background_grid() +
    theme(legend.position = 'none') +
    scale_fill_viridis_d()

standards$Sample <- factor(standards$Sample, levels=c(0, 25, 125, 250, 500, 750, 1000, 1500))
rsq_label <- paste('R^2 == ', round(curve$r2, 3))
form_label <- paste0('y = ', round(curve$`x^2`, 1), 'x^2 + ', round(curve$x, 1), 'x + ',  round(curve$intercept, 1))

# plot standard curve
p2 <- curve %>%
    ggplot() + 
    geom_function(fun = function(x) (curve$`x^2` * x^2 + curve$x * x + curve$intercept), color='orange', size=1.5) +
    geom_point(data=standards, aes(x=abs_595, y=as.numeric(as.character(Sample)), fill=Sample), shape=22, size=3) +
    theme_half_open(14) +
    labs(title='Bradford Assay Standard Curve', x=expression('OD'[595]), y='Concentration (ug/mL)', fill='Standard (ug/mL)') +
    background_grid() +
    theme(legend.position = 'bottom') +
    geom_text(aes(x=0.3, y=1400), label=form_label, size=6, hjust=0) +
    geom_text(aes(x=0.3, y=1300), label=rsq_label, size=6, hjust=0, parse=TRUE) 

#) save results
write_csv(samples, 'sample_concs.csv')
ggsave('bradford.pdf', p, width=9, height=6)
ggsave('std_curve.pdf', p2, width=9, height=7)

