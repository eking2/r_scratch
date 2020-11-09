library(tidyverse)
library(broom)
library(cowplot)
library(ggsci)

std <- read_csv('std_areas.csv')

# calculate response ratio to KIP
# remove standards below limit of detection
# map full compound names
ratios <- tibble(sample = std$sample,
                 conc = std$conc,
                 `AC` = (std$`S-AC_area` + std$`R-AC_area`) / std$KIP_area,
                 `M-BDO` = std$`M-BDO_area` / std$KIP_area,
                 `SS-BDO` = std$`SS-BDO_area` / std$KIP_area,
                 `RR-BDO` = std$`RR-BDO_area` / std$KIP_area) %>%
    pivot_longer(-c(sample, conc), names_to='compound', values_to ='resp_ratio') %>%
    filter(!sample %in% c('STD 1', 'STD 2')) %>%
    mutate(full_compound = case_when(
        compound == 'AC' ~ 'Acetoin',
        compound == 'SS-BDO' ~ '(S,S)-Butanediol',
        compound == 'RR-BDO' ~ '(R,R)-Butanediol',
        compound == 'M-BDO' ~ 'Meso-Butanediol'
    ))

# linear regression standard curve for each compound
lr <- ratios %>%
    group_by(full_compound) %>%
    nest() %>%
    mutate(model = map(data, ~lm(conc ~ resp_ratio, data=.x)),
           tidied = map(model, tidy)) %>%  # tidy gets params, glance get fit statistics, augment gets pred
    unnest(tidied) %>%
    select(full_compound, term, estimate) %>%
    pivot_wider(names_from=term, values_from=estimate) %>%
    rename(intercept = `(Intercept)`,
           slope = resp_ratio)

# set order of facets
ratios$full_compound = factor(ratios$full_compound, levels = c('Acetoin', 'Meso-Butanediol', '(S,S)-Butanediol', '(R,R)-Butanediol'))
lr$full_compound = factor(lr$full_compound, levels = c('Acetoin', 'Meso-Butanediol', '(S,S)-Butanediol', '(R,R)-Butanediol'))

# plot standard curves
p <- ratios %>%
    ggplot(aes(x = resp_ratio, y=conc, fill=full_compound)) +
    geom_point(size=4, pch=21, color='black') +
    geom_line(linetype='dashed', aes(color=full_compound)) +
    geom_abline(data=lr, aes(intercept = intercept, slope=slope), size=0.7) +
    facet_wrap(~full_compound) +
    theme_half_open(15) +
    panel_border() +
    background_grid() +
    theme(legend.position = 'none') +
    theme(strip.background = element_rect(fill='gray92')) + 
    labs(title='GC Standard Curves', subtitle = 'Areas relative to KIP', x='Response Ratio', y = 'Concentration (g/L)') +
    scale_fill_locuszoom() +
    scale_color_locuszoom() +
    geom_text(data=lr, aes(x=1.5, y=8.5, label=sprintf('Slope = %0.3f\nIntercept = %0.3f', slope, intercept)), hjust=0, size=4)

ggsave('gc_std_curve.pdf', p, width=8)
