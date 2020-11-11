library(tidyverse)
library(broom)
library(cowplot)
library(ggsci)
library(stringr)

std <- read_csv('std_areas.csv')
samples <- read_csv('sample_areas.csv')
condition_map <- read_csv('conditions_map.csv')

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


# function to fit response ratio to standard curve
fit_response_ratio <- function(response_ratio, standard) {
    
    slope <- filter(lr, full_compound == standard)$slope
    intercept <- filter(lr, full_compound == standard)$intercept
    
    return (response_ratio*slope + intercept)
}

# calculate concentrations for sample
# drop sample 26 due to pipetting error
# set standard to use for each compound
# limit of detection with response ratio > 0.01
sample_ratios <- tibble(
     sample_num = samples$sample_num,
     condition = samples$condition,
     timepoint = samples$timepoint,
    `R-AC` = samples$`R-AC_area` / samples$KIP_area,
    `S-AC` = samples$`S-AC_area` / samples$KIP_area,
    `SS-BDO` = samples$`SS-BDO_area` / samples$KIP_area,
    `RR-BDO` = samples$`RR-BDO_area` / samples$KIP_area,
    `M-BDO` = samples$`M-BDO_area` / samples$KIP_area) %>%
    pivot_longer(-c(sample_num, condition, timepoint), names_to='compound', values_to='response_ratio') %>%
    filter(sample_num != 26) %>%
    mutate(standard = case_when(
        str_detect(compound, 'AC') ~ 'Acetoin',
        compound == 'SS-BDO' ~ '(S,S)-Butanediol',
        compound == 'RR-BDO' ~ '(R,R)-Butanediol',
        compound == 'M-BDO' ~ 'Meso-Butanediol'
    )) %>%
    mutate(concentration = ifelse(response_ratio > 0.01, map2_dbl(response_ratio, standard, ~fit_response_ratio(.x, .y)),
                                  0))

# pivot to wide for concentrations and aggregate over conditions
concs <- sample_ratios %>%
    #select(sample_num, condition, timepoint, compound, concentration)
    pivot_wider(id_cols=c(sample_num, condition, timepoint), names_from=compound, values_from=concentration) %>%
    mutate(Total = `R-AC` + `S-AC` + `SS-BDO` + `RR-BDO` + `M-BDO`) %>%
    select(-sample_num) %>%
    group_by(condition, timepoint) %>%
    summarize(
        across(everything(), list(mean=mean, sd=sd)))

# setup figure titles
make_title <- function(condition) {
    
    row <- filter(condition_map, code == condition)
    title = sprintf('%s (%d uM)\n%s (%d uM)\n%s (%d uM)\n%s (%d uM)\n%s', 
                    row$enzyme_1, row$enzyme_1_conc, row$enzyme_2, row$enzyme_2_conc,
                    row$enzyme_3, row$enzyme_3_conc, row$enzyme_4, row$enzyme_4_conc,
                    row$cofactors)
    return (title)
}


# from wide back to long
concs_piv <- concs %>%
    pivot_longer(-c(condition, timepoint), values_to='concentration') %>%
    separate(name, sep='_', into=c('compound', 'type')) %>%
    pivot_wider(id_cols=c(condition, timepoint, compound), names_from=type, values_from=concentration)


copy_zeros <- function(){
    
    # set 0hr timepoint for each condition
    zero <- concs_piv %>%
        filter(condition == 'i')
    
    res = NULL
    conditions = as.character(concs_piv$condition) 
    for (c in unique(conditions)) {
        if (c != 'i') {
            subset <- zero %>%
                mutate(condition = recode(condition, i = c))
            res <- bind_rows(res, subset)
        }
    }
    return (res)
}

# copy zero data from no enzyme to all conditions, and drop no enzyme single read
zeros <- copy_zeros()
concs_piv <- bind_rows(concs_piv, zeros) %>%
    filter(condition != 'i')

# set facet order
#concs_piv$condition <- factor(concs_piv$condition, levels = c('h', 'g', 'b', 'a', 'c', 'd', 'f', 'e'))

# plot traces concentration vs time
p2 <- concs_piv %>%
    ggplot(aes(x = timepoint, y=mean, color=compound)) + 
    geom_point() +
    geom_line() +
    geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=0.7, color='black') +
    facet_wrap(~condition, labeller=labeller(condition=make_title)) +
    theme_half_open(15) +
    panel_border() +
    background_grid() +
    labs(title='Meso-Butanediol GC Feeding Experiment', x='Timepoint (hr)', y='Concentration (g/L)', color='Compound') +
    theme(strip.background = element_rect(fill='grey92')) +
    scale_color_npg()

ggsave('conc_trace.pdf', p2, width=12, height=13)

# save outputs
write_csv(lr, 'standard_curve_params.csv')
write_csv(sample_ratios, 'sample_calc_ratios.csv')
write_csv(concs, 'sample_calc_concs.csv')
