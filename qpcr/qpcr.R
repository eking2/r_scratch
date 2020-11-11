library(tidyverse)
library(readxl)
library(cowplot)
library(ggsci)

# raw data contains wbc comparison and preservation runs
df <- read_xls('qpcr_1.xls', sheet='Results', skip=34) %>%
    select(Well, `Well Position`, `Sample Name`, `Target Name`, CT) %>%
    rename(well = Well, 
           position = `Well Position`, 
           sample = `Sample Name`,
           target = `Target Name`)

# compare signal of each gene with increasing white blood cells
wbc <- df %>%
    filter(well > 193) %>%
    filter(target != 'GYPA') %>%
    mutate(CT = recode(CT, 'Undetermined' = '40')) %>%
    mutate(CT = as.numeric(CT)) %>%
    group_by(sample, target) %>%
    summarize(mean_ct = mean(CT),
              sd_ct = sd(CT))

wbc$target <- factor(wbc$target, levels=c('HER2', 'PTPRC', 'GAPDH'))
wbc$sample <- factor(wbc$sample, levels=c('WBC 100', 'WBC 1K', 'WBC 10K', 'WBC 100K', 'WBC 640K', 'Washbuf + cells', 'Washbuf only', 'water blank'))

# colorblind 
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

p <- wbc %>%
    mutate(sample = recode(sample, 'Washbuf + cells' = 'Buffer + CTC',
                                   'Washbuf only' = 'Buffer',
                                   'water blank' = 'Water')) %>%
    ggplot(aes(x=sample, y=40 - mean_ct, fill=sample)) +
    geom_bar(stat='identity', width=0.6, color='black') +
    geom_errorbar(aes(ymin=40 - mean_ct - sd_ct, ymax=40 - mean_ct + sd_ct), width=0.08, size=0.5) +
    facet_grid(target~.) +
    theme(legend.position = 'none') +
    labs(title='Detection of Circulating Tumor Cells with qPCR', x='', caption='HER2 measures CTC, PTPRC measures WBC, GAPDH loading control\nWBC 100K pipetting error', y='Expression (40 - Ct)', subtitle='100 CTC Spike') +
    theme_half_open(14) +
    panel_border() +
    background_grid() + 
    theme(legend.position = 'none') +
    theme(strip.background = element_rect(fill='gray92')) + 
    theme(axis.text.x = element_text(angle=45, hjust=1)) +
    scale_fill_manual(values=cbPalette) 
   # scale_y_discrete(expand=c(0, 0))

write_csv(wbc, 'wbc_calc.csv')
ggsave('wbc.pdf', p, width=8)

# evaluate magnetic beads pre-treatment, measure effect on WBC background and target CTC signal
magbeads <- df %>%
    #filter(well < 193 | well %in% c('304', '306')) %>%
    filter(well < 193) %>%
    mutate(CT = recode(CT, 'Undetermined' = '40')) %>%
    mutate(CT = as.numeric(CT)) %>%
    mutate(timepoint = str_extract(sample, '\\d+')) %>%
    mutate(timepoint = replace_na(timepoint, 2)) %>%
    #mutate(timepoint = as.numeric(timepoint)) %>%
    mutate(name = ifelse(str_detect(sample, 'EDTA'), str_split_fixed(sample, ' ', n=2)[, 2], sample)) %>%
    mutate(name = recode(name, 'EDTA' = 'Baseline',
                               'EDTA-MB' = 'Magnetic Beads',
                               'EDTA no cells' = 'No CTC',
                               'cDNA blank' = 'cDNA Blank',
                               'lysis blank' = 'Lysis Blank')) %>%
    group_by(target, timepoint, name) %>%
    summarise(mean_ct = mean(CT), sd_ct = sd(CT))

magbeads$target <- factor(magbeads$target, levels=c('HER2', 'PTPRC', 'GYPA', 'GAPDH'))
magbeads$name <- factor(magbeads$name, levels=c('Baseline', 'Magnetic Beads', 'No CTC', 'cDNA Blank', 'Lysis Blank'))

# ignore gapdh 2hr errors
magbeads <- magbeads %>%
    mutate(sd_ct = replace(sd_ct, target == 'GAPDH' & timepoint == 2, 0))

p2 <- magbeads %>%
    ggplot(aes(x = name, y=40-mean_ct, fill=timepoint)) +
    geom_bar(stat='identity', position='dodge', color='black') +
    geom_errorbar(aes(ymin=40 - mean_ct - sd_ct, ymax=40 - mean_ct + sd_ct), width=0.2, size=0.5, position=position_dodge(0.9)) +
    facet_grid(.~target) +
    theme_half_open(14) +
    panel_border() +
    background_grid() +
    theme(strip.background = element_rect(fill='grey92')) +
    theme(legend.position = 'bottom') +
    labs(title='Magnetic Beads White Blood Cell Separation', x='', y='Expression (40 - Ct)', fill='Timepoint (hr)') +
    scale_fill_d3() +
    theme(axis.text.x = element_text(angle=45, hjust=1)) 

write_csv(magbeads, 'magbeads_calc.csv')
ggsave('magbeads.pdf', p2, width=12, height=8)
