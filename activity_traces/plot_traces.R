library(tidyverse)
library(cowplot)
library(XML)
library(fs)
library(readxl)
library(ggsci)

xmls <- dir_ls(glob = '*.xml')

parse_xml <- function(xml_path) {
    
    # get kinetic absorbance data from each well
    
    xml <- xmlParse(xml_path)
    filename <- str_split(xml_path, '\\.')[[1]][1]
    
    abs <- xpathApply(xml, '//Well/RawData', xmlValue)
    time <- xpathApply(xml, '//Well/TimeData', xmlValue)[[1]][1]
    names <- xpathApply(xml, '//Well')
    names <- unlist(map(names, xmlGetAttr, 'Name'))
    
    abs <- str_split(abs, ' ')
    time <- str_split(time, ' ')
    
    df <- tibble(time = unlist(time))
    
    df_abs <- map(abs, tibble) %>%
        bind_cols() %>%
        set_names(names) %>%
        bind_cols(df) %>%
        pivot_longer(!time, values_to='abs', names_to='well') %>%
        mutate(time = as.numeric(time),
               abs = as.numeric(abs),
               trial = filename) 
    
    return (df_abs)
    
}

# combine all trials
df <- xmls %>%
    map(parse_xml) %>%
    bind_rows()


# save selected 
selected <- read_excel('labels.xlsx')
df_selected <- df %>%
    right_join(selected, by=c('trial', 'well')) %>%
    mutate(sample = factor(sample, levels = c('P253', 'A_graveolens', 'B_cereus', 'C_acetobutylicum', 'C_pasteurianum', 'N_meningitidis', 'P_sativum', 'S_agalactiae', 'S_pyogenes', 'S_solfataricus'))) %>%
    mutate(cofactor = recode(cofactor, 'no' = 'Water')) 
    

# plot dl-g3p 60c
dl_60 <- df_selected %>%
    filter(substrate == 'DL-G3P') %>%
    filter(temp == 60) %>%
    ggplot(aes(x = time, y=abs, color=cofactor, group=trial)) +
    geom_line(size=0.8) +
    facet_wrap(~sample) +
    coord_cartesian(ylim = c(0, 0.6)) +
    theme_half_open(12) +
    background_grid(color.major='gray94') +
    panel_border() +
    labs(title='gapN Traces', x='Time (s)', y=expression('OD'[340]), subtitle='2mM DL-G3P, 5mM BME, 60°C', color='Cofactor') +
    theme(strip.background = element_rect(fill='gray91')) +
    theme(legend.position = c(0.55, 0.18)) +
    theme(panel.spacing = unit(0.7, 'lines')) + 
    scale_color_manual(values = c('#0C5DA5', '#00B945', '#FF9500'))

save_plot('gapn_nad_traces.pdf', dl_60, base_width=9, base_height=7)

# nmn runs
# agg
nmn <- df_selected %>%
    filter(cofactor == 'NMN') %>%
    mutate(temp = factor(temp, c('25', '60'))) %>%
    group_by(time, sample, substrate, temp) %>%
    summarize(abs_mean = mean(abs), abs_std = sd(abs)) %>%
    unite(condition, c('substrate', 'temp')) %>%
    ggplot(aes(x=time, y=abs_mean, color=condition, fill=condition)) + 
    geom_line(size=0.8) +
    geom_ribbon(aes(ymin=abs_mean - abs_std, ymax=abs_mean+abs_std), alpha=0.3, linetype=0) +
    facet_wrap(~sample) +
    labs(title='gapN NMN Traces', x='Time (s)', y=expression('OD'[340]), color='Condition') +
    theme_half_open(12) +
    background_grid(color.major='gray94') +
    panel_border() +
    theme(strip.background = element_rect(fill='gray91')) +
    theme(legend.position = c(0.55, 0.18)) +
    theme(panel.spacing = unit(0.7, 'lines')) + 
    scale_color_manual(values = c('#FF2C00', '#845B97', '#00B945')) +
    scale_fill_manual(values = c('#FF2C00', '#845B97', '#00B945')) +
    guides(fill=FALSE) + 
    guides(color=guide_legend(override.aes=list(fill=NA)))
         
save_plot('gapn_nmn_traces.pdf', nmn, base_width=9, base_height=7)

