library(tidyverse)
library(cowplot)
library(XML)
library(scales)

parse_xml <- function(xml_path) {
  
  xml <- xmlParse(xml_path)
  abs <- xpathApply(xml, '//Well/RawData', xmlValue)
  wave <- xpathApply(xml, '//Well/WaveData', xmlValue)
  names <- xpathApply(xml, '//Well')
  names <- unlist(map(names, xmlGetAttr, 'Name'))
  
  # split single strings on space delim
  abs <- str_split(abs, ' ')
  wave <- str_split(wave, ' ')
  
  # combine waves and abs
  dfs = c()
  for (i in 1:length(abs)) {
    temp_df <- tibble(wave = wave[[i]],
                      abs = abs[[i]]) 
    temp_df <- rename(temp_df, !!names[[i]] := abs)
    dfs <- append(dfs, temp_df)
  }
  
  df <- bind_cols(dfs) %>%
    select(c('wave...1', 'H9', 'H11', 'H12')) %>%
    rename(wave = 'wave...1') %>%
    mutate_if(is.character, as.numeric) %>%
    rename('Tris-Cl blank' = 'H9') %>%
    rename('1mM NMNH' = 'H11') %>%
    rename('1mM NMN' = 'H12')
  
  return (df)
}

# fix exitation 382
fix_ex <- parse_xml('1-18-21_nmn_scan/T4.xml')
fix_ex$run = 'Fixed Excitation (382nm)'

# fix emission 445
fix_em <- parse_xml('1-18-21_nmn_scan/T5.xml')
fix_em$run = 'Fixed Emission (445nm)'

# concat and to tidy
df <- bind_rows(fix_ex, fix_em) %>%
  pivot_longer(c(-run, -wave))

# plot
# does not allow individual scaling of ylimits
g <- df %>%
  ggplot(aes(x = wave, y=value, color=name)) +
  geom_line(size=0.9) +
  facet_wrap(~run, scales='free') +
  labs(title='Wavelength Scan', color='Sample', x='Wavelength (nm)', y='Fluorescence (RFU)') +
  theme_half_open(14) +
  background_grid() +
  panel_border() +
  theme(legend.position='bottom') +
  theme(strip.background =element_rect(fill="grey92")) +
  scale_y_continuous(labels = scales::comma)

#ggsave('wavescan.png', width=10, height=5)

# grid separate plots
ex <- df %>%
    filter(run == 'Fixed Excitation (382nm)') %>%
    ggplot(aes(x=wave, y=value, color=name)) + 
    geom_line(size=0.9) +
    scale_color_manual(values = c('#0C5DA5', '#00B945', '#FF9500')) +
    theme_half_open(13) +
    labs(title='Fixed Excitation (382nm)', color='Sample', x='Wavelength (nm)', y='Fluorescence (RFU)') +
    background_grid() +
    panel_border() +
    coord_cartesian(ylim=c(0, 2000)) +
    scale_y_continuous(labels = scales::comma) +
    theme(plot.title = element_text(hjust=0.5))

em <- df %>%
    filter(run == 'Fixed Emission (445nm)') %>%
    ggplot(aes(x=wave, y=value, color=name)) +
    geom_line(size=0.9) +
    scale_color_manual(values = c('#0C5DA5', '#00B945', '#FF9500')) +
    theme_half_open(13) +
    labs(title='Fixed Emission (445nm)', x='Wavelength (nm)', y='') +
    background_grid() +
    panel_border() +
    theme(legend.position='none') +
    scale_y_continuous(labels = scales::comma) +
    theme(plot.title = element_text(hjust=0.5))

legend <- get_legend(ex + 
                         guides(color=guide_legend(nrow=1)) +
                         theme(legend.position ='bottom'))
prow <- plot_grid(ex + theme(legend.position = 'none'), em, nrow=1, labels=c('A', 'B'))
prow <- plot_grid(prow, legend, ncol=1, rel_heights = c(1, 0.1))

save_plot('wavescan_plots.pdf', prow, base_height=4, base_width=8.5)
