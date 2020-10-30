library(tidyverse)
library(cowplot)
library(ggrepel)

# combine all data
hids <- read_csv('all_hid.csv', col_names = c('Sample', 'All HID'), skip=1)
hips <- read_csv('all_hip.csv', col_names = c('Sample', 'All HIP'), skip=1)
small_hid <- read_csv('small_hid.csv', col_names = c('Sample', 'Small HID'), skip=1)
small_hip <- read_csv('small_hip.csv', col_names = c('Sample', 'Small HIP'), skip=1)
exp <- read_csv('experiment.csv', col_names = c('Sample', 'Experiment'), skip=1)

# drop sampling outlier
df <- exp %>%
    inner_join(hids) %>%
    inner_join(hips) %>%
    inner_join(small_hid) %>%
    inner_join(small_hip) %>%
    filter(Sample != '1GJD')

calc_r_rmse <- function(y_pred, y_true) {
    
    diff_sq <- (y_true - y_pred)^2
    rmse <- sqrt(mean(diff_sq))
    
    pearson <- cor(y_pred, y_true, method='pearson')
    
    return (tibble(rmse = rmse,
                   cor = pearson))
}

# metrics for each condition
stats <- df %>%
    select(-Experiment, -Sample) %>%
    map(~calc_r_rmse(., df$Experiment)) %>%
    bind_rows(.id = 'Condition')

df_tidy <- df %>%
    pivot_longer(!Sample, names_to='Condition', values_to='Simulation') %>%
    inner_join(exp) %>%
    filter(Condition != 'Experiment') %>%
    mutate(Size = ifelse(Sample %in% c('1C5X', '1C5Y', '1C5Z', '1GI7'), 'Small', 'Big'))

# set facet order
df_tidy$Condition <- factor(df_tidy$Condition, levels=c('All HIP', 'All HID', 'Small HIP', 'Small HID'))
df_tidy$Size <- factor(df_tidy$Size, levels=c('Small', 'Big'))
stats$Condition <- factor(stats$Condition, levels=c('All HIP', 'All HID', 'Small HIP', 'Small HID'))

p <- df_tidy %>%
    ggplot(aes(x = Experiment, y=Simulation, label=Sample)) + 
    geom_abline(intercept=0, slope=1, color='orange', size=1.5) +
    geom_point(aes(fill=Size), shape=21, size=3) + 
    facet_wrap(~Condition, ncol=2) +
    theme_half_open(15) +
    background_grid() +
    panel_border() + 
    theme(strip.background = element_rect(fill='gray92')) +
    labs(title='Binding Free Energy Predictions', subtitle='PBSA/MBAR', x='Experiment (kcal/mol)', y='Simulation (kcal/mol)') +
    scale_fill_manual(values=c('green1', 'deepskyblue1')) + 
    geom_text_repel(point.padding=0.1) +
    coord_cartesian(xlim=c(-12, -3.5), ylim=c(-12, -3.5)) +
    geom_text(data=stats, aes(x=-11.8, y=-4.6, label=sprintf('RMSE = %0.2f', round(rmse, 2))), size=6, hjust=0) +
    geom_text(data=stats, aes(x=-11.8, y=-5.3, label=sprintf('R = %0.2f', round(cor, 2))), size=6, hjust=0)
    
ggsave('binding_preds.pdf', p, width=12, height=10)
