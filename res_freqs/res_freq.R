library(tidyverse)
library(cowplot)

df <- read_csv('1qi1_msa_freq.csv')
entro <- read_csv('1qi1_entropy.csv')

# drop position numbering (use pdb)
# ignore gaps
# normalize freq each position

df_clean <- df %>%
    select(-c(position, '-')) %>%
    pivot_longer(-pdb_num, names_to='amino_acid',, values_to='count') %>%
    group_by(pdb_num) %>%
    mutate(freq = count / sum(count))

# amino acid distribution bars at selected residue positions

labels <- c('179' = 'P179K', '153' = 'F153S', '210' = 'G210Q', '214' = 'G214E', 
            '234' = 'I234E', '326' = 'I326KR', '330' = 'S330R')

p <- df_clean %>%
    filter(pdb_num %in% c(153, 179, 210, 214, 234, 326, 330)) %>%
    ggplot(aes(x = amino_acid, y=freq, fill=amino_acid)) +
    geom_bar(stat='identity', color='black') +
    facet_wrap(~pdb_num, nrow=2, labeller=labeller(pdb_num = labels), scales='free_x') +
    labs(x='', y='Frequency', title='Streptococcus mutans gapN', subtitle='MSA amino acid distributions') +
    theme_half_open(14) +
    background_grid() +
    panel_border() +
    scale_fill_viridis_d() +
    theme(legend.position = 'none') +
    theme(strip.background = element_rect(fill='grey92')) +
    theme(panel.spacing = unit(0.7, 'lines'))
    
ggsave('gapn_freqs.pdf', p, width=16, height=10)    

# entropy along length of protein

#close <- c(150,151,152,153,154,155,159,177,178,179,180,181,208,209,210,214,
#           215,217,218,228,229,230,231,232,234,237,238,241,250,251,252,284,377,379)

p2 <- entro %>%
    ggplot(aes(x = pdb_num, y = shannon)) +
    #geom_point(color='black') +
    #geom_segment(aes(x = pdb_num, xend = pdb_num, y=0, yend=shannon), alpha=0.2) +
    geom_line() + 
    theme_half_open(14) +
    background_grid() +
    labs(title='Sm gapN sequence conservation', y='Shannon Entropy (bits)' ,x='Residue') +
    annotate('rect', xmin=150, xmax=159, ymin=0, ymax=3.5, fill='green', alpha=0.3) +
    annotate('rect', xmin=177, xmax=181, ymin=0, ymax=3.5, fill='green', alpha=0.3) +
    annotate('rect', xmin=208, xmax=210, ymin=0, ymax=3.5, fill='green', alpha=0.3) +
    annotate('rect', xmin=214, xmax=218, ymin=0, ymax=3.5, fill='green', alpha=0.3) +
    annotate('rect', xmin=229, xmax=234, ymin=0, ymax=3.5, fill='green', alpha=0.3) +
    annotate('rect', xmin=237, xmax=241, ymin=0, ymax=3.5, fill='green', alpha=0.3) +
    annotate('rect', xmin=250, xmax=252, ymin=0, ymax=3.5, fill='green', alpha=0.3) +
    annotate('rect', xmin=283, xmax=284, ymin=0, ymax=3.5, fill='green', alpha=0.3) +
    annotate('rect', xmin=377, xmax=379, ymin=0, ymax=3.5, fill='green', alpha=0.3) 
    
    
ggsave('gapn_entropy.pdf', p2, width=10, height=6)

# sequence heatmap

df_clean$amino_acid <- factor(df_clean$amino_acid, levels = c('A', 'V', 'L', 'I', 'M', 'C',
                                                              'F', 'Y', 'W', 'H', 'S', 'T',
                                                              'N', 'Q', 'D', 'E', 'K', 'R',
                                                              'G', 'P'))

p3 <- df_clean %>%
    ggplot(aes(x = pdb_num, y=amino_acid, fill=freq)) +
    geom_tile() +
    labs(x='Residue', y='Amino Acid', fill='Frequency', title='Streptococcus mutans gapN') +
    scale_fill_viridis_c() +
    scale_y_discrete(limits = rev(levels(df_clean$amino_acid)), expand=c(0, 0)) +
    theme_half_open(14) +
    scale_x_continuous(expand=c(0, 0))

ggsave('gapn_seq_heat.pdf', p3, width=12, height=6)

