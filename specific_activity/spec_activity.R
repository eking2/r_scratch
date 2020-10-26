library(tidyverse)
library(cowplot)
library(ggsci)

slopes <- read.csv('slopes.csv')
enz <- read.csv('enz_conc.csv')

#######################################################################

# get specific activities from abs slopes, blank subtract

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


# separate no sub run
nosub <- slopes %>%
  filter(cofa == 'no') %>%
  select(sample, slope) %>%
  rename(blank_slope = slope)

# specific activity for each trace
df <- slopes %>%
  inner_join(nosub, by='sample') %>%
  inner_join(enz, by=c("sample" = 'sample_num')) %>%
  filter(cofa != 'no') %>%
  mutate(activity = calc_activity(slope, blank_slope, DF, conc_mg_ml))

write.csv(df, 'activities.csv', row.names=FALSE)

# agg data
df_grp <- df %>%
  group_by(sample, cofa, enzyme) %>%
  summarise(mean=round(mean(activity), 3), sd=round(sd(activity), 3)) %>%
  arrange(cofa, sample)

df_grp$cofa <- factor(df_grp$cofa, levels=c('NMN', 'NAD', 'NADP'))
df_grp$enzyme <- str_replace_all(df_grp$enzyme, '_', ' ')

write.csv(df_grp, 'activities_agg.csv', row.names=FALSE)

#######################################################################

# plot

p <- df_grp %>%
  ggplot(aes(x=enzyme, y=mean, fill=enzyme)) + 
  geom_bar(stat='identity', width=0.65) +
  geom_errorbar(aes(ymin=mean - sd, ymax=mean + sd), width=0.15) +
  facet_wrap(~cofa, scales='free') +
  theme_half_open(12) + 
  theme(legend.position = 'none') +
  theme(strip.background = element_rect(fill='grey92')) +
  theme(axis.text.x = element_text(angle=45, hjust=1)) + 
  background_grid() +
  labs(title='Phosphite dehydrogenase', x='', y='Specific activity (nmol/min/mg)') +
  scale_fill_npg()

ggsave('enz_plot.pdf', p, width=11, height=5)
