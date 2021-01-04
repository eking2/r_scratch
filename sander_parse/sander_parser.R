library(tidyverse)
library(stringr)
library(glue)
library(cowplot)

parseSander <- function(filepath, term) {

    connection <- file(filepath, 'r')
    line <- readLines(connection, n=1)
    res <- vector()
    
    while(length(line) > 0) {
        
        # find term = value, no match will return NA
        regex <- glue('\\s{term}\\s+=\\s+\\-?\\d+\\.?\\d+?')
        save_term <- str_extract(line, regex)
        
        # get energy value, do not need the label
        if (!is.na(save_term))  {
            split <- str_split(save_term, ' ')
            value <- split[[1]][length(split[[1]])]
            res <- append(res, value)
            
        # ignore averages
        } else if (str_detect(line, 'A V E R A G E')) {
            break
        }
        
        # step new line
        line <- readLines(connection, n=1)
    }
    
    close(connection)
    return(res)
}
    
fn = '05_equil_1C5Y.out'

# not very efficient, reads the file once per term
df <- tibble(
    NSTEP = parseSander(fn, 'NSTEP'),
    TEMP = parseSander(fn, 'TEMP\\(K\\)'),
    PRESS = parseSander(fn, 'PRESS'),
    Etot = parseSander(fn, 'Etot'),
    EKtot = parseSander(fn, 'EKtot'),
    EPtot = parseSander(fn, 'EPtot'),
    BOND = parseSander(fn, 'BOND'),
    ANGLE = parseSander(fn, 'ANGLE'),
    DIHED = parseSander(fn, 'DIHED'),
    NB14 = parseSander(fn, '1-4 NB'),
    EEL14 = parseSander(fn, '1-4 EEL'),
    VDWAALS = parseSander(fn, 'VDWAALS'),
    EELEC = parseSander(fn, 'EELEC'),
    EKCMT = parseSander(fn, 'EKCMT'),
    VIRIAL = parseSander(fn, 'VIRIAL'),
    VOLUME = parseSander(fn, 'VOLUME'),
    Density = parseSander(fn, 'Density')
) %>%
    mutate_if(is.character, as.numeric)

write.csv(df, 'sander_out.csv')

summary <- df %>%
    summarize_all(list(mean=mean, sd=sd))

p <- df %>%
    mutate(time = NSTEP / 500000)  %>% # to ns
    select(-NSTEP) %>%
    pivot_longer(!time, names_to='type', values_to='energy') %>% # to tidy
    ggplot(aes(x = time, y=energy)) +
    geom_line(alpha=0.2) +
    geom_smooth(aes(color=type), se=FALSE, size=1.5) + 
    facet_wrap(~type, scales='free_y') + 
    theme_half_open(13) +
    labs(x='Time (ns)', y='Energy (kcal/mol)', title='SANDER Energies') +
    theme(legend.position = 'none') + 
    background_grid() +
    theme(strip.background = element_rect(fill='grey95')) +
    panel_border()

ggsave('sander.pdf', p, width=14, height=9)
