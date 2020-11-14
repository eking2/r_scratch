library(tidyverse)
library(broom)
library(fs)
library(XML)
library(gtools)
library(cowplot)
library(pals)

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

# get slope for each well, 0-300s
# slope from abs/sec to mAbs/min
slopes <- df %>%
    filter(time < 300) %>%
    group_by(trial, well) %>%
    nest() %>%
    mutate(model = map(data, ~lm(abs ~ time, data=.x)),
           tidied = map(model, tidy)) %>%
    unnest(tidied) %>%
    select(well, trial, term, estimate) %>%
    pivot_wider(names_from=term, values_from=estimate) %>%
    rename(intercept = `(Intercept)`,
           slope = time) %>%
    mutate(slope = slope * 60 * 1000)

# average slope over all trials
avg_slopes <- slopes %>%
    group_by(well) %>%
    summarise(mean_slope = mean(slope))  %>%
    mutate(col = str_extract(well, '\\d+')) %>%
    mutate(row_name = str_sub(well, end=1)) %>%
    mutate(row = match(row_name, map_chr(letters, toupper))) %>%
    mutate(col = fct_relevel(col, mixedsort(unique(col)))) %>%
    mutate(well = fct_relevel(well, mixedsort(unique(well)))) %>%
    arrange(well)

# highlight active samples
active <- avg_slopes %>%
    filter(mean_slope > 1.5) %>%
    filter(!well %in% c('A1', 'A12', 'H1', 'H12'))

# look up char to position on plot
#plate_rows <- sort(unique(avg_slopes$row_name))
#plate_cols <- sort(levels(avg_slopes$col))
    
# plot average slopes
p <- avg_slopes %>%
    ggplot(aes(x = col, y = reorder(row_name, desc(row_name)), fill=mean_slope)) +
    geom_tile() +
    geom_tile(data = active, color='red', size=1) +
    geom_text(aes(label=round(mean_slope, 2), color=ifelse(mean_slope > 2.1, 'black', 'white')), size=4) + 
    #scale_x_discrete(position = 'top') +
    labs(title='PTDH Specific Activitiy Slopes', subtitle='GOR Selection Library 3', x='', y='', fill='Slope (mA/min)') +
    theme_half_open(14) +
    scale_fill_gradientn(colors = parula(100)) +
    scale_color_manual(values = c(black = 'black', white = 'white'), guide='none')


# save
write_csv(avg_slopes, 'avg_slopes.csv')
write_csv(slopes, 'slopes.csv')
ggsave('slope_heat.pdf', p, width = 10, height=6)
