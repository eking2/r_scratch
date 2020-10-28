library(tidyverse)
library(XML)
library(bio3d)
library(cowplot)

#####################################################################################

# get pdb ids for target protein

find_pdbs <- function(uniprot) {
  
  url <- paste0('https://www.uniprot.org/uniprot/', uniprot, '.xml')
  file_path <- paste0(uniprot, '.xml')
  
  # do not redownload
  if (!file.exists(file_path)) {
    download.file(url, file_path)
  }
  
  # filter from default uniprot namespace
  xml <- xmlParse(file_path)
  nd <- getNodeSet(xml, "//ns:dbReference[@type='PDBsum']/@id", namespaces=c(ns=getDefaultNamespace(xml)[[1]]$uri))
  pdbs <- unname(unlist(nd))
  
  return (pdbs)
}

#####################################################################################

# get distance between residue alpha carbons

calc_dist <- function(v1, v2) {
  return (sqrt( (v1$x - v2$x)^2 + (v1$y - v2$y)^2 + (v1$z - v2$z)^2))
}

get_dist <- function(pdb_id, resi_1, resi_2) {
  
  pdb <- read.pdb(pdb_id)
  
  # get atom indices and filter from all atoms
  sel_idx <- atom.select(pdb, elety='CA', resno=c(resi_1, resi_2)) 
  sel <- pdb$atom[sel_idx$atom, ]
  
  # subset xyz coords and calc norm
  v1 = select(sel[1, ], c(x, y, z))
  v2 = select(sel[2, ], c(x, y, z))
  dist <- calc_dist(v1, v2)
  
  return (dist)
}

#####################################################################################

# p450cam
# channel opening, N59 to S190
pdbs <- find_pdbs('P00183')
dists <- map(pdbs, get_dist, 59, 190)

df <- tibble(pdb = pdbs,
             dist = unlist(dists))

p <- df %>%
  ggplot(aes(x = dist)) + 
 # geom_histogram(color='blue', fill='lightblue') +
  geom_density(size=0.7, fill='lightblue', alpha=0.5) +
  theme_half_open(12) +
  background_grid() +
  labs(title='P450cam Channel Opening', x=expression(paste("Distance ", (ring(A)))), y='Density',
       subtitle='ASN-59 to SER-190') +
  geom_vline(xintercept=21, color='red', size=1.2) +
  annotate('text', x=26, y=0.21, label='Open', size=5) +
  annotate('text', x=18.2, y=0.89, label='Closed', size=5) +
  coord_cartesian(xlim=c(15, 29))

ggsave('p450cam_channel.pdf', p, width=6, height=4)

#####################################################################################

# plot rama

pdb <- read.pdb(pdbs[1])
tor <- torsion.pdb(pdb)
tor_res <- select(as.tibble(tor$tbl), c(phi, psi))

p <- tor_res %>%
  drop_na() %>%
  ggplot(aes(x = phi, y = psi)) + 
  stat_density_2d(aes(fill=..level..), geom='polygon', bins=100, alpha=0.1) + 
  geom_point(alpha=0.4) +  
  labs(title='Ramachandran Plot', subtitle='1AKD', x=expression(paste('Phi (', phi, ')')), y=expression(paste('Psi (', psi, ')'))) +
  theme_half_open(12) +
  background_grid() +
  theme(legend.position = 'none') +
  scale_fill_viridis_c()

ggsave('rama.pdf', p, width=6, height=5)

#####################################################################################

# plot bfac

ca_inds <- atom.select(pdb, "calpha")
bf <- as.tibble(pdb$atom[ca_inds$atom,])
  
p <- bf %>%
  ggplot(aes(x=resno, y=b)) +
  annotate('rect', xmin=165, xmax=219, ymin=12, ymax=50, alpha=0.4, fill='lightgreen') +
  geom_line() +
  theme_half_open(12) + 
  background_grid() +
  labs(title='B-factors', x='Residue', subtitle='1AKD', y=expression(paste("B ", (ring(A)^2)))) +
  annotate('text', x=260, y=45, label='FG Region')
  
ggsave('bfa.pdf', p, width=6.5, height=4)
  
  