library(tidyverse)
library(msa)
library(seqinr)
library(ape)
library(ggtree)
library(cowplot)

#####################################################################

# plot clustal aligned msa 

seqs <- readAAStringSet('bdh_sequences.fasta')
aln <- msa(seqs, 'ClustalOmega', order='input') 

msaPrettyPrint(aln, output='pdf', verbose=TRUE, askForOverwrite=FALSE, shadingMode = 'similar')

# save output
aln <- msaConvert(aln, type='seqinr::alignment')
write.fasta(as.list(aln$seq), aln$nam, 'bdh_aln.fasta', nbchar=60, as.string=TRUE)

#####################################################################

# sequence distance matrix

# dissimilarity
d <- dist.alignment(aln, 'identity')
mat <- 1 - as.matrix(d)
write.csv(mat, 'bdh_dist_mat.csv')

# cluster 
order <- hclust(d)$order
mat_ordered <- mat[order, order]
samp_names <- attr(mat_ordered, 'dimnames')[[1]]
samp_names = str_replace_all(samp_names, '_', ' ')

tib <- mat_ordered %>%
  as_tibble() %>%
  mutate(sample_y = samp_names) %>% 
  gather('sample', 'id', -sample_y) 

tib$sample <- factor(x=str_replace_all(tib$sample, '_', ' '), levels = samp_names, ordered=TRUE)
tib$sample_y <- factor(x=str_replace_all(tib$sample_y, '_', ' '), levels = rev(samp_names), ordered=TRUE)

p <- tib %>%
  ggplot(aes(x=sample, y=sample_y)) +
  geom_tile(aes(fill=id)) +
  labs(x='', y='', title='Sequence Similarity Matrix', fill='Identity') +
  theme_half_open(12) +
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  geom_text(aes(label = round(id, 2), color = ifelse(id < 0.5, 'low', 'high'))) + 
  scale_fill_viridis_c() +
  scale_color_manual(values=c(low = 'white', high='black'), guide='none')

ggsave('id_matrix.pdf', p, width=8)

#####################################################################

# phylogenetic tree

aln_tree <- nj(d)
aln_tree$tip.label <- str_replace_all(aln_tree$tip.label, '_', ' ')
p <- ggtree(aln_tree) +
  geom_tiplab(align=TRUE) +
  labs(title='BDO Phylogenetic Tree') +
  ggplot2::xlim(0, 0.7) +
  geom_text2(aes(subset=!isTip, label=node), hjust=-0.3) +
  geom_hilight(14, 'lightblue', alpha=0.2) +
  geom_hilight(18, 'lightgreen', alpha=0.2) +
  theme_tree2()

ggsave('example_tree.pdf', p)
