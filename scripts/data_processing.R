# =============================================================================
# Data Processing and Cox Analysis
# Neutrophil-Related Gene Signature for ccRCC Prognosis
# =============================================================================


setwd('C:\\SYH\\Tool\\R studio\\ccRCC\\KIRC\\1.TCGA')


install.packages('tidyverse')
library(tidyverse) ##f0f
.e$g

exp <- read.csv('TCGA-KIRC.htseq_counts.tsv', sep = '\t', header = T) %>% 
  rename(id = Ensembl_ID)

id <- read.table('gencode.v22.annotation.gene.probeMap', sep = '\t', header = T) %>% 
  select(id, gene)  

exp.1 <- id %>% 
  inner_join(exp) %>% 
  select(-id) %>% 
  filter(!duplicated(gene)) %>% 
  column_to_rownames(var = 'gene') 

exp.1 <- 2 ^ exp.1 - 1 
colnames(exp.1) <- gsub('\\.', '-', colnames(exp.1))

cln <- read.csv('TCGA-KIRC.GDC_phenotype.tsv', sep = '\t', header = T) %>% 
  select(submitter_id.samples,
         disease_type)

cln %>%
  group_by(disease_type) %>% 
  count()

TM <- filter(cln, submitter_id.samples == str_extract(cln$submitter_id.samples, '^.*(01A)'))

a <- exp.1[ ,colnames(exp.1) %in% str_extract(colnames(exp.1), '^.*(11A)')]
b <- exp.1[ ,colnames(exp.1) %in% TM$submitter_id.samples]
final.exp <- cbind(a, b)

group <- data.frame('sample' = colnames(final.exp),
                    'group' = c(rep('CT', 72), rep('tumor', 526)))
write.csv(group, 'group.csv')
exprm <- final.exp
write.csv(exprm, 'expmatrix.csv')
write.csv(b, 'expmatrix.TM.csv')

##e7.e<ef
setwd('C:\\SYH\\Tool\\R studio\\ccRCC\\KIRC\\1.TCGAA')
library(tidyverse)
exprm <- read.csv('expmatrix.csv', row.names = 1)
group <- read.csv('group.csv', row.names = 1)
Group <- factor(group$group, levels = c("CT","tumor"))
design <- model.matrix(~Group)
colnames(design) <- c("CT","tumorvsCT")

##edgeR
install.packages('BiocManager')
BiocManager::install('edgeR')

library(edgeR)
y <- DGEList(counts = exprm, group = Group)
keep <- filterByExpr(y)
y <- y[keep,,keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
y <- estimateDisp(y, design)
fit <- glmFit(y, design)
lrt <- glmLRT(fit, coef=2)
DEG <- topTags(lrt,
               n = Inf,
               adjust.method = 'BH')

DEG_list <- DEG$table %>%
  mutate(type = ifelse(abs(logFC) > 1.5 & FDR < 0.05, ifelse(logFC > 1.5, 'up', 'down'), 'normal'))

DEG_list %>% 
  group_by(type) %>% 
  count()

DEG2 <- DEG_list %>%
  filter(abs(logFC) > 1.5 & FDR < 0.05)

##g+e11e>
library(ggplot2)
ggplot(DEG_list)+
  geom_point(aes(logFC, -log10(FDR), 
                 colour = type), 
             show.legend = T)+
  theme_bw()+
  scale_colour_manual(values = c('down' = 'darkgreen', 'up' = 'darkred', 'normal' = 'black'))

write.csv(DEG2, 'DEG.1.5.csv')

##VENNe>
gene <- read.csv('GeneCards-SearchResults.csv') %>% 
  filter(Category == 'Protein Coding') %>% 
  filter(Relevance.score > 1)
DEG2 <- read.csv('DEG.1.5.csv', row.names = 1)

lista <- rownames(DEG2)
listb <- gene$Gene.Symbol
vennlist <- list(A = lista,
                 B = listb)
##venne>
library(VennDiagram)
P1 = venn.diagram(vennlist,
                  filename = 'overlap.tiff',
                  imagetype = 'tiff',
                  fill = c("yellow", "pink"),
                  category.names = c('DEG', 'Neutrophil'), 
                  cat.pos = c(-15, 15),
                  scale = F)

over <- VennDiagram::get.venn.partitions(vennlist)
overdf <- data.frame('Gene' = over$..values..$`1`)
write.csv(overdf, 'overlap_gene.csv')

###f0f
.e$g
mat.final <- b[rownames(b) %in% overdf$Gene,] %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column(var = 'sample')
mat.final$sample <- gsub('\\.', '-', mat.final$sample)

##ege-f0f
.ee96
surv <- read.table('TCGA-KIRC.survival.tsv', header = T, sep = '\t') %>% 
  rename(status = OS,
         time = OS.time)

merge <- surv %>% 
  inner_join(mat.final) %>% 
  select(-X_PATIENT) %>% 
  column_to_rownames(var = 'sample')

write.csv(merge, 'coxinput.csv') 

