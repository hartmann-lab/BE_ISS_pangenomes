######################################
### Comparative Pangenome Analysis ###
######################################

# Ryan Blaustein
# ryan.blaustein@northwestern.edu
# last edited: 12/04/2018

######################################

### set working directory to */BE_ISS_pangenomes
### call packages
library(ggplot2)
library(vegan)
library(reshape2)
library(minpack.lm)
library(ape)
library(phangorn)
library(gplots)
library(heatmap.plus)
library(VennDiagram)

### load pangenome matrices from roary output

## B. cereus
pangenome_bac = read.table("pangenome_roary_outputs/roary_out_B_cereus/gene_presence_absence.Rtab",
                           h=T, row.names=1)

## S. aureus
pangenome_staph = read.table("pangenome_roary_outputs/roary_out_S_aureus/gene_presence_absence.Rtab",
                             h=T, row.names=1)

### load metadata; for S. aureus, differentiate MRSA from other human-associated strains

## B. cereus
bac_meta = read.table("pangenomes_ad-hoc/R_files/B_cereus/Bac_metadata.txt",
                      sep='\t',h=T,row.names=1,check=F,comment='')

## S. aureus
staph_meta = read.table("pangenomes_ad-hoc/R_files/S_aureus/Staph_metadata.txt",
                        sep='\t',h=T,row.names=1,check=F,comment='')
# add MRSA column based on "health comments"
staph_meta$MRSA = rep("no", dim(pangenome_staph)[2])
staph_meta$MRSA[grep("MRSA", staph_meta$Health_comments)] = rep("yes", length(grep("MRSA", staph_meta$Health_comments)))
staph_meta$Category = as.character(staph_meta$Category)
staph_meta$Category[grep("MRSA", staph_meta$Health_comments)] = rep("MRSA", length(staph_meta$Category[grep("MRSA", staph_meta$Health_comments)]))
staph_meta$Category = as.factor(staph_meta$Category)

################################
### PANGENOME SUMMARY STATS ####
################################

###
### Total genes per genome ###
###

## Bacillus
totalgenes_bac = read.table("pangenome_roary_outputs/roary_out_B_cereus/number_of_genes_in_pan_genome.Rtab")
colnames(totalgenes_bac) = c(1:dim(totalgenes_bac)[2])
totalgenes_bac_m = melt(totalgenes_bac)
# plot
ggplot(data = totalgenes_bac_m, aes(x = variable, y = value)) +
  geom_boxplot(outlier.shape = NA, col = "blue") +
  xlab("No. of genomes") +
  ylab("Total genes") +
  scale_x_discrete(breaks=c(0, 10, 20, 30, 40, 50)) +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  theme(axis.text = element_text(size=13, color = "black"), 
        axis.title = element_text(size=14, color = "black"))
# quick stats (mean, sd, se)
mean(apply(pangenome_bac, 2, sum))
sd(apply(pangenome_bac, 2, sum))
sd(apply(pangenome_bac, 2, sum))/sqrt(dim(pangenome_bac)[2])

## Staphylococcus
totalgenes_staph = read.table("pangenome_roary_outputs/roary_out_S_aureus/number_of_genes_in_pan_genome.Rtab")
colnames(totalgenes_staph) = c(1:dim(totalgenes_staph)[2])
totalgenes_staph_m = melt(totalgenes_staph)
# plot
ggplot(data = totalgenes_staph_m, aes(x = variable, y = value)) +
  geom_boxplot(outlier.shape = NA, col = "red") +
  xlab("No. of genomes") +
  ylab("Total genes") +
  scale_x_discrete(breaks=c(0, 25, 50, 75, 100)) +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  theme(axis.text = element_text(size=13, color = "black"), 
        axis.title = element_text(size=14, color = "black"))
# quick stats (mean, sd, se)
mean(apply(pangenome_staph, 2, sum))
sd(apply(pangenome_staph, 2, sum))
sd(apply(pangenome_staph, 2, sum))/sqrt(dim(pangenome_staph)[2])

###
### Permute number of new genes per genome ###
###

## gene count function
new_gene_calc <- function (x) {
  gene = c(1:dim(x)[1])
  for (i in 1:dim(x)[1]) {
    gene[i] = which(x[i,] > 0)[1]
  }
  gene_count = c(1:dim(x)[2]) 
  for (z in 1:dim(x)[2]) {
    gene_count[z] = length(which(gene == z))
  }
  print(gene_count)
}

## 100x permutations

# Bacillus
perm_pan_bac = as.data.frame(matrix(nrow = 100, ncol = dim(pangenome_bac)[2]))
for (n in 1:100) {
  perm_pan_bac[n,] = new_gene_calc(pangenome_bac[,sample(dim(pangenome_bac)[2])])
  print(n)
} 
#write.table(perm_pan_bac, 'pangenomes_ad-hoc/R_files/B_cereus/permutation_new_genes_per_genome.txt', sep="\t", col.names=NA, quote = FALSE)

# Staphylococcus
perm_pan_staph = as.data.frame(matrix(nrow = 100, ncol = dim(pangenome_staph)[2]))
for (n in 1:100) {
  perm_pan_staph[n,] = new_gene_calc(pangenome_staph[,sample(dim(pangenome_staph)[2])])
  print(n)
} 
#write.table(perm_pan_staph, 'pangenomes_ad-hoc/R_files/S_aureus/permutation_new_genes_per_genome.txt', sep="\t", col.names=NA, quote = FALSE)

## load result (instead of re-running)
perm_pan_bac = read.table('pangenomes_ad-hoc/R_files/B_cereus/permutation_new_genes_per_genome.txt', h=T, row.names=1)
perm_pan_staph = read.table('pangenomes_ad-hoc/R_files/S_aureus/permutation_new_genes_per_genome.txt', h=T, row.names=1)

###
### Model the number of new genes per genome ###
###

## B. cereus
newgenes_bac = perm_pan_bac
colnames(newgenes_bac) = c(1:dim(newgenes_bac)[2])
newgenes_bac_m = melt(newgenes_bac)
newgenes_bac_stat = data.frame(genomes = c(1:dim(pangenome_bac)[2]),
                               ng_mean = tapply(newgenes_bac_m$value, newgenes_bac_m$variable, mean),
                               ng_median = tapply(newgenes_bac_m$value, newgenes_bac_m$variable, median),
                               ng_sd = tapply(newgenes_bac_m$value, newgenes_bac_m$variable, sd),
                               ng_se = tapply(newgenes_bac_m$value, newgenes_bac_m$variable, sd)/sqrt(100))
# Power-law fit 
ng_bac_power = nlsLM(log(ng_mean)~a*genomes^-b, 
                     data=newgenes_bac_stat, start = list(a=1,b=1))
summary(ng_bac_power)  
ng_bac_power_predict = data.frame(var = 1:dim(pangenome_bac)[2],
                                  val = predict(ng_bac_power))
# N50 value
exp(ng_bac_power_predict$val[50])

## S. aureus
newgenes_staph = perm_pan_staph
pangenome_staph = pangenome_staph
colnames(newgenes_staph) = c(1:dim(newgenes_staph)[2])
newgenes_staph_m = melt(newgenes_staph)
newgenes_staph_stat = data.frame(genomes = c(1:dim(pangenome_staph)[2]),
                                 ng_mean = tapply(newgenes_staph_m$value, newgenes_staph_m$variable, mean),
                                 ng_median = tapply(newgenes_staph_m$value, newgenes_staph_m$variable, median),
                                 ng_sd = tapply(newgenes_staph_m$value, newgenes_staph_m$variable, sd),
                                 ng_se = tapply(newgenes_staph_m$value, newgenes_staph_m$variable, sd)/sqrt(100))
# Power-law fit 
ng_staph_power = nlsLM(log(ng_mean)~a*genomes^-b,
                       data=newgenes_staph_stat, start = list(a=1,b=1)) 
summary(ng_staph_power)
ng_staph_power_predict = data.frame(var = 1:dim(pangenome_staph)[2],
                                    val = predict(ng_staph_power))
exp(ng_staph_power_predict$val[50])

## Plot both taxa
ggplot(data = newgenes_staph_m, aes(x = variable, y = value)) +
  geom_point(alpha = 0.1, shape = 1, col = "red", size = 0.9) +
  geom_point(data = newgenes_staph_stat, aes(x = genomes, y = ng_mean), 
             shape = 1, stroke = 1.0,
             size = 1.8, col = "red") +
  geom_line(data = ng_staph_power_predict, aes(x = var, y = exp(val)), 
            lty = 1, col = "red", size = 0.8) +
  geom_point(data = newgenes_bac_m, aes(x = variable, y = value),
             alpha = 0.1, shape = 1, col = "blue", size = 0.9) +
  geom_point(data = newgenes_bac_stat, aes(x = genomes, y = ng_mean), 
             shape = 1, stroke = 1.0,
             size = 1.8, col = "blue") +
  geom_line(data = ng_bac_power_predict, aes(x = var, y = exp(val)), 
            lty = 1, col = "blue", size = 0.8) +
  theme_classic() +
  xlab("No. of genomes") +
  ylab("New genes") +
  scale_x_discrete(breaks = c(0, 25, 50, 75, 100)) +
  scale_y_log10(breaks=c(1,5,10,50,100,500,1000,5000,10000), limits = c(1,10000)) +
  theme(axis.text = element_text(size=16, color = "black"), 
        axis.title = element_text(size=17, color = "black"))

###
### Distribution of pangenome components ###
###

### B. cereus ###

## Histogram 
genes_by_genome_bac <- apply(pangenome_bac, 1, sum)
gene_freq_bac = data.frame(genome = c(1:dim(pangenome_bac)[2]),
                           count = c(1:dim(pangenome_bac)[2]),
                           group = c(rep("Cloud", floor(dim(pangenome_bac)[2]*0.1)),
                                     rep("Shell", ceiling(dim(pangenome_bac)[2]*0.85)),
                                     rep("Core", round(dim(pangenome_bac)[2]*0.05))))
gene_freq_bac$group = factor(gene_freq_bac$group, levels = c("Cloud", "Shell", "Core"))
for (i in 1:dim(pangenome_bac)[2]) {
  gene_freq_bac[i,2] = length(which(genes_by_genome_bac == i))
}
# plot
ggplot(gene_freq_bac, aes(x = genome, y = count, fill = group)) +
  geom_bar(stat = "identity", width = 0.95) +
  xlab("Genome Number") +
  ylab("Total Genes") +
  xlim(0,60) +
  theme_classic() +
  scale_fill_manual(values=c("skyblue", "gold","magenta")) +
  theme(axis.text = element_text(size=15, color = "black"), 
        axis.title = element_text(size=16, color = "black")) +
  theme(legend.position="none")

## Pie chart
pie.df_bac <- data.frame(group = c("Cloud", "Shell", "Core"),
                         count = tapply(gene_freq_bac$count, gene_freq_bac$group, sum),
                         percent = 100*tapply(gene_freq_bac$count, gene_freq_bac$group, sum)/sum(gene_freq_bac$count))
pie.df_bac$percent = round(pie.df_bac$percent, 1)
pie.df_bac$group = factor(pie.df_bac$group, levels = c("Cloud", "Shell", "Core"))
# plot
ggplot(pie.df_bac, aes(x="", y=percent, fill=group)) +
  geom_bar(width = 1, stat = "identity") + 
  coord_polar("y", start=0) +
  scale_fill_manual(values=c("skyblue", "gold","magenta")) +
  #geom_text(aes(y = c(68, 21, 6), 
  #              label = count), size=9) +
  theme_minimal() +
  theme(axis.text.x=element_blank(), axis.title.x = element_blank(), 
        axis.title.y = element_blank(), panel.border = element_blank(), 
        panel.grid=element_blank(), axis.ticks = element_blank()) +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=22, color = "black"))
# total count/percent per fraction
pie.df_bac

### S. aureus ###

## Histogram
genes_by_genome_staph <- apply(pangenome_staph, 1, sum)
gene_freq_staph = data.frame(genome = c(1:dim(pangenome_staph)[2]),
                             count = c(1:dim(pangenome_staph)[2]),
                             group = c(rep("Cloud", floor(dim(pangenome_staph)[2]*0.1)),
                                       rep("Shell", ceiling(dim(pangenome_staph)[2]*0.85)),
                                       rep("Core", floor(dim(pangenome_staph)[2]*0.05))))
gene_freq_staph$group = factor(gene_freq_staph$group, levels = c("Cloud", "Shell", "Core"))
for (i in 1:dim(pangenome_staph)[2]) {
  gene_freq_staph[i,2] = length(which(genes_by_genome_staph == i))
}
# plot
ggplot(gene_freq_staph, aes(x = genome, y = count, fill = group)) +
  geom_bar(stat = "identity", width = 0.95) +
  xlab("Genome Number") +
  ylab("Total Genes") +
  theme_classic() +
  scale_fill_manual(values=c("skyblue", "gold","magenta")) +
  theme(axis.text = element_text(size=15, color = "black"), 
        axis.title.x = element_text(size=16, color = "black"),
        axis.title.y = element_text(size=16, color = "black")) +
  theme(legend.position="none")

## Pie chart
pie.df_staph <- data.frame(group = c("Cloud", "Shell", "Core"),
                           count = tapply(gene_freq_staph$count, gene_freq_staph$group, sum),
                           percent = 100*tapply(gene_freq_staph$count, gene_freq_staph$group, sum)/sum(gene_freq_staph$count))
pie.df_staph$percent = round(pie.df_staph$percent, 1)
pie.df_staph$group = factor(pie.df_staph$group, levels = c("Cloud", "Shell", "Core"))
# plot
ggplot(pie.df_staph, aes(x="", y=percent, fill=group)) +
  geom_bar(width = 1, stat = "identity") + 
  coord_polar("y", start=0) +
  scale_fill_manual(values=c("skyblue", "gold","magenta")) +
  #geom_text(aes(y = c(76, 40, 14), 
  #              label = count), size=9) +
  theme_minimal() +
  theme(axis.text.x=element_blank(), axis.title.x = element_blank(), 
        axis.title.y = element_blank(), panel.border = element_blank(), 
        panel.grid=element_blank(), axis.ticks = element_blank()) +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=22, color = "black"))
# total count/percent per fraction
pie.df_staph

###
### Pangenome fractions by strain origin
###

### B. cereus ###

# Assign genes to core/cloud/shell
gene_type = data.frame(genes_by_genome_bac,
                       category = rep("Core", length(genes_by_genome_bac)))
gene_type$category = as.character(gene_type$category)
gene_type$category[which(gene_type$genes_by_genome_bac < 54)] = 
  rep("Shell", length(gene_type$category[which(gene_type$genes_by_genome_bac < 54)]))
gene_type$category[which(gene_type$genes_by_genome_bac < 6)] = 
  rep("Cloud", length(gene_type$category[which(gene_type$genes_by_genome_bac < 6)]))

# Determine fraction count per genome                              
genes_per_genome_bac = matrix(nrow = 56, ncol = 3)
for (i in 1:dim(pangenome_bac)[2]) {
  genes_per_genome_bac[i,] = rbind(tapply(pangenome_bac[,i], gene_type$category, sum))
}
rownames(genes_per_genome_bac) = rownames(bac_meta)
colnames(genes_per_genome_bac) = names(tapply(pangenome_bac[,1], gene_type$category, sum))
genes_per_genome_bac = as.data.frame(genes_per_genome_bac)
# Compute average
GPG_bac_avg = data.frame(Cloud = tapply(genes_per_genome_bac$Cloud, bac_meta$Category, mean),
                         Shell = tapply(genes_per_genome_bac$Shell, bac_meta$Category, mean),
                         Core = tapply(genes_per_genome_bac$Core, bac_meta$Category, mean))
# Barplot
GPG_bac_gg = melt(t(GPG_bac_avg))
colnames(GPG_bac_gg) = c("gene_type", "sample_type", "gene_count")
ggplot(GPG_bac_gg, aes(x = sample_type, y = gene_count, fill = gene_type)) + 
  geom_bar(stat = "identity", alpha = 0.8) +
  theme_classic() + 
  scale_fill_manual(values = c("skyblue", "gold", "magenta")) +
  ylab("Gene Count") +
  scale_x_discrete(limits = rev(levels(GPG_bac_gg$sample_type))) +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size=24, color = "black", angle = 45, hjust = 1), 
        axis.title.x = element_text(size=24, color = "black"),
        legend.text = element_text(size=24, color = "black"), 
        legend.title = element_blank()) +
  coord_flip()

### S. aureus ###

# Assign genes to core/cloud/shell
gene_type = data.frame(genes_by_genome_staph,
                       category = rep("Core", length(genes_by_genome_staph)))
gene_type$category = as.character(gene_type$category)
gene_type$category[which(gene_type$genes_by_genome_staph < 101)] = 
  rep("Shell", length(gene_type$category[which(gene_type$genes_by_genome_staph < 101)]))
gene_type$category[which(gene_type$genes_by_genome_staph < 11)] = 
  rep("Cloud", length(gene_type$category[which(gene_type$genes_by_genome_staph <11)]))

# Determine fraction count per genome                              
genes_per_genome_staph = matrix(nrow = 105, ncol = 3)
for (i in 1:dim(pangenome_staph)[2]) {
  genes_per_genome_staph[i,] = rbind(tapply(pangenome_staph[,i], gene_type$category, sum))
}
rownames(genes_per_genome_staph) = rownames(staph_meta)
colnames(genes_per_genome_staph) = names(tapply(pangenome_staph[,1], gene_type$category, sum))
genes_per_genome_staph = as.data.frame(genes_per_genome_staph)
# Compute average
GPG_staph_avg = data.frame(Cloud = tapply(genes_per_genome_staph$Cloud, staph_meta$Category, mean),
                           Shell = tapply(genes_per_genome_staph$Shell, staph_meta$Category, mean),
                           Core = tapply(genes_per_genome_staph$Core, staph_meta$Category, mean))

# Barplot
GPG_staph_gg = melt(t(GPG_staph_avg))
colnames(GPG_staph_gg) = c("gene_type", "sample_type", "gene_count")
ggplot(GPG_staph_gg, aes(x = sample_type, y = gene_count, fill = gene_type)) + 
  geom_bar(stat = "identity", alpha = 0.8) +
  theme_classic() + 
  scale_fill_manual(values = c("skyblue", "gold", "magenta")) +
  ylab("Gene Count") +
  scale_x_discrete(limits = rev(levels(GPG_staph_gg$sample_type))) +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size=24, color = "black", angle = 45, hjust = 1), 
        axis.title.x = element_text(size=24, color = "black"),
        legend.text = element_text(size=24, color = "black"), 
        legend.title = element_blank()) +
  coord_flip()

#########################################
### Ordination/Clustering: B. cereus ####
#########################################

###
### GENES ###
###

### PCoA ###

# beta-diversity measure
beta <- vegdist(t(pangenome_bac), 'jaccard', binary = T)
beta.mat <- as.matrix(beta)

# projection
pcoa <- cmdscale(beta, k = 4, eig = TRUE)

# cleanup
ord <- as.data.frame(pcoa$points)
names(ord) <- c('pcoa1', 'pcoa2', 'pco3', 'pco4') ## rename coordinates

# add metadata
ord$Category = bac_meta$Category
ord$Assembler = bac_meta$Assembler2
ord$Seq_Tech = bac_meta$Seq_Tech2
ord$Reference = bac_meta$Reference
ord$Location = bac_meta$Sample_Site_2
ord$Culture_Media = bac_meta$Culture_Media 

# Percent explained variation
eig <- eigenvals(pcoa)
100*head(eig/sum(eig))

# plot PCoA (color by sample category)
ggplot(data = ord, aes(x = pcoa1, y = pcoa2, fill = Category)) +
  geom_point(alpha = 0.65, 
             shape = 21, 
             stroke = 0.5, size = 6) +
  theme_bw() +
  xlab('PCo1 (26.7%)') +
  ylab('PCo2 (10.5%)') +
  theme(axis.text = element_text(size=15, color = "black"),
        axis.title = element_text(size=16, color = "black"),
        legend.text = element_text(size=15, color = "black"),
        legend.title = element_blank()) +
  scale_fill_manual(values=c("gray", "skyblue", "green", "darkgreen", "yellow", "brown"),
                    #guide_legend(title = "Category"),
                    labels=c("BE-E (n=1)", "BE-SC (n=8)", "Cul-E (n=1)", "Cul-SC (n=2)", "Human (n=33)", "Soil (n=11)")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

### PERMANOVA (may need to re-compute distance due to "NA's") ###

# Strain origin
beta <- vegdist(t(pangenome_bac), 'jaccard', binary = T)
adonis(beta ~ Category, data = bac_meta, permutations = 999)

# Assembler
beta <- vegdist(t(pangenome_bac[,-c(which(is.na(bac_meta$Assembler2) == T))]), 'jaccard', binary = T)
adonis(beta ~ Assembler2, data = bac_meta, permutations = 999) 

# Sequencing technology
beta <- vegdist(t(pangenome_bac[,-c(which(is.na(bac_meta$Seq_Tech2) == T))]), 'jaccard', binary = T)
adonis(beta ~ Seq_Tech2, data = bac_meta, permutations = 999) 

# Reference
beta <- vegdist(t(pangenome_bac), 'jaccard', binary = T)
adonis(beta ~ Reference, data = bac_meta, permutations = 999) 

# Location
beta <- vegdist(t(pangenome_bac[,-c(which(is.na(bac_meta$Sample_Site_2) == T))]), 'jaccard', binary = T)
adonis(beta ~ Sample_Site_2, data = bac_meta, permutations = 999) 

# Culture Media
beta <- vegdist(t(pangenome_bac[,-c(which(is.na(bac_meta$Culture_Media) == T))]), 'jaccard', binary = T)
adonis(beta ~ Culture_Media, data = bac_meta, permutations = 999) 

### PCoA: ISS samples only ###

# beta-diversity measure
beta <- vegdist(t(pangenome_bac[,grep("ISS", bac_meta$Category)]), 'jaccard', binary = T)
beta.mat <- as.matrix(beta)

# projection
pcoa <- cmdscale(beta, k = 4, eig = TRUE)

# cleanup
ord <- as.data.frame(pcoa$points)
names(ord) <- c('pcoa1', 'pcoa2', 'pco3', 'pco4') ## rename coordinates

# add metadata
ord$Site = bac_meta[grep("ISS", bac_meta$Category),]$Sample_Site_2

# Percent explained variation
eig <- eigenvals(pcoa)
100*head(eig/sum(eig))

# plot PCoA (color by sample category)
ggplot(data = ord, aes(x = pcoa1, y = pcoa2, shape = Site)) +
  geom_point(alpha = 0.8, 
             stroke = 0.9, 
             size = 7) +
  theme_bw() +
  xlab('PCo1 (45.2%)') +
  ylab('PCo2 (24.5%)') +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_text(size=20, color = "black"),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank()) + 
  #theme(legend.position = "none") +
  theme(legend.text = element_text(size=14, color = "black"),
        legend.title = element_blank(),
        legend.position = "bottom",
        legend.background = element_rect(size = 0.5, colour = 1)) +
  guides(shape=guide_legend(nrow = 4, byrow=F))

### PERMANOVA (may need to re-compute distance due to "NA's") ###

# Culture Media
adonis(beta ~ Culture_Media, data = bac_meta[grep("ISS", bac_meta$Category),], permutations = 999) 

###
### FUNCTIONS ###
###

### Prep data ###

## load "annotated" pangenome matrix
# **prepared from "Gene, Non-unique Gene name, Annotation" columns in gene_presence_absence.csv file
pangenome_anno = read.table("pangenomes_ad-hoc/R_files/B_cereus/pangenome_gene_annotations.txt",
                            sep="\t",header=T,row.names=1, quote="\"")

## bin annotations for overall presence/absence
fxn_df = as.data.frame(matrix(nrow = length(levels(pangenome_anno$Annotation)), ncol = ncol(pangenome_bac)))
for (i in 1:ncol(pangenome_bac)) {
  fxn_df[,i] = tapply(pangenome_bac[,i], pangenome_anno$Annotation,sum)
}
rownames(fxn_df) = names(tapply(pangenome_bac[,1], pangenome_anno$Annotation,sum))
colnames(fxn_df) = colnames(pangenome_bac)

## convert fxn_df to binary presence/absence data
fxn_df_bin = ifelse(fxn_df>0, 1, 0)

## set data frame for B. cereus fxns (used later)
fxn_df_bin_bac = fxn_df_bin

### PCoA ###

# beta-diversity measure
beta <- vegdist(t(fxn_df_bin), 'jaccard', binary = T)
beta.mat <- as.matrix(beta)

# projection
pcoa <- cmdscale(beta, k = 4, eig = TRUE)

# cleanup
ord <- as.data.frame(pcoa$points)
names(ord) <- c('pcoa1', 'pcoa2', 'pco3', 'pco4') ## rename coordinates

# add metadata
ord$Category = bac_meta$Category
ord$Assembler = bac_meta$Assembler2
ord$Seq_Tech = bac_meta$Seq_Tech2
ord$Reference = bac_meta$Reference
ord$Location = bac_meta$Sample_Site_2
ord$Culture_Media = bac_meta$Culture_Media

# Percent explained variation
eig <- eigenvals(pcoa)
100*head(eig/sum(eig))

# plot PCoA (color by sample category)
ggplot(data = ord, aes(x = pcoa1, y = pcoa2, fill = Category)) +
  geom_point(alpha = 0.65, 
             shape = 21, 
             stroke = 0.5, size = 5) +
  theme_bw() +
  xlab('PCo1 (24.0%)') +
  ylab('PCo2 (16.1%)') +
  theme(axis.text = element_text(size=15, color = "black"),
        axis.title = element_text(size=16, color = "black"),
        legend.text = element_text(size=15, color = "black"),
        legend.title = element_blank()) +
  scale_fill_manual(values=c("gray", "skyblue", "green", "darkgreen", "yellow", "brown"),
                    #guide_legend(title = "Category"),
                    labels=c("BE-E (n=1)", "BE-SC (n=8)", "Cul-E (n=1)", "Cul-SC (n=2)", "Human (n=33)", "Soil (n=11)")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

### PERMANOVA (may need to re-compute distance due to "NA's") ###

# Strain origin
beta <- vegdist(t(fxn_df_bin), 'jaccard', binary = T)
adonis(beta ~ Category, data = bac_meta, permutations = 999) 

# Assembler
beta <- vegdist(t(fxn_df_bin[,-c(which(is.na(bac_meta$Assembler2) == T))]), 'jaccard', binary = T)
adonis(beta ~ Assembler2, data = bac_meta, permutations = 999) 

# Sequencing technology
beta <- vegdist(t(fxn_df_bin[,-c(which(is.na(bac_meta$Seq_Tech2) == T))]), 'jaccard', binary = T)
adonis(beta ~ Seq_Tech2, data = bac_meta, permutations = 999) 

# Reference
beta <- vegdist(t(fxn_df_bin), 'jaccard', binary = T)
adonis(beta ~ Reference, data = bac_meta, permutations = 999) 

# Location
beta <- vegdist(t(fxn_df_bin[,-c(which(is.na(bac_meta$Sample_Site_2) == T))]), 'jaccard', binary = T)
adonis(beta ~ Sample_Site_2, data = bac_meta, permutations = 999) 

# Culture Media
beta <- vegdist(t(fxn_df_bin[,-c(which(is.na(bac_meta$Culture_Media) == T))]), 'jaccard', binary = T)
adonis(beta ~ Culture_Media, data = bac_meta, permutations = 999) 

#########################################
### Ordination/Clustering: S. aureus ####
#########################################

###
### GENES ###
###

### PCoA ###

# beta-diversity measure
beta <- vegdist(t(pangenome_staph), 'jaccard', binary = T)
beta.mat <- as.matrix(beta)

# projection
pcoa <- cmdscale(beta, k = 4, eig = TRUE)

# cleanup
ord <- as.data.frame(pcoa$points)
names(ord) <- c('pcoa1', 'pcoa2', 'pco3', 'pco4') ## rename coordinates

# add metadata
ord$Category = staph_meta$Category
ord$MRSA = staph_meta$MRSA
ord$Assembler = staph_meta$Assembler2
ord$Seq_Tech = staph_meta$Seq_Tech2
ord$Reference = staph_meta$Reference
ord$Location = staph_meta$Sample_Site_2
ord$Culture_Media = staph_meta$Culture_Media

# Percent explained variation
eig <- eigenvals(pcoa)
100*head(eig/sum(eig))

# plot PCoA (color by sample category, shape for MRSA strains)
ggplot(data = ord, aes(x = pcoa1, y = pcoa2, fill = Category, shape = Category)) +
  geom_point(alpha = 0.65, 
             #shape = 21,
             stroke = 0.5, size = 5) +
  theme_bw() +
  xlab('PCo1 (26.4%)') +
  ylab('PCo2 (15.4%)') +
  theme(#axis.text = element_text(size=15, color = "black"),
        axis.title = element_text(size=17, color = "black"),
        legend.text = element_text(size=15, color = "black"),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.title = element_blank()) +
  scale_shape_manual(values=c(21,21,21,24,21),
                     guide = "none") +
  scale_fill_manual(values=c("gray", "skyblue", "yellow", "yellow", "brown"),
                    labels=c("BE-E (n=3)", "BE-SC (n=21)", "Human (n=24)", "H-MRSA (n=55)", "Soil (n=2)")) +
  guides(fill = guide_legend(override.aes=list(shape=c(21,21,21,24,21)))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

### Compute avg. "distance" between human, MRSA, and ISS strains ###

dist_test = data.frame(Human = apply(pangenome_staph[,grep("Human",staph_meta$Category)],1,mean),
                       MRSA = apply(pangenome_staph[,grep("MRSA",staph_meta$Category)],1,mean),
                       ISS = apply(pangenome_staph[,grep("BE_ISS",staph_meta$Category)],1,mean),
                       BE_soil = apply(pangenome_staph[,c(grep("BE_Earth",staph_meta$Category),
                                                          grep("Soil",staph_meta$Category))],1,mean),
                       BE_E = apply(pangenome_staph[,grep("BE_Earth",staph_meta$Category)],1,mean),
                       Soil = apply(pangenome_staph[,grep("Soil",staph_meta$Category)],1,mean))

vegdist(t(dist_test), 'jaccard', binary = T)

### PERMANOVA (may need to re-compute distance due to "NA's") ###

# Strain origin
beta <- vegdist(t(pangenome_staph), 'jaccard', binary = T)
adonis(beta ~ Category, data = staph_meta, permutations = 999)

# Assembler
adonis(beta ~ Assembler2, data = staph_meta, permutations = 999) 

# Sequencer
adonis(beta ~ Seq_Tech2, data = staph_meta, permutations = 999) 

# Reference
adonis(beta ~ Reference, data = staph_meta, permutations = 999) 

# Culture Media
beta <- vegdist(t(pangenome_staph[,-c(which(is.na(staph_meta$Culture_Media) == T))]), 'jaccard', binary = T)
adonis(beta ~ Culture_Media, data = staph_meta, permutations = 999) 

# MRSA (Human vs. Human only)
beta <- vegdist(t(pangenome_staph[,c(grep("Human", staph_meta$Category), 
                                     grep("MRSA", staph_meta$Category))]), 'jaccard', binary = T)
adonis(beta ~ MRSA, data = staph_meta[c(grep("Human", staph_meta$Category), 
                                         grep("MRSA", staph_meta$Category)),], permutations = 999) 

###
### ISS ONLY: Genes ###
###

# beta-diversity measure
beta <- vegdist(t(pangenome_staph[,grep("ISS", staph_meta$Category)]), 'jaccard', binary = T)
beta.mat <- as.matrix(beta)

# projection
pcoa <- cmdscale(beta, k = 4, eig = TRUE)

# cleanup
ord <- as.data.frame(pcoa$points)
names(ord) <- c('pcoa1', 'pcoa2', 'pco3', 'pco4') ## rename coordinates

# add metadata
ord$Site = staph_meta[grep("ISS", staph_meta$Category),]$Sample_Site_2

# Percent explained variation
eig <- eigenvals(pcoa)
100*head(eig/sum(eig))

## plot PCoA (color by sample category)
ggplot(data = ord, aes(x = pcoa1, y = pcoa2, fill = Site, shape = Site)) +
  geom_point(alpha = 0.7, 
             shape = 24, 
             stroke = 0.9, size = 7) +
  theme_bw() +
  xlab('PCo1 (57.3%)') +
  ylab('PCo2 (12.3%)') +
  scale_fill_manual(values=c("#00FFCC", "skyblue", "#0000CC", "purple", "pink", "red", 
                             "#663300", "#FF6600", "#FFFF99", "#CCFF00", "green", "#003300")) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_text(size=20, color = "black"),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank()) + 
  #theme(legend.position = "none") +
  theme(legend.text = element_text(size=16, color = "black"),
         legend.title = element_blank(),
         legend.background = element_rect(size = 0.5, colour = 1),
         legend.position = "bottom") +
  guides(fill=guide_legend(nrow = 6, byrow=F)) 

###
### Wallace and Voorhies: Dataset Comparisons ###
###

# assign data
WV_pan = pangenome_staph[,grep("Wallace and Voorhies", staph_meta$Reference)]
WV_meta = staph_meta[grep("Wallace and Voorhies", staph_meta$Reference),]

# BE_Earth vs. BE_ISS (R2=0.134, p=0.177)
beta <- vegdist(t(WV_pan[,grep("BE", WV_meta$Category)]), 
                'jaccard', binary = T)
adonis(beta ~ Category, data = WV_meta[grep("BE", WV_meta$Category),], 
       permutations = 999)

# BE_Earth vs. Human (R2=0.089, p=0.972)
beta <- vegdist(t(WV_pan[,grep("BE_Earth|Human", WV_meta$Category)]), 
                'jaccard', binary = T)
adonis(beta ~ Category, data = WV_meta[grep("BE_Earth|Human", WV_meta$Category),], 
       permutations = 999)

# BE_ISS vs. Human (R2=0.139, p=0.097)
beta <- vegdist(t(WV_pan[,grep("BE_ISS|Human", WV_meta$Category)]), 
                'jaccard', binary = T)
adonis(beta ~ Category, data = WV_meta[grep("BE_ISS|Human", WV_meta$Category),], 
       permutations = 999)

###
### FUNCTIONS ###
###

### Prep data ###

## load "annotated" pangenome matrix
# **prepared from "Gene, Non-unique Gene name, Annotation" columns in gene_presence_absence.csv file
pangenome_anno = read.table("pangenomes_ad-hoc/R_files/S_aureus/pangenome_gene_annotations.txt", 
                            sep="\t",header=T,row.names=1, quote="\"")

## bin annotations for overall presence/absence
fxn_df = as.data.frame(matrix(nrow = length(levels(pangenome_anno$Annotation)), ncol = ncol(pangenome_staph)))
for (i in 1:ncol(pangenome_staph)) {
  fxn_df[,i] = tapply(pangenome_staph[,i], pangenome_anno$Annotation,sum)
}
rownames(fxn_df) = names(tapply(pangenome_staph[,1], pangenome_anno$Annotation,sum))
colnames(fxn_df) = colnames(pangenome_staph)

## convert fxn_df to binary presence/absence data
fxn_df_bin = ifelse(fxn_df>0, 1, 0)

## set data frame for S. aureus fxns (used later)
fxn_df_bin_staph = fxn_df_bin

### PCoA ###

# beta-diversity measure
beta <- vegdist(t(fxn_df_bin), 'jaccard', binary = T)
beta.mat <- as.matrix(beta)

# projection
pcoa <- cmdscale(beta, k = 4, eig = TRUE)

# cleanup
ord <- as.data.frame(pcoa$points)
names(ord) <- c('pcoa1', 'pcoa2', 'pco3', 'pco4') ## rename coordinates

# add metadata
ord$Category = staph_meta$Category
ord$MRSA = staph_meta$MRSA
ord$Assembler = staph_meta$Assembler2
ord$Seq_Tech = staph_meta$Seq_Tech2
ord$Reference = staph_meta$Reference
ord$Location = staph_meta$Sample_Site_2
ord$Culture_Media = staph_meta$Culture_Media

# Percent explained variation
eig <- eigenvals(pcoa)
100*head(eig/sum(eig))

### plot PCoA (color by sample category)
ggplot(data = ord, aes(x = pcoa1, y = pcoa2, fill = Category, shape = Category)) +
  geom_point(alpha = 0.65, 
             stroke = 0.5, size = 5) +
  theme_bw() +
  xlab('PCo1 (35.3%)') +
  ylab('PCo2 (16.4%)') +
  theme(axis.text = element_text(size=15, color = "black"),
        axis.title = element_text(size=16, color = "black"),
        legend.text = element_text(size=15, color = "black"),
        legend.title = element_blank()) +
  scale_shape_manual(values=c(21,21,21,24,21),
                     guide = "none") +
  scale_fill_manual(values=c("gray", "skyblue", "yellow", "yellow", "brown"),
                    labels=c("BE-E (n=3)", "BE-SC (n=21)", "Human (n=24)", "H-MRSA (n=55)", "Soil (n=2)")) +
  guides(fill = guide_legend(override.aes=list(shape=c(21,21,21,24,21)))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

### PERMANOVA (may need to re-compute distance due to "NA's") ###

# Strain origin
beta <- vegdist(t(fxn_df_bin), 'jaccard', binary = T)
adonis(beta ~ Category, data = staph_meta, permutations = 999)

# Assembler
adonis(beta ~ Assembler2, data = staph_meta, permutations = 999) 

# Sequencer
adonis(beta ~ Seq_Tech2, data = staph_meta, permutations = 999) 

# Reference
adonis(beta ~ Reference, data = staph_meta, permutations = 999) 

# Culture Media
beta <- vegdist(t(fxn_df_bin[,-c(which(is.na(staph_meta$Culture_Media) == T))]), 'jaccard', binary = T)
adonis(beta ~ Culture_Media, data = staph_meta, permutations = 999) 

# MRSA (Human vs. Human only)
beta <- vegdist(t(fxn_df_bin[,c(grep("Human", staph_meta$Category),
                                grep("MRSA", staph_meta$Category))]), 'jaccard', binary = T)
adonis(beta ~ MRSA, data = staph_meta[c(grep("Human", staph_meta$Category),
                                        grep("MRSA", staph_meta$Category)),], permutations = 999) 

###
### ISS ONLY: Functions ###
###

# beta-diversity measure
beta <- vegdist(t(fxn_df_bin[,grep("ISS", staph_meta$Category)]), 'jaccard', binary = T)
beta.mat <- as.matrix(beta)

# subset metadata table
staph_meta_ISS = staph_meta[grep("ISS", staph_meta$Category),]

# Reference
adonis(beta ~ Reference, data = staph_meta_ISS, permutations = 999) 

# Year
staph_meta_ISS$Year = substr(staph_meta[grep("ISS", staph_meta$Category),]$Sample_Site, 
                             start = nchar(as.character(staph_meta[grep("ISS", staph_meta$Category),]$Sample_Site))-4, 
                             stop = nchar(as.character(staph_meta[grep("ISS", staph_meta$Category),]$Sample_Site))-1)
adonis(beta ~ Year, data = staph_meta_ISS, permutations = 999) 


####################################
### Trees -- core gene alignment ###
####################################

### B. cereus ###

## load core gene alignment
core_tree_bac <- read.tree("pangenome_roary_outputs/roary_out_B_cereus/core_gene_aln_tree.newick")

## set tip fill colors (by strain origin)
col.tree = c(1:56)
col.tree[grep("Human", bac_meta$Category)] = rep("yellow", length(grep("Human", bac_meta$Category)))
col.tree[grep("ISS", bac_meta$Category)] = rep("skyblue", length(grep("ISS", bac_meta$Category)))
col.tree[grep("Soil", bac_meta$Category)] = rep("brown", length(grep("Soil", bac_meta$Category)))
col.tree[grep("BE_Earth", bac_meta$Category)] = rep("gray", length(grep("BE_Earth", bac_meta$Category)))
col.tree[grep("Cul_Earth", bac_meta$Category)] = rep("green", length(grep("Cul_Earth", bac_meta$Category)))
col.tree[grep("Cul_SC", bac_meta$Category)] = rep("darkgreen", length(grep("Cul_SC", bac_meta$Category)))
bac_meta$Color = col.tree
for (i in 1:56){
  col.tree[i] = bac_meta$Color[grep(core_tree_bac$tip.label[i], rownames(bac_meta))]
}

## set labels
label.tree = c(1:56)
for (i in 1:56){
  label.tree[i] = as.character(bac_meta$Sample_Site_2)[grep(core_tree_bac$tip.label[i], rownames(bac_meta))]
}
label.tree = paste0("   ", label.tree)
core_tree_bac$tip.label = label.tree

## plot phylogram
plot(midpoint(core_tree_bac),
     cex = 0.9,
     font = 1,
     no.margin = F,
     direction = "downwards",
     edge.width = 1.2)
tiplabels(bg = col.tree, pch = 21, adj = c(0.501, 0.496), cex = 1.3)
add.scale.bar(ask = TRUE, cex = 1.2)

## correlation with accessory gene diversity

# re-load core gene alignment (for matching column/row names)
core_tree_bac <- read.tree("pangenome_roary_outputs/roary_out_B_cereus/core_gene_aln_tree.newick")

# prep core gene distance matrix
bac_core_distance <- cophenetic.phylo(midpoint(core_tree_bac))
bac_core_distance <- bac_core_distance[order(rownames(bac_core_distance)), order(colnames(bac_core_distance))]

# prep accessory gene distance matrix
bac_accessory_distance <- vegdist(t(pangenome_bac[which(apply(pangenome_bac, 1, sum)<0.95*dim(pangenome_bac)[2]),]), 
                                  'jaccard', binary = T)
bac_accessory_distance <- as.matrix(bac_accessory_distance)
bac_accessory_distance <- bac_accessory_distance[order(rownames(bac_accessory_distance)), order(colnames(bac_accessory_distance))]

# check for same order
colnames(bac_core_distance) == colnames(bac_accessory_distance)
rownames(bac_core_distance) == rownames(bac_accessory_distance)

# run mantel test
mantel(bac_core_distance, bac_accessory_distance)

### S. aureus ###

## load core gene alignment
core_tree_staph <- read.tree("pangenome_roary_outputs/roary_out_S_aureus/core_gene_aln_tree.newick")

## set tip fill colors (by strain origin)
col.tree = c(1:105)
col.tree[grep("Human", staph_meta$Category)] = rep("yellow", length(grep("Human", staph_meta$Category)))
col.tree[grep("MRSA", staph_meta$Category)] = rep("pink", length(grep("MRSA", staph_meta$Category)))
col.tree[grep("ISS", staph_meta$Category)] = rep("skyblue", length(grep("ISS", staph_meta$Category)))
col.tree[grep("Soil", staph_meta$Category)] = rep("brown", length(grep("Soil", staph_meta$Category)))
col.tree[grep("BE_Earth", staph_meta$Category)] = rep("gray", length(grep("BE_Earth", staph_meta$Category)))
staph_meta$Color = col.tree
for (i in 1:105){
  col.tree[i] = staph_meta$Color[grep(core_tree_staph$tip.label[i], rownames(staph_meta))]
}

## add shape-triangles for MRSA strains
shape.tree = rep(21, 105)
shape.tree[grep("pink", col.tree)] = rep(24, length(grep("pink", col.tree)))
col.tree[grep("pink", col.tree)] = rep("yellow", length(grep("pink", col.tree)))

## set labels
label.tree = c(1:105)
for (i in 1:105){
  label.tree[i] = as.character(staph_meta$Sample_Site_2)[grep(core_tree_staph$tip.label[i], rownames(staph_meta))]
}
label.tree = paste0("    ", label.tree)
core_tree_staph$tip.label = label.tree

#core_tree_staph$tip.label = paste0("    ",substr(core_tree_staph$tip.label, start = 1, stop = nchar(core_tree_staph$tip.label) - 5))

## plot phylogram
plot(midpoint(core_tree_staph), 
     #type = "unrooted",
     cex = 0.6,
     font = 1,
     no.margin = F,
     direction = "downward",
     edge.width = 1.2)
tiplabels(bg = col.tree, pch = shape.tree, adj = c(0.501, 0.4995), cex = 0.9)
add.scale.bar(ask = TRUE, cex = 1.2)

## correlation with accessory gene diversity

# re-load core gene alignment (for matching column/row names)
core_tree_staph <- read.tree("pangenome_roary_outputs/roary_out_S_aureus/core_gene_aln_tree.newick")

# prep core gene distance matrix
staph_core_distance <- cophenetic.phylo(midpoint(core_tree_staph))
staph_core_distance <- staph_core_distance[order(rownames(staph_core_distance)), order(colnames(staph_core_distance))]

# prep accessory gene distance matrix
staph_accessory_distance <- vegdist(t(pangenome_staph[which(apply(pangenome_staph, 1, sum)<0.95*dim(pangenome_staph)[2]),]), 
                                    'jaccard', binary = T)
staph_accessory_distance <- as.matrix(staph_accessory_distance)
staph_accessory_distance <- staph_accessory_distance[order(rownames(staph_accessory_distance)), order(colnames(staph_accessory_distance))]

# check for same order
colnames(staph_core_distance) == colnames(staph_accessory_distance)
rownames(staph_core_distance) == rownames(staph_accessory_distance)

# run mantel test
mantel(staph_core_distance, staph_accessory_distance)


#################################
### Function Enrichments: GLM ###
#################################

### B. cereus ###

## number of strains per origin encoding each function
fxn_sample_type_bac = as.data.frame(matrix(nrow = dim(fxn_df_bin_bac)[1], 
                                           ncol = length(levels(bac_meta$Category))))
for (i in 1:dim(fxn_df_bin_bac)[1]) {
  fxn_sample_type_bac[i,] = tapply(fxn_df_bin_bac[i,], bac_meta$Category, sum)
}
rownames(fxn_sample_type_bac) = rownames(fxn_df_bin_bac)
colnames(fxn_sample_type_bac) = levels(bac_meta$Category)

## proportions per origin
fxn_sample_type_bac$BE_Earth_prop = 100*fxn_sample_type_bac$BE_Earth/length(grep("BE_Earth", bac_meta$Category))
fxn_sample_type_bac$BE_ISS_prop = 100*fxn_sample_type_bac$BE_ISS/length(grep("BE_ISS", bac_meta$Category))
fxn_sample_type_bac$Soil_prop = 100*fxn_sample_type_bac$Soil/length(grep("Soil", bac_meta$Category))
fxn_sample_type_bac$Human_prop = 100*fxn_sample_type_bac$Human/length(grep("Human", bac_meta$Category))
fxn_sample_type_bac$Cul_Earth_prop = 100*fxn_sample_type_bac$Cul_Earth/length(grep("Cul_Earth", bac_meta$Category))
fxn_sample_type_bac$Cul_SC_prop = 100*fxn_sample_type_bac$Cul_SC/length(grep("Cul_SC", bac_meta$Category))
head(fxn_sample_type_bac)

## create vector for samples to omit for statistical analysis (i.e., those with n<3)
omit = c(which(bac_meta$Category == "Culture"), which(bac_meta$Category == "BE_Earth"))

## GLM p-values for function presence/absence by origin
p.glm.fxn_bac = c(1:dim(fxn_df_bin_bac)[1])
for (i in 1:dim(fxn_df_bin_bac)[1]) {
  p.glm.fxn_bac[i] = unlist(anova(glm(as.numeric(fxn_df_bin_bac[i,-c(omit)]) ~ bac_meta$Category[-c(omit)], 
                                      family = binomial(link='logit')), test = "Chisq"))[10]
} 

## set table for p-values and corrected p-values
fxn_assoc_bac = data.frame(fxn = rownames(fxn_df_bin_bac),
                           p.val = p.glm.fxn_bac,
                           p.adj.BH = p.adjust(p.glm.fxn_bac, method = "BH"))

## combine function totals and p-values
fxn_all_bac = data.frame(fxn_sample_type_bac, fxn_assoc_bac)

## re-order based on significance
fxn_all_bac = fxn_all_bac[order(fxn_all_bac$p.adj.BH),]

## write table
#write.table(fxn_all_bac, 'pangenomes_ad-hoc/R_files/B_cereus/fxns_proportions_glm.txt', sep="\t", col.names=NA, quote = FALSE)

### S. aureus ###

## number of strains per origin encoding each function
fxn_sample_type_staph = as.data.frame(matrix(nrow = dim(fxn_df_bin_staph)[1], 
                                             ncol = length(levels(staph_meta$Category))))
for (i in 1:dim(fxn_df_bin_staph)[1]) {
  fxn_sample_type_staph[i,] = tapply(fxn_df_bin_staph[i,], staph_meta$Category, sum)
}
rownames(fxn_sample_type_staph) = rownames(fxn_df_bin_staph)
colnames(fxn_sample_type_staph) = levels(staph_meta$Category)

## proportions per origin
fxn_sample_type_staph$BE_Earth_prop = 100*fxn_sample_type_staph$BE_Earth/length(grep("BE_Earth", staph_meta$Category))
fxn_sample_type_staph$BE_ISS_prop = 100*fxn_sample_type_staph$BE_ISS/length(grep("BE_ISS", staph_meta$Category))
fxn_sample_type_staph$Soil_prop = 100*fxn_sample_type_staph$Soil/length(grep("Soil", staph_meta$Category))
fxn_sample_type_staph$Human_prop = 100*fxn_sample_type_staph$Human/length(grep("Human", staph_meta$Category))
fxn_sample_type_staph$MRSA_prop = 100*fxn_sample_type_staph$MRSA/length(grep("MRSA", staph_meta$Category))
head(fxn_sample_type_staph)

## create vector for samples to omit for statistical analysis (i.e., those with n<3)
omit = c(which(staph_meta$Category == "Soil"))

## compute glm p-values for presence-absence of functions by sample category
p.glm.fxn_staph = c(1:dim(fxn_df_bin_staph)[1])
for (i in 1:dim(fxn_df_bin_staph)[1]) {
  p.glm.fxn_staph[i] = unlist(anova(glm(as.numeric(fxn_df_bin_staph[i,-c(omit)]) ~ staph_meta$Category[-c(omit)], 
                                        family = binomial(link='logit')), test = "Chisq"))[10]
} 

## set table for p-values and corrected p-values
fxn_assoc_staph = data.frame(fxn = rownames(fxn_df_bin_staph),
                             p.val = p.glm.fxn_staph,
                             p.adj.BH = p.adjust(p.glm.fxn_staph, method = "BH"))

## combine function totals and p-values
fxn_all_staph = data.frame(fxn_sample_type_staph, fxn_assoc_staph)

## re-order based on significance
fxn_all_staph = fxn_all_staph[order(fxn_all_staph$p.adj.BH),]

## write table
#write.table(fxn_all_staph, 'pangenomes_ad-hoc/R_files/S_aureus/fxns_proportions_glm.txt', sep="\t", col.names=NA, quote = FALSE)


##################################
### Function Enrichments: ARGs ###
##################################

bac_sig_fxns = fxn_all_bac[which(fxn_all_bac$p.adj.BH < 0.1),c(8,9,10,15)]
bac_sig_fxns_res = bac_sig_fxns[grep("esistance|lactam|acrolide|transferase|etracyclin", rownames(bac_sig_fxns)),]
bac_sig_fxns_res # screen functions in Uniprot to confirm whether biological process is "antibiotic resistance"

staph_sig_fxns = fxn_all_staph[which(fxn_all_staph$p.adj.BH < 0.1),c(6,7,9,10,13)]
staph_sig_fxns_res = staph_sig_fxns[grep("esistance|lactam|acrolide|transferase|etracyclin", rownames(staph_sig_fxns)),]
staph_sig_fxns_res # screen functions in Uniprot to confirm whether biological process is "antibiotic resistance"


######################################
### Function Enrichments: Heatmaps ###
######################################

### B. cereus ###

## load function categories (from Uniprot data in Supp. Table 2)
bac_fxn_cat = read.table("pangenomes_ad-hoc/R_files/B_cereus/fxns_biological_process_group.txt", sep="\t",header=T,row.names=1, quote="\"")
bac_fxn_cat$Gene.Product == rownames(fxn_all_bac[which(fxn_all_bac$p.adj.BH < 0.001),])

## set row side colors by biological process
tapply(bac_fxn_cat$GO.Bio.2, bac_fxn_cat$GO.Bio.2, length)

bac_fxn_cat$colors = rep("beige", length(bac_fxn_cat$Gene.Product))
bac_fxn_cat$colors[grep("antibiotic_biosynthesis", bac_fxn_cat$GO.Bio.2)] = rep("blue", length(grep("antibiotic_biosynthesis", bac_fxn_cat$GO.Bio.2)))
#bac_fxn_cat$colors[grep("antibiotic_resistance", bac_fxn_cat$GO.Bio.2)] = rep("pink", length(grep("antibiotic_resistance", bac_fxn_cat$GO.Bio.2)))
bac_fxn_cat$colors[grep("biosynthesis_other", bac_fxn_cat$GO.Bio.2)] = rep("skyblue", length(grep("biosynthesis_other", bac_fxn_cat$GO.Bio.2)))
bac_fxn_cat$colors[grep("catabolism", bac_fxn_cat$GO.Bio.2)] = rep("darkgreen", length(grep("catabolism", bac_fxn_cat$GO.Bio.2)))
bac_fxn_cat$colors[grep("metabolism", bac_fxn_cat$GO.Bio.2)] = rep("brown", length(grep("metabolism", bac_fxn_cat$GO.Bio.2)))
bac_fxn_cat$colors[grep("stress", bac_fxn_cat$GO.Bio.2)] = rep("purple", length(grep("stress", bac_fxn_cat$GO.Bio.2)))
bac_fxn_cat$colors[grep("toxin-antitoxin", bac_fxn_cat$GO.Bio.2)] = rep("gray", length(grep("toxin-antitoxin", bac_fxn_cat$GO.Bio.2)))
bac_fxn_cat$colors[grep("transcription", bac_fxn_cat$GO.Bio.2)] = rep("yellow", length(grep("transcription", bac_fxn_cat$GO.Bio.2)))
bac_fxn_cat$colors[grep("transport", bac_fxn_cat$GO.Bio.2)] = rep("orange", length(grep("transport", bac_fxn_cat$GO.Bio.2)))
bac_fxn_cat$colors[grep("virulence", bac_fxn_cat$GO.Bio.2)] = rep("red", length(grep("virulence", bac_fxn_cat$GO.Bio.2)))

## set color palette
mypalette <- colorRampPalette(c("black","red","yellow"))(n = 100)
mycolors = c(0:100)

## plot key for heat colors
hist(c(0:100), breaks = 100, col = mypalette, 
     border = mypalette, ylim = c(0,0.5),
     xlab = NULL, ylab = NULL, 
     axes = FALSE,
     main = NA)
axis(side = 3, at = c(0,25,50,75,100), las = 2, srt = 45,
     labels = c("0%", "25%", "50%", "75%", "100%"))

## plot key for biological processes
plot(1:50, col = "white")
legend(1,49, legend=c("Other", "Antibiotic Biosynthesis", "Antibiotic Resistance", "Biosynthesis other",
                      "Catabolism", "DNA Repair", "Metabolism", "Stress Response", "Toxin-antitoxin",
                      "Transcription", "Transport", "Virulence"),
       c("beige", "blue", "pink", "skyblue", "darkgreen", "lightgreen",
         "brown", "purple", "gray", "yellow", "orange", "red"), 
       #pch=15, 
       cex=1)

## plot heatmap
dev.off()
heatmap.2(as.matrix(fxn_all_bac[which(fxn_all_bac$p.adj.BH < 0.001),
                                c(grep("BE_ISS_prop", colnames(fxn_all_bac)),
                                  grep("Soil_prop", colnames(fxn_all_bac)),
                                  grep("Human_prop", colnames(fxn_all_bac)))]),
          dendrogram = "row",
          Colv = FALSE,
          distfun = function(x) vegdist(x, method = 'jaccard'),
          key=FALSE, 
          symkey=FALSE, 
          density.info='none',
          trace = "none",
          cexCol = 1,
          labCol = F,
          cexRow = 1.1,
          col = mypalette,
          breaks = mycolors,
          lhei = c(0.5,20),
          lwid = c(0.8,5),
          offsetRow = 0,
          offsetCol = 0,
          RowSideColors = as.character(bac_fxn_cat$colors),
          margins = c(1, 35))

### S. aureus ###

## load function categories (from Uniprot data in Supp. Table 2)
staph_fxn_cat = read.table("pangenomes_ad-hoc/R_files/S_aureus/fxns_biological_process_group.txt", sep="\t",header=T,row.names=1, quote="\"")
staph_fxn_cat$Gene.Product == rownames(fxn_all_staph[which(fxn_all_staph$p.adj.BH < 0.001),])

## set row side colors by biological process
tapply(staph_fxn_cat$GO.Bio.2, staph_fxn_cat$GO.Bio.2, length)

staph_fxn_cat$colors = rep("beige", length(staph_fxn_cat$Gene.Product))
staph_fxn_cat$colors[grep("antibiotic_biosynthesis", staph_fxn_cat$GO.Bio.2)] = rep("blue", length(grep("antibiotic_biosynthesis", staph_fxn_cat$GO.Bio.2)))
staph_fxn_cat$colors[grep("antibiotic_resistance", staph_fxn_cat$GO.Bio.2)] = rep("pink", length(grep("antibiotic_resistance", staph_fxn_cat$GO.Bio.2)))
staph_fxn_cat$colors[grep("biosynthesis_other", staph_fxn_cat$GO.Bio.2)] = rep("skyblue", length(grep("biosynthesis_other", staph_fxn_cat$GO.Bio.2)))
staph_fxn_cat$colors[grep("catabolism", staph_fxn_cat$GO.Bio.2)] = rep("darkgreen", length(grep("catabolism", staph_fxn_cat$GO.Bio.2)))
staph_fxn_cat$colors[grep("DNA_repair", staph_fxn_cat$GO.Bio.2)] = rep("lightgreen", length(grep("DNA_repair", staph_fxn_cat$GO.Bio.2)))
staph_fxn_cat$colors[grep("metabolism", staph_fxn_cat$GO.Bio.2)] = rep("brown", length(grep("metabolism", staph_fxn_cat$GO.Bio.2)))
staph_fxn_cat$colors[grep("stress", staph_fxn_cat$GO.Bio.2)] = rep("purple", length(grep("stress", staph_fxn_cat$GO.Bio.2)))
staph_fxn_cat$colors[grep("toxin-antitoxin", staph_fxn_cat$GO.Bio.2)] = rep("gray", length(grep("toxin-antitoxin", staph_fxn_cat$GO.Bio.2)))
staph_fxn_cat$colors[grep("transcription", staph_fxn_cat$GO.Bio.2)] = rep("yellow", length(grep("transcription", staph_fxn_cat$GO.Bio.2)))
staph_fxn_cat$colors[grep("transport", staph_fxn_cat$GO.Bio.2)] = rep("orange", length(grep("transport", staph_fxn_cat$GO.Bio.2)))
staph_fxn_cat$colors[grep("virulence", staph_fxn_cat$GO.Bio.2)] = rep("red", length(grep("virulence", staph_fxn_cat$GO.Bio.2)))

## plot heatmap
dev.off()
heatmap.2(as.matrix(fxn_all_staph[which(fxn_all_staph$p.adj.BH < 0.001),
                                  c(grep("BE_ISS_prop", colnames(fxn_all_staph)),
                                    grep("BE_Earth_prop", colnames(fxn_all_staph)),
                                    #grep("Soil_prop", colnames(fxn_all_staph)),
                                    grep("Human_prop", colnames(fxn_all_staph)),
                                    grep("MRSA_prop", colnames(fxn_all_staph)))]),
          dendrogram = "row",
          Colv = FALSE,
          distfun = function(x) vegdist(x, method = 'jaccard'),
          key=FALSE, 
          symkey=FALSE, 
          density.info='none',
          trace = "none",
          cexCol = 0.9,
          labCol = F,
          cexRow = 0.8,
          col = mypalette,
          breaks = mycolors,
          lhei = c(0.5,20),
          lwid = c(0.5,3),
          offsetRow = 0,
          RowSideColors = as.character(staph_fxn_cat$colors),
          offsetCol = 0,
          margins = c(1, 28))

### S. aureus: Space Only ###

# set matrix
ISS_staph_heat = fxn_df_bin_staph[which(data.frame(fxn_sample_type_staph,
                                                   fxn_assoc_staph)$p.adj.BH < 0.001),
                                  grep("ISS", staph_meta$Category)]
colnames(ISS_staph_heat) = substr(colnames(ISS_staph_heat), start = 1, stop = 13)

# prerpare metadata layers: study (Sielfaff, Wallace) and year (2002, 2008, 2010, 2011, 2015)
staph_ISS = cbind(1:21, 1:21)

staph_meta[grep("ISS", staph_meta$Category),]$Sample_Site
staph_ISS[grep("2002", staph_meta[grep("ISS", staph_meta$Category),]$Sample_Site), 1] = 
  rep("orange",length(grep("2002", staph_meta[grep("ISS", staph_meta$Category),]$Sample_Site)))
staph_ISS[grep("2008", staph_meta[grep("ISS", staph_meta$Category),]$Sample_Site), 1] = 
  rep("green",length(grep("2008", staph_meta[grep("ISS", staph_meta$Category),]$Sample_Site)))
staph_ISS[grep("2010", staph_meta[grep("ISS", staph_meta$Category),]$Sample_Site), 1] = 
  rep("brown",length(grep("2010", staph_meta[grep("ISS", staph_meta$Category),]$Sample_Site)))
staph_ISS[grep("2011", staph_meta[grep("ISS", staph_meta$Category),]$Sample_Site), 1] = 
  rep("purple",length(grep("2011", staph_meta[grep("ISS", staph_meta$Category),]$Sample_Site)))
staph_ISS[grep("2015", staph_meta[grep("ISS", staph_meta$Category),]$Sample_Site), 1] = 
  rep("pink",length(grep("2015", staph_meta[grep("ISS", staph_meta$Category),]$Sample_Site)))

staph_meta[grep("ISS", staph_meta$Category),]$Reference
staph_ISS[grep("Sielaff", staph_meta[grep("ISS", staph_meta$Category),]$Reference), 2] = 
  rep("red",length(grep("Sielaff", staph_meta[grep("ISS", staph_meta$Category),]$Reference)))
staph_ISS[grep("Wallace", staph_meta[grep("ISS", staph_meta$Category),]$Reference), 2] = 
  rep("blue",length(grep("Wallace", staph_meta[grep("ISS", staph_meta$Category),]$Reference)))

# plot heatmap
heatmap.plus(as.matrix(ISS_staph_heat),
          cexCol = 0.8,
          cexRow = 0.6,
          scale = "none",
          col = c("black", "yellow"),
          ColSideColors = staph_ISS,
          margins = c(9, 17))


##########################################################
### KEGG pathway correlation with phylogeny: S. aureus ###
##########################################################

## load data and subset pathways with differential ontology counts
kegg_s_aureus_subset = read.table("kegg_functions_subset/s_aureus_7branches_clevel.txt",row.names=1, sep="\t",header=T)

## subset pathways with different ontology counts
kegg_s_aureus_subset_clean = kegg_s_aureus_subset[which(apply(kegg_s_aureus_subset,1,sd) > 0),]

## PCoA

# beta-diversity measure
beta <- vegdist(t(kegg_s_aureus_subset_clean))
beta.mat <- as.matrix(beta)

# projection
pcoa <- cmdscale(beta, k = 4, eig = TRUE)

# cleanup
ord <- as.data.frame(pcoa$points)
names(ord) <- c('pco1', 'pco2', 'pco3', 'pco4') ## rename coordinates

# add metadata
ord$Origin = substr(colnames(kegg_s_aureus_subset_clean), start = 10, stop = nchar(colnames(kegg_s_aureus_subset_clean)))
ord$Branch = substr(colnames(kegg_s_aureus_subset_clean), start = 1, stop = 8)

# percent explained variation
eig <- eigenvals(pcoa)
100*head(eig/sum(eig))

# plot 
ggplot(data = ord, aes(x = pco1, y = pco2, col = Branch, shape = Origin)) +
  geom_point(alpha = 0.9, 
             stroke = 0.5, size = 6.5) +
  theme_bw() +
  xlab('PCo1 (24.8%)') +
  ylab('PCo2 (20.4%)') +
  scale_color_manual(values = c("blue", "red", "brown", "orange", "purple", "green", "pink")) +
  theme(axis.text = element_text(size=14, color = "black"),
        axis.title = element_text(size=15, color = "black"),
        legend.text = element_text(size=15, color = "black"),
        legend.title = element_text(size = 16))

##PERMANOVA

# tree branch
adonis(beta ~ Branch, data = ord, permutations = 999)

# sample origin
adonis(beta ~ Origin, data = ord, permutations = 999)


######################################
### Functional overlaps cross taxa ###
######################################

## Venn diagrams: pan, core, accessory
venn.diagram(list(rownames(fxn_df_bin_bac), 
                  rownames(fxn_df_bin_staph)), 
             category.names = c("B. cereus", "S. aureus"),
             filename = 'pangenomes_ad-hoc/R_files/taxon_overlaps/fxn_overlap_all.png',
             fill = c('blue', 'red'),
             cex = 4,
             scaled = 0,
             cat.cex = 0,
             cat.pos = 0) 
venn.diagram(list(rownames(fxn_df_bin_bac[which(apply(fxn_df_bin_bac,1,sum) < dim(fxn_df_bin_bac)[2]*.95),]), 
                  rownames(fxn_df_bin_staph[which(apply(fxn_df_bin_staph,1,sum) < dim(fxn_df_bin_staph)[2]*.95),])), 
             category.names = c("B. cereus", "S. aureus"),
             filename = 'pangenomes_ad-hoc/R_files/taxon_overlaps/fxn_overlap_accessory.png',
             fill = c('blue', 'red'),
             cex = 4,
             scaled = 0,
             cat.cex = 0,
             cat.pos = 0) 
venn.diagram(list(rownames(fxn_df_bin_bac[which(apply(fxn_df_bin_bac,1,sum) > dim(fxn_df_bin_bac)[2]*.95),]), 
                  rownames(fxn_df_bin_staph[which(apply(fxn_df_bin_staph,1,sum) > dim(fxn_df_bin_staph)[2]*.95),])), 
             category.names = c("B. cereus", "S. aureus"),
             filename = 'pangenomes_ad-hoc/R_files/taxon_overlaps/fxn_overlap_core.png',
             fill = c('blue', 'red'),
             cex = 4,
             scaled = 0,
             cat.cex = 0,
             cat.pos = 0) 

## overlapping enrichments
overlap = intersect(as.character(fxn_assoc_staph[which(fxn_assoc_staph$p.adj.BH < 0.1),1]),
                    as.character(fxn_assoc_bac[which(fxn_assoc_bac$p.adj.BH < 0.1),1]))
overlap_call = data.frame(fxn_all_bac[grep(paste(overlap, collapse= "|"),rownames(fxn_all_bac)),
                                      c(grep("BE_ISS_prop", colnames(fxn_all_bac)),
                                        grep("Soil_prop", colnames(fxn_all_bac)),
                                        grep("Human_prop", colnames(fxn_all_bac)))],
                          fxn_all_staph[grep(paste(overlap, collapse= "|"),rownames(fxn_all_staph)),
                                        c(grep("BE_ISS_prop", colnames(fxn_all_staph)),
                                          grep("Soil_prop", colnames(fxn_all_staph)),
                                          grep("Human_prop", colnames(fxn_all_staph)),
                                          grep("MRSA_prop", colnames(fxn_all_staph)))])

# subset data differences between Human and ISS strains
lolly_overlap = data.frame(Human_vs_ISS = c((overlap_call[,3] - overlap_call[,1]),
                                            (overlap_call[,5] - overlap_call[,4])),
                           Species = c(rep("Bacillus", 7), rep("Staphylococcus", 7)),
                           Fxn = rep(rownames(overlap_call), 2))

# plot differene in percentage strains encoding each enriched function
ggplot(lolly_overlap) +
  geom_linerange(aes(x = Fxn, ymin = 0, ymax = Human_vs_ISS, color = Species),
                 size = 2, alpha = 0.7) +
  geom_point(aes(x = Fxn, y = Human_vs_ISS, color = Species), 
             size = 12, alpha = 0.7) +
  coord_flip() +
  theme_bw() +
  scale_color_manual(values = c("blue", "red"),
                     labels = c("B. cereus", "S. aureus")) +
  ylab("Difference in Percentage of Strains Encoding Function") +
  ylim(-100,100) +
  scale_x_discrete(position = "top") +
  geom_hline(yintercept = c(0)) +
  theme(#axis.title.x = element_text(size=20, color = "black"),
        axis.title.x = element_blank(), 
        #axis.text.x = element_text(size=20, color = "black"),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(), 
        #axis.text.y = element_text(size=12, color = "black"),
        axis.text.y = element_blank(),
        #axis.ticks.y = element_blank(),
        legend.position = "none",
        panel.border = element_rect(color = "black"))


##########################################################################
### Gene presence/absence distance: original vs. standardized assembly ###
##########################################################################

## load pan-genome of B. cereus de novo assemblies with spades
pangenome_bac_spades = read.table("pangenome_roary_outputs/roary_out_B_cereus_de_novo_assembly/gene_presence_absence.Rtab",
                                  h=T, row.names=1)

## load associated metadata
bac_meta_spades = read.table("pangenomes_ad-hoc/R_files/B_cereus_de_novo_assembly/Bac_metadata_spades.txt",
                             sep='\t',h=T,row.names=1,check=F,comment='')

## subset B. cereus pangenome (from all genomes) to appropriate counterparts from de novo assemblies
pangenome_bac_all_subset <- pangenome_bac[,grep(paste(bac_meta_spades$GCA_Counterpart, collapse = "|"), colnames(pangenome_bac))]
bac_all_subset_meta <- bac_meta[grep(paste(bac_meta_spades$GCA_Counterpart, collapse = "|"), colnames(pangenome_bac)),]

## subset B. cereus de novo assemblies pangenome to match counterparts from original assemblies
pangenome_bac_spades_subset <- pangenome_bac_spades[,grep(paste(colnames(pangenome_bac_all_subset), collapse = "|"), bac_meta_spades$GCA_Counterpart)]
bac_spades_subset_meta <- bac_meta_spades[grep(paste(colnames(pangenome_bac_all_subset), collapse = "|"), bac_meta_spades$GCA_Counterpart),]

## confirm matching sample names across datasets
bac_spades_subset_meta$GCA_Counterpart[order(bac_spades_subset_meta$GCA_Counterpart)] == 
  colnames(pangenome_bac_all_subset)[order(colnames(pangenome_bac_all_subset))]

## distance stats on de novo assemblies (mean, se)
mean(vegdist(t(pangenome_bac_spades_subset), 'jaccard', binary = T)) #avg
sd(vegdist(t(pangenome_bac_spades_subset), 'jaccard', binary = T))/sqrt(120) #se

## distance stats on original assemblies (mean, se)
mean(vegdist(t(pangenome_bac_all_subset), 'jaccard', binary = T))
sd(vegdist(t(pangenome_bac_all_subset), 'jaccard', binary = T))/sqrt(120)

## wilcox test 
wilcox.test(vegdist(t(pangenome_bac_spades_subset), 'jaccard', binary = T),
            vegdist(t(pangenome_bac_all_subset), 'jaccard', binary = T))



