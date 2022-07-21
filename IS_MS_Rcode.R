#### IS MS figure and stat R script ####

# The paper this R code was written for can be accessed at: 
# Greenrod STE, Stoycheva M, Elphinstone J, Friman V-P. 2022. Influence of insertion sequences on population structure of phytopathogenic bacteria in the Ralstonia solanacearum species complex. bioRxiv: 2022.07.16.500299.

### Dependencies ###

require(ggplot2)
require(ggbeeswarm)
require(kableExtra)
require(tidyr)
require(lme4)
require(ggfortify)
require(MASS)
require(reshape2)
require(gridExtra)
require(corrplot)
require(ggtree)
require(phytools)
require(tidyr)
require(pheatmap)
require(scico)
require(vegan)
require(phangorn)
require(paco)
require(factoextra)
require(FSA)
require(car)
require(treeio)
require(Polychrome)


### RSSC phylogeny generation

RSSC_tree <- read.tree("core_variant_tree_356subset.treefile")
ggtree(RSSC_tree) + geom_tiplab() + geom_text(aes(label=node), hjust=-.3)

RSSC_tree_ladderize <- ladderize(RSSC_tree, right=FALSE)
is_tip <- RSSC_tree_ladderize$edge[,2] <= length(RSSC_tree_ladderize$tip.label)
ordered_tips <- RSSC_tree_ladderize$edge[is_tip, 2]
ordered_RSSC_tree_tips <- data.frame(x$tip.label[ordered_tips])
write.csv(ordered_RSSC_tree_tips,"Whole_phylogeny_tiplabels_ref.csv")

RSSC_tree_fig <- ggtree(Isolate_tree,ladderize=TRUE)+ theme(plot.margin=margin(0,0,30,0))


# Figure S1 - Ralstonia tree + references

RSSC_tree_ref <- read.tree("Samuels_phylo_march_2022.treefile")
RSSC_tree_ref_rooted <- phytools::midpoint.root(RSSC_tree_ref)

RSSC_tree_ref_rooted_ladderize <- ladderize(RSSC_tree_ref_rooted, right=FALSE)
is_tip <- RSSC_tree_ref_rooted_ladderize$edge[,2] <= length(RSSC_tree_ref_rooted_ladderize$tip.label)
ordered_tips <- RSSC_tree_ref_rooted_ladderize$edge[is_tip, 2]
ordered_RSSC_tree_ref_rooted_tips <- data.frame(RSSC_tree_ref_rooted_ladderize$tip.label[ordered_tips])
write.csv(ordered_RSSC_tree_ref_rooted_tips,"Whole_phylogeny_tiplabels_ref.csv")

RSSC_tree_ref_fig <- ggtree(RSSC_tree_rooted,ladderize=TRUE)+ theme(plot.margin=margin(0,0,30,0))


### IS sorting and filtering

## Determine IS clusters

is_mash_matrix <- read.csv("All_Nanopore_IS_filtered_mashmatrix_heatmap.csv")
rownames(is_mash_matrix) <- is_mash_matrix[, 1];
is_mash_matrix2 <- is_mash_matrix[, -1];
is_mash_matrix3 <- as.matrix(is_mash_matrix2)
is_mash_heatmap <- pheatmap(is_mash_matrix3, kmeans_k = 66,show_rownames = FALSE, show_colnames = FALSE)

fviz_nbclust(is_mash_matrix3, FUN = kmeans,k.max = 80, method = "wss")
fviz_nbclust(is_mash_matrix3, FUN = kmeans,k.max = 80, method = "silhouette") # 66

kmean_IS_clusters <- kmeans(is_mash_matrix3, centers = 66, nstart = 25)

kmean_IS_clusters_labels <- data.frame(kmean_IS_clusters$cluster)
write.csv(kmean_IS_clusters_labels,"IS_clusters.csv")


## Figure S2A - Compare ISEScan and ISMapper IS prediction number

isescan_v_ismapper <- read.csv("isescan_v_ismapper.csv",fileEncoding="UTF-8-BOM")

ggplot(isescan_v_ismapper,aes(x=Num_IS_ISEScan,y=Num_IS_ISMapper)) + geom_point() +
  geom_abline(intercept = 0, slope = 1,linetype="dotted")+
  theme_bw()+
  theme(axis.title=element_text(size=10,face="bold"))+
  ylab("IS copy number (ISMapper)")+
  xlab("IS copy number (ISEScan)")

hist(IStotal$Num_IS_ISMapper)
hist(IStotal$Num_IS_ISEScan)

test <- cor.test(IStotal$Num_IS_ISEScan, IStotal$Num_IS_ISMapper, method = "kendall",use = "complete.obs")
test


## Figure S2B - Determine read coverage across RSSC phylogeny

Read_num <- read.csv("Isolate_read_depth.csv",fileEncoding="UTF-8-BOM")
Read_num$Isolate <- factor(Read_num$Isolate, levels = Read_num$Isolate)

safe_colorblind_palette <- c("#CC6677", "#DDCC77", "#117733","#6699CC", "#888888")

ggplot(data=Read_num, aes(x = Phylotype, y= Average_depth_perbase, fill = Region)) + geom_violin(alpha=0.5, position = position_dodge(width = .75),size=1,color=NA)+
  geom_boxplot(notch = FALSE,  outlier.size = -1, color="black",lwd=0.6, alpha = 0.7,show.legend = F, varwidth=TRUE,position = position_dodge(width = .75))+
  ggbeeswarm::geom_quasirandom(shape = 21,size=2, dodge.width = .75, color="black",alpha=.3,show.legend = F)+
  ylab("Average read depth per bp") +
  xlab("Phylotype") +
  theme_bw()+
  geom_smooth(method = "lm")+
  scale_fill_manual(values=c("#88CCEE","#0033FF"),labels = c("Chromosome", "Megaplasmid"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title=element_text(size=12,face="bold")) +
  theme(legend.title = element_text(colour="black", size=12,face="bold"),legend.text = element_text(size=12)) +
  labs(fill="Region")

leveneTest(log(Average_depth_perbase) ~ Phylotype * Region, data = Read_num)

kruskal.test(Average_depth_perbase ~ Phylotype, data = Read_num)
dunnTest(Average_depth_perbase ~ Phylotype, data = Read_num)


# Figure S3 - Comparing IS number between chromosome and megaplasmid

# Load IS data and gather variables
IStotal <- read.csv("IS_data.csv",fileEncoding="UTF-8-BOM")
IStotal_long <- gather(IStotal, IS_region, IS_freq, Chr_IS:MP_IS, factor_key=TRUE)
IStotal_long <- gather(IStotal_long, IS_prox_region, IS_prox_freq, Chr_prox_100bp:MP_prox_100bp, factor_key=TRUE)
IStotal_long <- gather(IStotal_long, IS_disrupt_region, IS_disrupt_freq, Chr_disruptions:MP_disruptions, factor_key=TRUE)
IStotal_long <- gather(IStotal_long, Prop_IS_disrupt_region, Prop_IS_disrupt_freq, Proportion_Chr_disruptions:Proportion_MP_disruptions, factor_key=TRUE)

# Plot chromosome v megaplasmid
ggplot(data=IStotal_long,aes(x=IS_region,y=IS_freq,fill=IS_region)) + geom_violin(alpha=0.5, position = position_dodge(width = .75),size=1,color=NA)+
  geom_boxplot(notch = FALSE,  outlier.size = -1, color="black",lwd=1, alpha = 0.3,show.legend = F, varwidth=TRUE,position = position_dodge(width = .75))+
  ggbeeswarm::geom_quasirandom(shape = 21,size=2, dodge.width = .75, color="black",alpha=.3,show.legend = F)+
  ylab("IS number per genome") +
  xlab("Genomic region") +
  theme_bw()+
  geom_smooth(method = "lm")+
  scale_fill_manual(values=c("#88CCEE","#0033FF"),labels = c("Chromosome", "Megaplasmid"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title=element_text(size=12,face="bold")) +
  labs(fill="Genomic Region") +
  scale_x_discrete(labels=c("Chromosome", "Megaplasmid"))

# Differences are normally distributed
d <- with(IStotal_long, 
          IS_freq[IS_region == "Chr_IS"] - IS_freq[IS_region =="MP_IS"])
hist(d)

# Stats to compare IS number between regions
t.test(IStotal$MP_IS , IStotal$Chr_IS, paired = TRUE, alternative = "two.sided")


## Figure 1 - Compare IS family abundance

ISfam <- read.csv("IS_family_abundance.csv",fileEncoding="UTF-8-BOM")
ISfam_long <- gather(ISfam, IS_family, IS_family_freq, IS110:ISL3, factor_key=TRUE)

Fig_1A <-ggplot(data=ISfam_long, aes(x = IS_family, y= IS_family_freq, fill = Region)) + geom_violin(alpha=0.5, position = position_dodge(width = .75),size=1,color=NA)+
  geom_boxplot(notch = FALSE,  outlier.size = -1, color="black",lwd=0.6, alpha = 0.7,show.legend = F, varwidth=TRUE,position = position_dodge(width = .75))+
  #ggbeeswarm::geom_quasirandom(shape = 21,size=2, dodge.width = .75, color="black",alpha=.1,show.legend = F)+
  ylab("IS number per genome") +
  xlab("IS family") +
  theme_bw()+
  geom_smooth(method = "lm")+
  scale_fill_manual(values=c("#88CCEE","#0033FF"),labels = c("Chromosome", "Megaplasmid"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title=element_text(size=12,face="bold")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(legend.title = element_text(colour="black", size=12,face="bold"),legend.text = element_text(size=12)) +
  labs(fill="Genomic region")


## Compare IS abundance between phylotypes

Fig_1B <- ggplot(data=IStotal_long, aes(x = Phylotype, y= IS_freq, fill = IS_region)) + geom_violin(alpha=0.5, position = position_dodge(width = .75),size=1,color=NA)+
  geom_boxplot(notch = FALSE,  outlier.size = -1, color="black",lwd=0.6, alpha = 0.7,show.legend = F, varwidth=TRUE,position = position_dodge(width = .75))+
  #ggbeeswarm::geom_quasirandom(shape = 21,size=2, dodge.width = .75, color="black",alpha=.3,show.legend = F)+
  ylab("IS number per genome") +
  xlab("Phylotype") +
  theme_bw()+
  geom_smooth(method = "lm")+
  scale_fill_manual(values=c("#88CCEE","#0033FF"),labels = c("Chromosome", "Megaplasmid"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title=element_text(size=12,face="bold")) +
  theme(legend.title = element_text(colour="black", size=12,face="bold"),legend.text = element_text(size=12)) +
  labs(fill="Genomic region")

kruskal.test(Num_IS ~ Phylotype, data = IStotal)
dunnTest(Num_IS ~ Phylotype, data = IStotal)

Fig_1A/Fig_1B + plot_layout(guides="collect")

## Compare IS family abundances between phylotypes

ISfam <- read.csv("IS_family_abundance_phylogeny.csv",fileEncoding="UTF-8-BOM")
ISfam_long <- gather(ISfam, IS_family, IS_family_freq, IS110:ISL3, factor_key=TRUE)

ISfam_long$Isolate <- as.character(ISfam_long$Isolate)
ISfam_long$Isolate <- factor(ISfam_long$Isolate, levels=unique(ISfam_long$Isolate))

Fig_1C <-ggplot(ISfam_long,aes(x=Isolate,y=IS_family_freq, fill=IS_family)) +
  geom_bar(stat="identity", position = "stack",width=1.1)+
  theme_bw()+
  theme(axis.ticks.y = element_blank()) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title=element_text(size=12,face="bold"),axis.text.y=element_blank(),axis.title.y=element_blank()) + 
  coord_flip()+ 
  scale_x_discrete(limits = rev(levels(ISfam_long$Isolate)))+
  scale_fill_viridis_d()+
  ylab("IS copy number") +
  theme(legend.title = element_text(colour="black", size=12,face="bold"),legend.text = element_text(size=12)) +
  theme(plot.margin=margin(0,0,0,20))+
  labs(fill="IS family")

grid.arrange(
  arrangeGrob(
    RSSC_tree_fig,Fig_1C,
    ncol=2 ,widths=c(0.4,1)))


## Figure 2A - Determine IS subfamily abundances across RSSC phylogeny

ISsubfam <- read.csv("IS_subfamily_abundance_phylogeny_total.csv",fileEncoding="UTF-8-BOM")
ISsubfam <- subset(ISsubfam, select = -c(Phylotype))

ISsubfam_matrix <- ISsubfam[, -1];
row.names(ISsubfam_matrix) <- ISsubfam$Isolate
ISsubfam_matrix2 <- as.matrix(ISsubfam_matrix)
ISsubfam_matrix2_melted <- melt(ISsubfam_matrix2)

Fig_2A <- ggplot(data = Chr_IS_melted_mat, aes(x=Var2, y=Var1, fill=value)) + 
  geom_tile()+
  scale_y_discrete(limits = rev)+
  theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.y=element_blank(),axis.text.x=element_blank(),legend.position="right")+
  scale_fill_gradient(low="grey98",high="black", limits=c(0,80)) + 
  theme(plot.margin=margin(0,0,0,20)) +
  theme(legend.title = element_text(colour="black", size=12,face="bold"),legend.text = element_text(size=12))+
  labs(fill="IS copy number")

grid.arrange(
  arrangeGrob(
    RSSC_tree_fig,Fig_2A,
    ncol=2 ,widths=c(0.4,1)))


## Determine IS subfamily abundances across RSSC phylogeny split between chromosome and megaplasmid (not included in MS)

# Chromosome
ISsubfam_Chr <- ISsubfam[- grep("Chr", ISsubfam$Region),]
ISsubfam_Chr <- subset(ISsubfam_Chr, select = -c(Region))

ISsubfam_Chr_matrix <- ISsubfam_Chr[, -1];
row.names(ISsubfam_Chr_matrix) <- ISsubfam_Chr$Isolate
ISsubfam_Chr_matrix2 <- as.matrix(ISsubfam_Chr_matrix)
ISsubfam_Chr_matrix2_melted <- melt(ISsubfam_Chr_matrix2)

ISsubfam_Chr_matrix_fig <- ggplot(data = Chr_IS_melted_mat, aes(x=Var2, y=Var1, fill=value)) + 
  geom_tile()+
  scale_y_discrete(limits = rev)+
  theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.y=element_blank(),axis.text.x=element_blank(),legend.position="right")+
  scale_fill_gradient(low="grey98",high="black", limits=c(0,90)) + 
  theme(plot.margin=margin(0,0,0,20)) +
  theme(legend.title = element_text(colour="black", size=12,face="bold"),legend.text = element_text(size=12))+
  labs(fill="IS copy number")

grid.arrange(
  arrangeGrob(
    RSSC_tree_fig,ISsubfam_Chr_matrix_fig,
    ncol=2 ,widths=c(0.4,1)))

# Megaplasmid

ISsubfam_MP <- ISsubfam[- grep("MP", ISsubfam$Region),]
ISsubfam_MP <- subset(ISsubfam_MP, select = -c(Region))

ISsubfam_MP_matrix <- ISsubfam_MP[, -1];
row.names(ISsubfam_MP_matrix) <- ISsubfam_MP$Isolate
ISsubfam_MP_matrix2 <- as.matrix(ISsubfam_MP_matrix2)
ISsubfam_MP_matrix2_melted <- melt(v)

ISsubfam_MP_matrix_fig <- ggplot(data = MP_IS_melted_mat, aes(x=Var2, y=Var1, fill=value)) + 
  geom_tile()+
  scale_y_discrete(limits = rev)+
  theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.y=element_blank(),axis.text.x=element_blank(),legend.position="right")+
  scale_fill_gradient(low="grey98",high="black", limits=c(0,45)) + 
  theme(plot.margin=margin(0,0,0,20)) +
  theme(legend.title = element_text(colour="black", size=12,face="bold"),legend.text = element_text(size=12))+
  labs(fill="IS copy number")

grid.arrange(
  arrangeGrob(
    RSSC_tree_fig,ISsubfam_MP_matrix_fig,
    ncol=2 ,widths=c(0.4,1)))


# Figure 2B - PCoA analysis of phylotype IS content

bray_curtis_ISsubfam_matrix <- vegdist(ISsubfam_matrix2, method="bray")

PcoA_phylotype <- pcoa(bray_curtis_ISsubfam_matrix)
PcoA_data <- data.frame(PcoA_phylotype$vectors)
Best_axes <- PcoA_data[,c("Axis.1","Axis.2")]
Best_axes$Phylotype <- matrix$Phylotype

ggplot(data=Best_axes, aes(x=Axis.1,y=Axis.2, colour = Phylotype))+ 
  geom_point(size=2)+
  theme_bw()+
  stat_ellipse(size=1)+
  labs(x="PCo1 (53.9%)",y="PCo2 (9.0%)")+
  scale_colour_manual(values=safe_colorblind_palette)+
  theme(axis.title=element_text(size=12,face="bold"))+
  theme(legend.title = element_text(colour="black", size=12,face="bold"),legend.text = element_text(size=12))+
  labs(fill="Phylotype")

ano = anosim(mat, matrix$Phylotype, distance = "bray", permutations = 9999)
ano

# Figure 2C - PACo analysis to see if IS content reflects host genetic distance

# Splitting RSSC phylogeny into clades

ggtree(RSSC_tree) + geom_text(aes(label=node), hjust=-.3)
phyloI <- tree_subset(RSSC_tree, node=643, levels_back=0)
ggtree(phyloI)+geom_tiplab()

phyloIIA <- tree_subset(RSSC_tree, node=627, levels_back=0)
ggtree(phyloIIA)+geom_tiplab()

phyloIIB <- tree_subset(RSSC_tree, node=359, levels_back=0)
ggtree(phyloIIB)+geom_tiplab()

phyloIII <- tree_subset(RSSC_tree, node=701, levels_back=0)
ggtree(phyloIII)+geom_tiplab()

phyloIV <- tree_subset(RSSC_tree, node=708, levels_back=0)
ggtree(phyloIV)+geom_tiplab()


## PACo analysis of phylotype I clade

phyloI_subfam_matrix <- read.csv("IS_subfamily_abundance_phylogeny_I.csv")
rownames(phyloI_subfam_matrix) <- phyloI_subfam_matrix[, 1];
phyloI_subfam_matrix2 <- phyloI_subfam_matrix[, -1];
phyloI_subfam_matrix3 <- as.matrix(phyloI_subfam_matrix2)

phyloI_subfam_BC_matrix <- vegdist(phyloI_subfam_matrix3, method="bray")
phyloI_subfam_BC_matrix2 <- as.matrix(phyloI_subfam_BC_matrix)
phyloI_subfam_BC_matrix2[is.na(phyloI_subfam_BC_matrix2)] <- 0
PhyloI_BC_UPGMA <- upgma(phyloI_subfam_BC_matrix2)

# Start running PACo
host.D <- cophenetic(phyloI)
BC.D <- cophenetic(PhyloI_BC_UPGMA)

HP <- read.csv("I_PACo_binarymatrix.csv",fileEncoding="UTF-8-BOM")
row.names(HP) <- HP$X
HP <- HP[, -1];
HP2 <- as.matrix(HP)
host.D <- host.D[rownames(HP2),colnames(HP2)]
BC.D <- BC.D[colnames(HP2),colnames(HP2)] 

D <- prepare_paco_data(H=host.D,P=BC.D,HP=HP2)
D <- add_pcoord(D,correction="cailliez")
D <- PACo(D,nperm=1000,seed=12,method="r0")

D <- paco_links(D)
res <- residuals_paco(D$proc)
D$gof

# Make tanglegram
assoc <- data.frame(pol=rownames(HP2)[which(HP2==1, arr.ind=TRUE)[,'row']], pla=colnames(HP2)[which(HP2==1, arr.ind=TRUE)[,'col']])

cophyloplot(I, PhyloI_UPGMA, assoc, show.tip.label=FALSE, use.edge.length=FALSE,
            lwd=1, col='steelblue', length.line=0, gap=-15, space=70,rotate=TRUE)


## PACo analysis of phylotype IIB clade

phyloIIB_subfam_matrix <- read.csv("IS_subfamily_abundance_phylogeny_IIB.csv")
rownames(phyloIIB_subfam_matrix) <- phyloIIB_subfam_matrix[, 1];
phyloIIB_subfam_matrix2 <- phyloIIB_subfam_matrix[, -1];
phyloIIB_subfam_matrix3 <- as.matrix(phyloIIB_subfam_matrix2)

phyloIIB_subfam_BC_matrix <- vegdist(phyloIIB_subfam_matrix3, method="bray")
phyloIIB_subfam_BC_matrix2 <- as.matrix(phyloIIB_subfam_BC_matrix)
phyloIIB_subfam_BC_matrix2[is.na(phyloIIB_subfam_BC_matrix2)] <- 0
PhyloIIB_BC_UPGMA <- upgma(phyloIIB_subfam_BC_matrix2)

# Start running PACo
host.D <- cophenetic(phyloIIB)
BC.D <- cophenetic(PhyloIIB_BC_UPGMA)

HP <- read.csv("IIB_PACo_binarymatrix.csv",fileEncoding="UTF-8-BOM")
row.names(HP) <- HP$X
HP <- HP[, -1];
HP2 <- as.matrix(HP)
host.D <- host.D[rownames(HP2),colnames(HP2)]
BC.D <- BC.D[colnames(HP2),colnames(HP2)] 

D <- prepare_paco_data(H=host.D,P=BC.D,HP=HP2)
D <- add_pcoord(D,correction="cailliez")
D <- PACo(D,nperm=1000,seed=12,method="r0")

D <- paco_links(D)
res <- residuals_paco(D$proc)
D$gof

# Make tanglegram
assoc <- data.frame(pol=rownames(HP2)[which(HP2==1, arr.ind=TRUE)[,'row']], pla=colnames(HP2)[which(HP2==1, arr.ind=TRUE)[,'col']])

cophyloplot(phyloIIB, PhyloIIB_BC_UPGMA, assoc, show.tip.label=FALSE, use.edge.length=FALSE,
            lwd=1, col='steelblue', length.line=0, gap=0, space=100,rotate=TRUE)


## PACo analysis of phylotype IIA clade

phyloIIA_subfam_matrix <- read.csv("IS_subfamily_abundance_phylogeny_IIA.csv")
rownames(phyloIIA_subfam_matrix) <- phyloIIA_subfam_matrix[, 1];
phyloIIA_subfam_matrix2 <- phyloIIA_subfam_matrix[, -1];
phyloIIA_subfam_matrix3 <- as.matrix(phyloIIA_subfam_matrix2)

phyloIIA_subfam_BC_matrix <- vegdist(phyloIIA_subfam_matrix3, method="bray")
phyloIIA_subfam_BC_matrix2 <- as.matrix(phyloIIA_subfam_BC_matrix)
phyloIIA_subfam_BC_matrix2[is.na(phyloIIA_subfam_BC_matrix2)] <- 0
PhyloIIA_BC_UPGMA <- upgma(phyloIIA_subfam_BC_matrix2)

# Start running PACo
host.D <- cophenetic(phyloIIA)
BC.D <- cophenetic(PhyloIIA_BC_UPGMA)

HP <- read.csv("IIA_PACo_binarymatrix.csv",fileEncoding="UTF-8-BOM")
row.names(HP) <- HP$X
HP <- HP[, -1];
HP2 <- as.matrix(HP)
host.D <- host.D[rownames(HP2),colnames(HP2)]
BC.D <- BC.D[colnames(HP2),colnames(HP2)] 

D <- prepare_paco_data(H=host.D,P=BC.D,HP=HP2)
D <- add_pcoord(D,correction="cailliez")
D <- PACo(D,nperm=1000,seed=12,method="r0")


D <- paco_links(D)
res <- residuals_paco(D$proc)
D$gof

# Make tanglegram
assoc <- data.frame(pol=rownames(HP2)[which(HP2==1, arr.ind=TRUE)[,'row']], pla=colnames(HP2)[which(HP2==1, arr.ind=TRUE)[,'col']])

phyloIIA_rotate <- rotateNodes(IIA,"all")

cophyloplot(phyloIIA_rotate, PhyloIIA_BC_UPGMA, assoc, show.tip.label=FALSE, use.edge.length=FALSE,
            lwd=1, col='steelblue', length.line=0, gap=0, space=70,rotate=TRUE)


## PACo analysis of phylotype III clade

phyloIII_subfam_matrix <- read.csv("IS_subfamily_abundance_phylogeny_III.csv")
rownames(phyloIII_subfam_matrix) <- phyloIII_subfam_matrix[, 1];
phyloIII_subfam_matrix2 <- phyloIII_subfam_matrix[, -1];
phyloIII_subfam_matrix3 <- as.matrix(phyloIII_subfam_matrix2)

phyloIII_subfam_BC_matrix <- vegdist(phyloIII_subfam_matrix3, method="bray")
phyloIII_subfam_BC_matrix2 <- as.matrix(phyloIII_subfam_BC_matrix)
phyloIII_subfam_BC_matrix2[is.na(phyloIII_subfam_BC_matrix2)] <- 0
PhyloIII_BC_UPGMA <- upgma(phyloIII_subfam_BC_matrix2)

# Start running PACo
host.D <- cophenetic(phyloIII)
BC.D <- cophenetic(PhyloIII_BC_UPGMA)

HP <- read.csv("III_PACo_binarymatrix.csv",fileEncoding="UTF-8-BOM")
row.names(HP) <- HP$X
HP <- HP[, -1];
HP2 <- as.matrix(HP)
host.D <- host.D[rownames(HP2),colnames(HP2)]
BC.D <- BC.D[colnames(HP2),colnames(HP2)] 

D <- prepare_paco_data(H=host.D,P=BC.D,HP=HP2)
D <- add_pcoord(D,correction="cailliez")
D <- PACo(D,nperm=1000,seed=12,method="r0")

D <- paco_links(D)
res <- residuals_paco(D$proc)
D$gof

# Make tanglegram
assoc <- data.frame(pol=rownames(HP2)[which(HP2==1, arr.ind=TRUE)[,'row']], pla=colnames(HP2)[which(HP2==1, arr.ind=TRUE)[,'col']])

cophyloplot(phyloIII, PhyloIII_BC_UPGMA, assoc, show.tip.label=FALSE, use.edge.length=FALSE,
            lwd=1, col='steelblue', length.line=0, gap=-20, space=70,rotate=TRUE)


## PACo analysis of phylotype IV clade

phyloIV_subfam_matrix <- read.csv("IS_subfamily_abundance_phylogeny_IV.csv")
rownames(phyloIV_subfam_matrix) <- phyloIV_subfam_matrix[, 1];
phyloIV_subfam_matrix2 <- phyloIV_subfam_matrix[, -1];
phyloIV_subfam_matrix3 <- as.matrix(phyloIV_subfam_matrix2)

phyloIIV_subfam_BC_matrix <- vegdist(phyloIV_subfam_matrix3, method="bray")
phyloIIV_subfam_BC_matrix2 <- as.matrix(phyloIIV_subfam_BC_matrix)
phyloIIV_subfam_BC_matrix2[is.na(phyloIIV_subfam_BC_matrix2)] <- 0
PhyloIV_BC_UPGMA <- upgma(phyloIIV_subfam_BC_matrix2)

# Start running PACo
host.D <- cophenetic(phyloIV)
BC.D <- cophenetic(PhyloIV_BC_UPGMA)

HP <- read.csv("IV_PACo_binarymatrix.csv",fileEncoding="UTF-8-BOM")
row.names(HP) <- HP$X
HP <- HP[, -1];
HP2 <- as.matrix(HP)
host.D <- host.D[rownames(HP2),colnames(HP2)]
BC.D <- BC.D[colnames(HP2),colnames(HP2)] 

D <- prepare_paco_data(H=host.D,P=BC.D,HP=HP2)
D <- add_pcoord(D,correction="cailliez")
D <- PACo(D,nperm=100,seed=12,method="r0")

D <- paco_links(D)
res <- residuals_paco(D$proc)
D$gof

# Make tanglegram
assoc <- data.frame(pol=rownames(HP2)[which(HP2==1, arr.ind=TRUE)[,'row']], pla=colnames(HP2)[which(HP2==1, arr.ind=TRUE)[,'col']])

cophyloplot(phyloIV, PhyloIV_BC_UPGMA, assoc, show.tip.label=FALSE, use.edge.length=FALSE,
            lwd=1, col='steelblue', length.line=0, gap=-20, space=70,rotate=TRUE)


#PACo analysis of all isolates together

RSSC_subfam_matrix <- read.csv("IS_subfamily_abundance_phylogeny_total.csv")
rownames(RSSC_subfam_matrix) <- RSSC_subfam_matrix[, 1];
RSSC_subfam_matrix2 <- RSSC_subfam_matrix[, -1];
RSSC_subfam_matrix3 <- as.matrix(RSSC_subfam_matrix2)

RSSC_subfam_BC_matrix <- vegdist(RSSC_subfam_matrix3, method="bray")
RSSC_subfam_BC_matrix2 <- as.matrix(RSSC_subfam_BC_matrix)
RSSC_subfam_BC_matrix2[is.na(RSSC_subfam_BC_matrix2)] <- 0
RSSC_BC_UPGMA <- upgma(RSSC_subfam_BC_matrix2)

# Start running PACo
host.D <- cophenetic(RSSC_tree)
BC.D <- cophenetic(RSSC_BC_UPGMA)

HP <- read.csv("IV_PACo_binarymatrix.csv",fileEncoding="UTF-8-BOM")
row.names(HP) <- HP$X
HP <- HP[, -1];
HP2 <- as.matrix(HP)
host.D <- host.D[rownames(HP2),colnames(HP2)]
BC.D <- BC.D[colnames(HP2),colnames(HP2)] 

# Make tanglegram
assoc <- data.frame(pol=rownames(HP2)[which(HP2==1, arr.ind=TRUE)[,'row']], pla=colnames(HP2)[which(HP2==1, arr.ind=TRUE)[,'col']])

RSSC_tree_rotated <- rotateNodes(Isolate_tree,"all")

cophyloplot(RSSC_tree_rotated, AllBC_UPGMA, assoc, show.tip.label=FALSE, use.edge.length=FALSE,
            lwd=1, col='steelblue', length.line=0, gap=-20, space=70,rotate=TRUE)


## Figure S5 - determine proportion of each IS subfamily in chromosome or megaplasmid by phylotype

Proportion_IS_chr <- read.csv("Proportion_IS_in_Chr_by_phylotype.csv",fileEncoding="UTF-8-BOM")

Proportion_IS_chr2 <- Proportion_IS_chr[, -1];
row.names(Proportion_IS_chr2) <- Proportion_IS_chr$Phylotype
Proportion_IS_chr3 <- as.matrix(Proportion_IS_chr2)
Proportion_IS_chr3_melted <- melt(Proportion_IS_chr3)

ggplot(data = Proportion_IS_chr3_melted, aes(x=Var2, y=Var1, fill=value)) + 
  geom_tile(colour="black")+
  scale_y_discrete(limits = rev)+
  theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.y=element_blank(),axis.text.x=element_blank(),legend.position="right")+
  scale_fill_scico(palette = 'cork',na.value="grey") +
  theme(plot.margin=margin(0,0,0,20)) +
  theme(legend.title = element_text(colour="black", size=12,face="bold"),legend.text = element_text(size=12))+
  labs(fill="Proportion of IS in chromosome")


# Figure S6 - Checking the difference in IS prox between chromosome and megaplasmid

ggplot(data=IStotal_long,aes(x=IS_prox_region,y=IS_prox_freq,fill=IS_prox_region)) + geom_violin(alpha=0.5, position = position_dodge(width = .75),size=1,color=NA)+
  geom_boxplot(notch = FALSE,  outlier.size = -1, color="black",lwd=1, alpha = 0.3,show.legend = F, varwidth=TRUE,position = position_dodge(width = .75))+
  ggbeeswarm::geom_quasirandom(shape = 21,size=2, dodge.width = .75, color="black",alpha=.3,show.legend = F)+
  ylab("Number of gene-proximate ISs per genome") +
  xlab("Genomic region") +
  theme_bw()+
  scale_fill_manual(values=c("#88CCEE","#0033FF"),labels = c("Chromosome", "Megaplasmid"))+
  geom_smooth(method = "lm")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title=element_text(size=12,face="bold")) +
  labs(fill="Genomic Region") +
  scale_x_discrete(labels=c("Chromosome", "Megaplasmid"))

t.test(IStotal$MP_prox_100bp , IStotal$Chr_prox_100bp, paired = TRUE, alternative = "two.sided")

# Differences are normally distributed
d <- with(IStotal_long, 
          IS_prox_freq[IS_prox_region == "Chr_prox_100bp"] - IS_prox_freq[IS_prox_region =="MP_prox_100bp"])
hist(d)


# Figure S7 - Proportion of insertions that cause disruptions caused by each IS 

ISsubtype_disrupt <- read.csv("Proportion_ISsubtype_disruptions_by_phylotype.csv",fileEncoding="UTF-8-BOM")

ISsubtype_disrupt2 <- ISsubtype_disrupt[, -1];
row.names(ISsubtype_disrupt2) <- ISsubtype_disrupt$Phylotype
ISsubtype_disrupt_mat <- as.matrix(ISsubtype_disrupt2)
ISsubtype_disrupt_melted_mat <- melt(ISsubtype_disrupt_mat)

ggplot(data = ISsubtype_disrupt_melted_mat, aes(x=Var2, y=Var1, fill=value)) + 
  geom_tile(colour="black")+
  scale_y_discrete(limits = rev)+
  scale_fill_gradient(low="grey98",high="firebrick", na.value="grey",limits=c(0,1)) + 
  theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.y=element_blank(),axis.text.x=element_blank(),legend.position="right")+
  theme(plot.margin=margin(0,0,0,20)) +
  theme(legend.title = element_text(colour="black", size=12,face="bold"),legend.text = element_text(size=12))+
  labs(fill="Proportion of IS disruptions")


## Figure 3 - Gene disruption analysis

# Figure 3 - inset

IS_distance_histogram <- read.csv("IS_distancetogene_histogram.csv",fileEncoding="UTF-8-BOM")

ggplot(IS_distance_histogram, aes(x=Shortest_distance_to_gene,fill = Region,color=Region)) + 
  geom_histogram(binwidth=10,alpha=0.3,position="identity")+
  theme(legend.position="top")+
  scale_fill_manual(values=c("#88CCEE","#0033FF"),labels = c("Chromosome", "Megaplasmid"))+
  scale_colour_manual(values=c("#88CCEE","#0033FF"))+
  theme_bw()+
  ylab("IS copy number") +
  xlab("Shortest distance to adjacent gene (bp)")+
  theme(legend.title = element_text(colour="black", size=12,face="bold"),legend.text = element_text(size=12))+
  theme(axis.title=element_text(size=12,face="bold"))+
  geom_vline(aes(xintercept=100),
             color="black", linetype="dashed", size=0.5)


# Figure 3 - bottom

IS_disrupt_morethan5_matrix <- read.csv("IS_disruption_morethan5_phylogeny.csv",fileEncoding="UTF-8-BOM")

IS_disrupt_morethan5_matrix2 <- IS_disrupt_morethan5_matrix[, -1];
row.names(IS_disrupt_morethan5_matrix2) <- IS_disrupt_morethan5_matrix$Isolate
IS_disrupt_morethan5_matrix3 <- as.matrix(IS_disrupt_morethan5_matrix2)
IS_disrupt_morethan5_matrix3_melted <- melt(IS_disrupt_morethan5_matrix3)

Fig_3_bottom <- ggplot(data = IS_disrupt_morethan5_matrix3_melted, aes(x=Var2, y=Var1, fill=value)) + 
  geom_tile()+
  scale_y_discrete(limits = rev)+
  theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.y=element_blank(),axis.text.x=element_blank(),legend.position="right")+
  scale_fill_gradient(low="grey98",high="black", limits=c(0,10)) + 
  theme(plot.margin=margin(0,0,0,20)) +
  theme(legend.title = element_text(colour="black", size=12,face="bold"),legend.text = element_text(size=12))+
  labs(fill="IS disruptions")

grid.arrange(
  arrangeGrob(
    RSSC_tree,Fig_3_bottom,
    ncol=2 ,widths=c(0.4,1)))

# Figure 3 - top

data(alphabet)
names(alphabet) <- NULL

All_subtype_disrupt <- read.csv("ISsubtype_disruption_morethan5_phylogeny.csv",fileEncoding="UTF-8-BOM")
All_subtype_disrupt_long <- gather(All_subtype_disrupt, IS_subfamily, IS_subfamily_freq, ISBma3:ISRso15, factor_key=TRUE)

All_subtype_disrupt_long$IS <- factor(All_subtype_disrupt_long$IS, levels = unique(All_subtype_disrupt_long$IS))

ggplot(All_subtype_disrupt_long, aes(fill=IS_subfamily, y=IS_subfamily_freq, x=IS)) + 
  geom_bar(position="stack", stat="identity")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme_classic()+
  scale_fill_manual(values = alphabet)+
  theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),legend.position="none")