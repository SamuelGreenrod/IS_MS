#### IS MS figure and stat R script ####

# The paper this R code was written for can be accessed at: 
# Greenrod STE, Stoycheva M, Elphinstone J, Friman V-P. 2022. Influence of insertion sequences on population structure of phytopathogenic bacteria in the Ralstonia solanacearum species complex. bioRxiv: 2022.07.16.500299.

### Dependencies ###

require("ggplot2")
require("ggbeeswarm")
require("pheatmap")
require("kableExtra")
require("tidyr")
require("lme4")
require("ggfortify")
require("MASS")
require("reshape2")
require("gridExtra")
require("corrplot")
require("treeio")
require("ggtree")
require("phytools")
require("scico")
require("vegan")
require("phangorn")
require("paco")
require("factoextra")
require("FSA")
require("dplyr")
require("cowplot")
require("rcartocolor")
require("car")
require("viridisLite")

#### RSSC phylogeny

## RSSC isolate maximum-likelihood tree

safe_colorblind_palette <- c("#CC6677", "#DDCC77", "#117733","#6699CC", "#888888")

RSSC_tree_backbone <- read.tree("core_variant_tree_356subset.treefile")
RSSC_tree <- ggtree(RSSC_tree_backbone,ladderize=TRUE)+ theme(plot.margin=margin(0,0,30,0))


## Figure S1

# Load RSSC isolate maximum-likelihood tree with representative strains

RSSC_tree_wrep <- read.tree("Samuels_phylo_march_2022.treefile")
RSSC_tree_wrep_rooted <- phytools::midpoint.root(RSSC_tree_wrep)

ggtree(rooted_tree,layout = "rectangular")+
  theme_tree2() +
  geom_tiplab(aes(subset=(grepl('GMI1000',label,fixed=TRUE)==TRUE)), size = 3)+
  geom_tiplab(aes(subset=(grepl('K60',label,fixed=TRUE)==TRUE)), size = 3)+
  geom_tiplab(aes(subset=(grepl('UY031',label,fixed=TRUE)==TRUE)), size = 3)+
  geom_tiplab(aes(subset=(grepl('CMR15',label,fixed=TRUE)==TRUE)), size = 3)+
  geom_tiplab(aes(subset=(grepl('PSI07',label,fixed=TRUE)==TRUE)), size = 3)
  


#### Identify IS clusters

# Group Nanopore genome IS into clusters
matrix <- read.csv("All_Nanopore_IS_filtered_mashmatrix_heatmap.csv")
rownames(matrix) <- matrix[, 1];
matrix2 <- matrix[, -1];
mat <- as.matrix(matrix2)

#Determine optimum number of IS clusters

fviz_nbclust(mat, FUN = kmeans,k.max = 80, method = "silhouette") # 66
pheatmap(mat, kmeans_k = 66,show_rownames = FALSE, show_colnames = FALSE)

# Sort IS into clusters
k2 <- kmeans(mat, centers = 66, nstart = 25)
IS_clusters <- data.frame(k2$cluster)
write.csv(cluster,"IS_clusters.csv")


#### Compare ISEScan detected IS (whole genome) to ISMapper detected IS (read mapping)

# Figure S2A

isescan_v_ismapper <- read.csv("isescan_v_ismapper.csv",fileEncoding="UTF-8-BOM")

Figure_S2A <- ggplot(isescan_v_ismapper,aes(x=Num_IS_ISEScan,y=Num_IS_ISMapper)) + geom_point() +
  geom_abline(intercept = 0, slope = 1,linetype="dotted")+
  theme_bw()+
  theme(axis.title=element_text(size=10,face="bold"))+
  ylab("IS copy number (ISMapper)")+
  xlab("IS copy number (ISEScan)")

# Determine correlation between IS detection methods

hist(IStotal$Num_IS_ISMapper)
hist(IStotal$Num_IS_ISEScan)
cor.test(IStotal$Num_IS_ISEScan, IStotal$Num_IS_ISMapper, method = "kendall",use = "complete.obs")


## Figure S2B (read depth)

Read_num <- read.csv("Isolate_read_depth.csv",fileEncoding="UTF-8-BOM")
Read_num$Isolate <- factor(Read_num$Isolate, levels = Read_num$Isolate)

Figure_S2B <- ggplot(data=Read_num, aes(x = Phylotype, y= Average_depth_perbase, fill = Region)) + geom_violin(alpha=0.5, position = position_dodge(width = .75),size=1,color=NA)+
  geom_boxplot(notch = FALSE,  outlier.size = -1, color="black",lwd=0.6, alpha = 0.7,show.legend = F, varwidth=TRUE,position = position_dodge(width = .75))+
  ggbeeswarm::geom_quasirandom(shape = 21,size=2, dodge.width = .75, color="black",alpha=.3,show.legend = F)+
  ylab("Average read depth per bp") +
  xlab("Phylotype") +
  theme_bw()+
  geom_smooth(method = "lm")+
  scale_fill_manual(values=c("#88CCEE","#0033FF"),labels = c("Chromosome", "Megaplasmid"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title=element_text(size=10,face="bold")) +
  theme(legend.title = element_text(colour="black", size=10,face="bold"),legend.text = element_text(size=9)) +
  labs(fill="Genomic element")

kruskal.test(Average_depth_perbase ~ Phylotype, data = Read_num)
dunnTest(Average_depth_perbase ~ Phylotype, data = Read_num)


plot_grid(Figure_S2A, Figure_S2B, align = "h", ncol=2, rel_widths = c(0.5,1))




## Figure S3 (Difference in IS number between chromosome and megaplasmid)

IStotal <- read.csv("IS_data.csv",fileEncoding="UTF-8-BOM")
IStotal_long <- gather(IStotal, IS_region, IS_freq, Chr_IS:MP_IS, factor_key=TRUE)
IStotal_long <- gather(IStotal_long, IS_prox_region, IS_prox_freq, Chr_prox_100bp:MP_prox_100bp, factor_key=TRUE)
IStotal_long <- gather(IStotal_long, IS_disrupt_region, IS_disrupt_freq, Chr_disruptions:MP_disruptions, factor_key=TRUE)
IStotal_long <- gather(IStotal_long, Prop_IS_disrupt_region, Prop_IS_disrupt_freq, Proportion_Chr_disruptions:Proportion_MP_disruptions, factor_key=TRUE)

ggplot(data=IStotal_long,aes(x=IS_region,y=IS_freq,fill=IS_region)) + geom_violin(alpha=0.5, position = position_dodge(width = .75),size=1,color=NA)+
  geom_boxplot(notch = FALSE,  outlier.size = -1, color="black",lwd=1, alpha = 0.3,show.legend = F, varwidth=TRUE,position = position_dodge(width = .75))+
  ggbeeswarm::geom_quasirandom(shape = 21,size=2, dodge.width = .75, color="black",alpha=.3,show.legend = F)+
  ylab("IS number per genome") +
  xlab("Genomic element") +
  theme_bw()+
  geom_smooth(method = "lm")+
  scale_fill_manual(values=c("#88CCEE","#0033FF"),labels = c("Chromosome", "Megaplasmid"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title=element_text(size=12,face="bold")) +
  labs(fill="Genomic element") +
  scale_x_discrete(labels=c("Chromosome", "Megaplasmid"))

wilcox.test(IStotal$MP_IS , IStotal$Chr_IS, paired = TRUE)


## Figure 1

# Figure 1A (IS family abundance in chromosome and megaplasmid)

ISfam <- read.csv("IS_family_abundance.csv",fileEncoding="UTF-8-BOM")
ISfam_long <- gather(ISfam, IS_family, IS_family_freq, IS110:ISL3, factor_key=TRUE)

a <- ggplot(data=ISfam_long, aes(x = IS_family, y= IS_family_freq, fill = Region)) +
  geom_boxplot(notch = FALSE,  outlier.size = -1, color="black",lwd=0.6, alpha = 0.9,show.legend = F, varwidth=TRUE,position = position_dodge(width = .75))+
  ggbeeswarm::geom_quasirandom(shape = 21,size=2, dodge.width = .75, color="black",alpha=.03,show.legend = F)+
  ylab("IS number per genome") +
  xlab("IS family") +
  theme_bw()+
  geom_smooth(method = "lm")+
  scale_fill_manual(values=c("#88CCEE","#0033FF"),labels = c("Chromosome", "Megaplasmid"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title=element_text(size=10,face="bold")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(legend.title = element_text(colour="black", size=12,face="bold"),legend.text = element_text(size=11)) +
  labs(fill="Genomic region")


# Figure 1B (IS abundance by phylotype)

b <- ggplot(data=IStotal_long, aes(x = Phylotype, y= IS_freq, fill = IS_region)) +
  geom_boxplot(notch = FALSE,  outlier.size = -1, color="black",lwd=0.6, alpha = 0.9,show.legend = F, varwidth=TRUE,position = position_dodge(width = .75))+
  ggbeeswarm::geom_quasirandom(shape = 21,size=2, dodge.width = .75, color="black",alpha=.03,show.legend = F)+
  ylab("IS number per genome") +
  xlab("Phylotype") +
  theme_bw()+
  geom_smooth(method = "lm")+
  scale_fill_manual(values=c("#88CCEE","#0033FF"),labels = c("Chromosome", "Megaplasmid"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title=element_text(size=10,face="bold")) +
  theme(legend.title = element_text(colour="black", size=12,face="bold"),legend.text = element_text(size=11)) +
  labs(fill="Genomic region")

kruskal.test(Num_IS ~ Phylotype, data = IStotal)
dunnTest(Num_IS ~ Phylotype, data = IStotal)

a /b + plot_layout(guides="collect")

# Figure 1C (IS family abundance across phylogeny)

ISfam <- read.csv("IS_family_abundance_phylogeny.csv",fileEncoding="UTF-8-BOM")
ISfam_long <- gather(ISfam, IS_family, IS_family_freq, IS110:ISL3, factor_key=TRUE)
ISfam_long$Isolate <- as.character(ISfam_long$Isolate)
ISfam_long$Isolate <- factor(ISfam_long$Isolate, levels=unique(ISfam_long$Isolate))

IS_fam_cols <- c("IS110" = "#5D69B1",  "IS1182" ="#52BCA3", "IS1595" = "#CC61B0", "IS21" =  "#2F8AC4", "IS256" = "#764E9F", "IS3" =  "#99C945",  "IS4"  = "#ED645A", "IS5"  = "#E58606", "IS630" = "#24796C", "IS701" = "#DAA51B", "ISL3" = "#A5AA99")


c <- ggplot(ISfam_long,aes(x=Isolate,y=IS_family_freq, fill=IS_family)) +
  geom_bar(stat="identity", position = "stack",width=1.1)+
  theme_bw()+
  theme(axis.ticks.y = element_blank()) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title=element_text(size=10,face="bold"),axis.text.y=element_blank(),axis.title.y=element_blank()) + 
  coord_flip()+ 
  scale_x_discrete(limits = rev(levels(ISfam_long$Isolate)))+
  scale_fill_manual(values=IS_fam_cols)+
  ylab("IS copy number") +
  theme(legend.title = element_text(colour="black", size=12,face="bold"),legend.text = element_text(size=12)) +
  theme(plot.margin=margin(0,0,0,20))+
  labs(fill="IS family")

plot_grid(tree, c, align = "h", ncol=2, rel_widths = c(0.25,1))



## Figure 2

# Figure 2A (IS subfamily abundance in chromosome and megaplasmid)

ISsubfam <- read.csv("IS_subfamily_abundance_phylogeny.csv",fileEncoding="UTF-8-BOM")
tree <- ggtree(RSSC_tree_backbone,ladderize=TRUE)+ theme(plot.margin=margin(0,0,0,0))

# Chromosome

Chr <- ISsubfam[(ISsubfam$Region == "Chr"),]
Chr <- subset(Chr, select = -c(Region))

Chr2 <- Chr[, -1];
row.names(Chr2) <- Chr$Isolate

Chr_mat <- as.matrix(Chr2)
Chr_mat_melt <- melt(Chr_mat)

Chr_matrix <-  ggplot(data = Chr_mat_melt, aes(x=Var2, y=Var1, fill=ifelse(value != 0, value, NA))) + 
  geom_tile()+
  scale_fill_viridis(option="inferno",limits=c(1,50), na.value="white",direction=-1) + 
  scale_y_discrete(limits = rev)+
  theme(axis.title.y=element_blank(),axis.text.y=element_blank(),legend.position="none")+
  theme(plot.margin=margin(0,10,5,20)) +
  theme(legend.title = element_text(colour="black", size=12,face="bold"),legend.text = element_text(size=12))+
  labs(fill="IS copy number")+
  xlab("IS subgroups (Chromosome)")+
  theme(axis.title=element_text(size=11,face="bold")) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size = 0.2))+
  theme(axis.ticks.y = element_blank())+
  theme(axis.text.x = element_text(angle = 55, hjust = 1))+
  geom_vline(xintercept = 1.5, color = "grey60") +
  geom_vline(xintercept = 2.5, color = "grey60") +
  geom_vline(xintercept = 3.5, color = "grey60") +
  geom_vline(xintercept = 4.5, color = "grey60") +
  geom_vline(xintercept = 5.5, color = "grey60") +
  geom_vline(xintercept = 6.5, color = "grey60") +
  geom_vline(xintercept = 7.5, color = "grey60") +
  geom_vline(xintercept = 8.5, color = "grey60") +
  geom_vline(xintercept = 9.5, color = "grey60") +
  geom_vline(xintercept = 10.5, color = "grey60") +
  geom_vline(xintercept = 11.5, color = "grey60") +
  geom_vline(xintercept = 12.5, color = "grey60") +
  geom_vline(xintercept = 13.5, color = "grey60") +
  geom_vline(xintercept = 14.5, color = "grey60") +
  geom_vline(xintercept = 15.5, color = "grey60") +
  geom_vline(xintercept = 16.5, color = "grey60") +
  geom_vline(xintercept = 17.5, color = "grey60") +
  geom_vline(xintercept = 18.5, color = "grey60") +
  geom_vline(xintercept = 19.5, color = "grey60") +
  geom_vline(xintercept = 20.5, color = "grey60") +
  geom_vline(xintercept = 21.5, color = "grey60") +
  geom_vline(xintercept = 22.5, color = "grey60") +
  geom_vline(xintercept = 23.5, color = "grey60") +
  geom_vline(xintercept = 24.5, color = "grey60") +
  geom_vline(xintercept = 25.5, color = "grey60")

# Megaplasmid

MP <- ISsubfam[(ISsubfam$Region == "MP"),]
MP <- subset(MP, select = -c(Region))

MP2 <- MP[, -1];
row.names(MP2) <- MP$Isolate
MP_mat <- as.matrix(MP2)
MP_mat_melt <- melt(MP_mat)
MP_mat_melt <- na.omit(MP_mat_melt)

MP_matrix <- ggplot(data = MP_mat_melt, aes(x=Var2, y=Var1, fill=ifelse(value != 0, value, NA))) + 
  geom_tile()+
  scale_fill_viridis(option="inferno",limits=c(1,50), na.value="white",direction=-1) + 
  scale_y_discrete(limits = rev)+
  theme(axis.title.y=element_blank(),axis.text.y=element_blank(),legend.position="right")+
  theme(plot.margin=margin(0,0,0,10)) +
  theme(legend.title = element_text(colour="black", size=12,face="bold"),legend.text = element_text(size=12))+
  labs(fill="IS copy number")+
  theme(axis.title=element_text(size=11,face="bold")) +
  xlab("IS subgroups (Megaplasmid)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border = element_rect(colour = "black", fill=NA, size = 0.2))+
  theme(axis.ticks.y = element_blank())+
  theme(axis.text.x = element_text(angle = 55, hjust = 1))+
  geom_vline(xintercept = 1.5, color = "grey60") +
  geom_vline(xintercept = 2.5, color = "grey60") +
  geom_vline(xintercept = 3.5, color = "grey60") +
  geom_vline(xintercept = 4.5, color = "grey60") +
  geom_vline(xintercept = 5.5, color = "grey60") +
  geom_vline(xintercept = 6.5, color = "grey60") +
  geom_vline(xintercept = 7.5, color = "grey60") +
  geom_vline(xintercept = 8.5, color = "grey60") +
  geom_vline(xintercept = 9.5, color = "grey60") +
  geom_vline(xintercept = 10.5, color = "grey60") +
  geom_vline(xintercept = 11.5, color = "grey60") +
  geom_vline(xintercept = 12.5, color = "grey60") +
  geom_vline(xintercept = 13.5, color = "grey60") +
  geom_vline(xintercept = 14.5, color = "grey60") +
  geom_vline(xintercept = 15.5, color = "grey60") +
  geom_vline(xintercept = 16.5, color = "grey60") +
  geom_vline(xintercept = 17.5, color = "grey60") +
  geom_vline(xintercept = 18.5, color = "grey60") +
  geom_vline(xintercept = 19.5, color = "grey60") +
  geom_vline(xintercept = 20.5, color = "grey60") +
  geom_vline(xintercept = 21.5, color = "grey60") +
  geom_vline(xintercept = 22.5, color = "grey60") +
  geom_vline(xintercept = 23.5, color = "grey60") +
  geom_vline(xintercept = 24.5, color = "grey60") +
  geom_vline(xintercept = 25.5, color = "grey60")

tree + Chr_matrix + MP_matrix + plot_layout(guides = "collect",widths = c(1, 3,3))



## Figure 2B (PCoA analysis of phylotype IS content)

matrix <- read.csv("IS_subfamily_abundance_phylogeny_total.csv",fileEncoding="UTF-8-BOM")
matrix2 <- matrix
matrix2$Phylotype <- NULL
rownames(matrix2) <- matrix2[, 1];
matrix3 <- matrix2[, -1];
mat <- as.matrix(matrix3)

x <- vegdist(mat, method="bray")

PcoA_phylotype <- pcoa(x)
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


##Figure 2C/Figure S4 (PACo analysis looking at congruence between IS content and host genetic distance)

# Splitting tree

I <- tree_subset(RSSC_tree, node=643, levels_back=0)
IIA <- tree_subset(RSSC_tree, node=627, levels_back=0)
IIB <- tree_subset(RSSC_tree, node=359, levels_back=0)
III <- tree_subset(RSSC_tree, node=701, levels_back=0)
IV <- tree_subset(RSSC_tree, node=708, levels_back=0)


## PACo analysis of phylotype I clade

matrix <- read.csv("IS_subfamily_abundance_phylogeny_I.csv")
rownames(matrix) <- matrix[, 1];
matrix2 <- matrix[, -1];
mat <- as.matrix(matrix2)

x <- vegdist(mat, method="bray")
PhyloI_BC <- as.matrix(x)
PhyloI_BC[is.na(PhyloI_BC)] <- 0
PhyloI_UPGMA <- upgma(PhyloI_BC)

# Start running PACo
host.D <- cophenetic(I)
BC.D <- cophenetic(PhyloI_UPGMA)

HP <- read.csv("I_PACo_binarymatrix.csv")
row.names(HP) <- HP$X
HP <- HP[, -1];
HP2 <- as.matrix(HP)
host.D <- host.D[rownames(HP2),colnames(HP2)]
BC.D <- BC.D[colnames(HP2),colnames(HP2)] 

D <- prepare_paco_data(H=host.D,P=BC.D,HP=HP2)
D <- add_pcoord(D,correction="cailliez")
D <- PACo(D,nperm=10,seed=12,method="r0")

D <- paco_links(D)
res <- residuals_paco(D$proc)

# Congruence result
D$gof

assoc <- data.frame(pol=rownames(HP2)[which(HP2==1, arr.ind=TRUE)[,'row']], pla=colnames(HP2)[which(HP2==1, arr.ind=TRUE)[,'col']])
cophyloplot(I, PhyloI_UPGMA, assoc, show.tip.label=FALSE, use.edge.length=FALSE,
            lwd=1, col='steelblue', length.line=0, gap=-15, space=70,rotate=TRUE)



## PACo analysis of phylotype IIB clade

matrix <- read.csv("IS_subfamily_abundance_phylogeny_IIB.csv")
rownames(matrix) <- matrix[, 1];
matrix2 <- matrix[, -1];
mat <- as.matrix(matrix2)

x <- vegdist(mat, method="bray")
PhyloIIB_BC <- as.matrix(x)
PhyloIIB_BC[is.na(PhyloIIB_BC)] <- 0
PhyloIIB_UPGMA <- upgma(PhyloIIB_BC)

# Start running PACo
host.D <- cophenetic(IIB)
BC.D <- cophenetic(PhyloIIB_UPGMA)

HP <- read.csv("IIB_PACo_binarymatrix.csv")
row.names(HP) <- HP$X
HP <- HP[, -1];
HP2 <- as.matrix(HP)
host.D <- host.D[rownames(HP2),colnames(HP2)]
BC.D <- BC.D[colnames(HP2),colnames(HP2)] 

D <- prepare_paco_data(H=host.D,P=BC.D,HP=HP2)
D <- add_pcoord(D,correction="cailliez")
D <- PACo(D,nperm=10,seed=12,method="r0")

D <- paco_links(D)
res <- residuals_paco(D$proc)

# Congruence result
D$gof

assoc <- data.frame(pol=rownames(HP2)[which(HP2==1, arr.ind=TRUE)[,'row']], pla=colnames(HP2)[which(HP2==1, arr.ind=TRUE)[,'col']])
cophyloplot(IIB, PhyloIIB_UPGMA, assoc, show.tip.label=FALSE, use.edge.length=FALSE,
            lwd=1, col='steelblue', length.line=0, gap=0, space=100,rotate=TRUE)



## PACo analysis of phylotype IIA clade

matrix <- read.csv("IS_subfamily_abundance_phylogeny_IIA.csv")
rownames(matrix) <- matrix[, 1];
matrix2 <- matrix[, -1];
mat <- as.matrix(matrix2)

x <- vegdist(mat, method="bray")
PhyloIIA_BC <- as.matrix(x)
PhyloIIA_BC[is.na(PhyloIIA_BC)] <- 0
PhyloIIA_UPGMA <- upgma(PhyloIIA_BC)

# Start running PACo
host.D <- cophenetic(IIA)
BC.D <- cophenetic(PhyloIIA_UPGMA)

HP <- read.csv("IIA_PACo_binarymatrix.csv")
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

# Congruence result
D$gof

assoc <- data.frame(pol=rownames(HP2)[which(HP2==1, arr.ind=TRUE)[,'row']], pla=colnames(HP2)[which(HP2==1, arr.ind=TRUE)[,'col']])
IIA_2 <- rotateNodes(IIA,"all")
cophyloplot(IIA_2, PhyloIIA_UPGMA, assoc, show.tip.label=FALSE, use.edge.length=FALSE,
            lwd=1, col='steelblue', length.line=0, gap=0, space=70,rotate=TRUE)



## PACo analysis of phylotype III clade

matrix <- read.csv("IS_subfamily_abundance_phylogeny_III.csv")
rownames(matrix) <- matrix[, 1];
matrix2 <- matrix[, -1];
mat <- as.matrix(matrix2)

x <- vegdist(mat, method="bray")
PhyloIII_BC <- as.matrix(x)
PhyloIII_BC[is.na(PhyloIII_BC)] <- 0
PhyloIII_UPGMA <- upgma(PhyloIII_BC)

# Start running PACo
host.D <- cophenetic(III)
BC.D <- cophenetic(PhyloIII_UPGMA)

HP <- read.csv("III_PACo_binarymatrix.csv")
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

# Congruence result
D$gof

assoc <- data.frame(pol=rownames(HP2)[which(HP2==1, arr.ind=TRUE)[,'row']], pla=colnames(HP2)[which(HP2==1, arr.ind=TRUE)[,'col']])
cophyloplot(III, PhyloIII_UPGMA, assoc, show.tip.label=FALSE, use.edge.length=FALSE,
            lwd=1, col='steelblue', length.line=0, gap=-20, space=70,rotate=TRUE)



## PACo analysis of phylotype IV clade

matrix <- read.csv("IS_subfamily_abundance_phylogeny_IV.csv")
rownames(matrix) <- matrix[, 1];
matrix2 <- matrix[, -1];
mat <- as.matrix(matrix2)

x <- vegdist(mat, method="bray")
PhyloIV_BC <- as.matrix(x)
PhyloIV_BC[is.na(PhyloIV_BC)] <- 0
PhyloIV_UPGMA <- upgma(PhyloIV_BC)

# Start running PACo
host.D <- cophenetic(IV)
BC.D <- cophenetic(PhyloIV_UPGMA)

HP <- read.csv("IV_PACo_binarymatrix.csv")
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

# Congruence result
D$gof

assoc <- data.frame(pol=rownames(HP2)[which(HP2==1, arr.ind=TRUE)[,'row']], pla=colnames(HP2)[which(HP2==1, arr.ind=TRUE)[,'col']])
cophyloplot(IV, PhyloIV_UPGMA, assoc, show.tip.label=FALSE, use.edge.length=FALSE,
            lwd=1, col='steelblue', length.line=0, gap=-20, space=70,rotate=TRUE)




#PACo analysis of all isolates together

matrix <- read.csv("IS_subfamily_abundance_phylogeny_total.csv",fileEncoding="UTF-8-BOM")
matrix2 <- matrix
matrix2$Phylotype <- NULL
rownames(matrix2) <- matrix2[, 1];
matrix3 <- matrix2[, -1];
mat <- as.matrix(matrix3)

x <- vegdist(mat, method="bray")
All_BC <- as.matrix(x)
All_BC[is.na(All_BC)] <- 0
AllBC_UPGMA <- upgma(All_BC)

# Start running PACo
host.D <- cophenetic(RSSC_tree_backbone)
BC.D <- cophenetic(AllBC_UPGMA)

HP <- read.csv("All_isolate_PACo_binarymatrix.csv")
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

# Congruence result
D$gof

assoc <- data.frame(pol=rownames(HP2)[which(HP2==1, arr.ind=TRUE)[,'row']], pla=colnames(HP2)[which(HP2==1, arr.ind=TRUE)[,'col']])
Isolate_tree2 <- rotateNodes(Isolate_tree,"all")
cophyloplot(Isolate_tree2, AllBC_UPGMA, assoc, show.tip.label=FALSE, use.edge.length=FALSE,
            lwd=1, col='steelblue', length.line=0, gap=-20, space=70,rotate=TRUE)



## Figure S5

Proportion_chr <- read.csv("Proportion_IS_in_Chr_by_phylotype.csv",fileEncoding="UTF-8-BOM")

Proportion_chr2 <- Proportion_chr[, -1];
row.names(Proportion_chr2) <- Proportion_chr$Phylotype
Proportion_chr_mat <- as.matrix(Proportion_chr2)
Proportion_chr_melted_mat <- melt(Proportion_chr_mat)


Proportion_chr_label <- read.csv("IS_subgroup_chromosome_labels.csv",fileEncoding="UTF-8-BOM")

Proportion_chr_label2 <- Proportion_chr_label[, -1];
row.names(Proportion_chr_label2) <- Proportion_chr_label$Phylotype
Proportion_chr_mat_label <- as.matrix(Proportion_chr_label2)
Proportion_chr_melted_mat_label <- melt(Proportion_chr_mat_label)

b <- c(0, 0.5, 1)

ggplot(data = Proportion_chr_melted_mat, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile(colour="black")+
  geom_text(aes(label=Proportion_chr_melted_mat_label$value),size=3) +
  scale_y_discrete(limits = rev)+
  #xlim(-0.2,1.2)+
  xlab("Phylotype")+
  theme(axis.title.y=element_blank(),axis.title.x=element_blank(),legend.position="right")+
  scale_fill_scico(palette = 'cork',na.value="grey", limits=c(-0.3,1.3),breaks=b, labels=format(b)) +
  theme(plot.margin=margin(0,0,0,20)) +
  theme(legend.title = element_text(colour="black", size=12,face="bold"),legend.text = element_text(size=12))+
  labs(fill="Proportion of IS in chromosome")



## Figure 3 (IS positions distributed across phylogeny)

IS_info <- read.csv("All_ISMapper_ISFinder_ISIncl_filtered_intraISremoved_overlapsremoved.csv", fileEncoding="UTF-8-BOM")

# Phylotype I isolates

I_tree <- ggtree(I,ladderize=TRUE)+ theme_tree2(fontsize-0.8)+theme(plot.margin=margin(0,0,30,10))

x <- ladderize(I, right=FALSE)
is_tip <- x$edge[,2] <= length(x$tip.label)

# Extract tip order to order IS position matrix
ordered_tips <- x$edge[is_tip, 2]
ordered_tips_df <- data.frame(x$tip.label[ordered_tips])

ordered_tips_df$x.tip.label.ordered_tips. <- substr(ordered_tips_df$x.tip.label.ordered_tips., 0, 5)
ordered_tipsdf_list <- ordered_tips_df$x.tip.label.ordered_tips.
ordered_tipsdf_listrev <- rev(ordered_tipsdf_list)

## Phylotype I chromosome positions

phyloI_chr_IS <- IS_info[(IS_info$Phylotype=="I" & IS_info$Region=="Chr"),]

phyloI_chr_IS_order <- phyloI_chr_IS %>%
  mutate(Isolate =  factor(Isolate, levels = ordered_tipsdf_list)) %>%
  arrange(Isolate) 

phyloI_ispos_chr_fam <- ggplot(data=phyloI_chr_IS_order, aes(x=x, y=Isolate, colour = IS_family))+
  geom_point(size=1)+
  theme_bw()+
  theme(axis.ticks.y = element_blank()) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title=element_text(size=12,face="bold"),axis.text.y=element_blank(),axis.title.y=element_blank()) + 
  scale_colour_manual(values=IS_fam_cols,
                      name="IS family",
                      drop=F)+
  xlab("Chromosome position (bp)") +
  theme(legend.title = element_text(colour="black", size=12,face="bold"),legend.text = element_text(size=12)) +
  theme(plot.margin=margin(0,0,0,10))+
  guides(color = guide_legend(
    override.aes=list(shape = 15, size=2)))

plot_grid(I_tree, phyloI_ispos_chr_fam, align = "h", ncol=2, rel_widths = c(0.25,0.75))


## Phylotype I megaplasmid positions

phyloI_mp_IS <- IS_info[(IS_info$Phylotype=="I" & IS_info$Region=="MP"),]

phyloI_mp_IS_order <- phyloI_mp_IS %>%
  mutate(Isolate =  factor(Isolate, levels = ordered_tipsdf_list)) %>%
  arrange(Isolate) 

phyloI_ispos_mp_fam <- ggplot(data=phyloI_mp_IS_order, aes(x=x, y=Isolate, colour = IS_family))+
  geom_point(size=1)+
  theme_bw()+
  theme(axis.ticks.y = element_blank()) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title=element_text(size=12,face="bold"),axis.text.y=element_blank(),axis.title.y=element_blank()) + 
  scale_colour_manual(values=IS_fam_cols,name="IS family")+
  xlab("Megaplasmid position (bp)") +
  theme(legend.title = element_text(colour="black", size=12,face="bold"),legend.text = element_text(size=12)) +
  theme(plot.margin=margin(0,0,0,10))+
  guides(color = guide_legend(
    override.aes=list(shape = 15, size=2)))

plot_grid(I_tree, phyloI_ispos_mp_fam, align = "h", ncol=2, rel_widths = c(0.25,0.75))


#### Phylotype I IS position congruence analysis

matrix <- read.csv("IS_location_presabs_phylotype_I.csv")
rownames(matrix) <- matrix[, 1];
matrix2 <- matrix[, -1];
mat <- as.matrix(matrix2)

x <- vegdist(mat, method="jaccard")
PhyloI_BC <- as.matrix(x)
PhyloI_BC[is.na(PhyloI_BC)] <- 0
PhyloI_UPGMA <- upgma(PhyloI_BC)

# Start running PACo
host.D <- cophenetic(I)
BC.D <- cophenetic(PhyloI_UPGMA)

HP <- read.csv("I_PACo_binarymatrix.csv")
row.names(HP) <- HP$X
HP <- HP[, -1];
HP2 <- as.matrix(HP)
host.D <- host.D[rownames(HP2),colnames(HP2)]
BC.D <- BC.D[colnames(HP2),colnames(HP2)] 

D <- prepare_paco_data(H=host.D,P=BC.D,HP=HP2)
D <- add_pcoord(D,correction="cailliez")
D <- PACo(D,nperm=10000,seed=12,method="r0")

D <- paco_links(D)
res <- residuals_paco(D$proc)
D$gof

assoc <- data.frame(pol=rownames(HP2)[which(HP2==1, arr.ind=TRUE)[,'row']], pla=colnames(HP2)[which(HP2==1, arr.ind=TRUE)[,'col']])
I_2 <- rotateNodes(I,"all")
cophyloplot(I_2, PhyloI_UPGMA, assoc, show.tip.label=FALSE, use.edge.length=FALSE,
            lwd=1, col='steelblue', length.line=0, gap=-15, space=70,rotate=TRUE)


# Phylotype IIB isolates

IIB_tree <- ggtree(IIB,ladderize=TRUE)+ theme_tree2()+theme(plot.margin=margin(0,0,30,10))

x <- ladderize(IIB, right=FALSE)
is_tip <- x$edge[,2] <= length(x$tip.label)

# Extract tip order to order IS position matrix
ordered_tips <- x$edge[is_tip, 2]
ordered_tips_df <- data.frame(x$tip.label[ordered_tips])
ordered_tips_df$x.tip.label.ordered_tips. <- substr(ordered_tips_df$x.tip.label.ordered_tips., 0, 5)
ordered_tipsdf_list <- ordered_tips_df$x.tip.label.ordered_tips.
ordered_tipsdf_listrev <- rev(ordered_tipsdf_list)


## Phylotype IIB chromosome positions

phyloIIB_chr_IS <- IS_info[(IS_info$Phylotype=="IIB" & IS_info$Region=="Chr"),]

phyloIIB_chr_IS_order <- phyloIIB_chr_IS %>%
  mutate(Isolate =  factor(Isolate, levels = ordered_tipsdf_list)) %>%
  arrange(Isolate) 

phyloIIB_ispos_chr_fam <- ggplot(data=phyloIIB_chr_IS_order, aes(x=x, y=Isolate, colour = IS_family))+
  geom_point(size=1)+
  theme_bw()+
  theme(axis.ticks.y = element_blank()) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title=element_text(size=12,face="bold"),axis.text.y=element_blank(),axis.title.y=element_blank()) + 
  scale_colour_manual(values=IS_fam_cols,name="IS family")+
  xlab("Chromosome position (bp)") +
  theme(legend.title = element_text(colour="black", size=12,face="bold"),legend.text = element_text(size=12)) +
  theme(plot.margin=margin(0,0,0,10))+
  guides(color = guide_legend(
    override.aes=list(shape = 15, size=2)))

plot_grid(IIB_tree, phyloIIB_ispos_chr_fam, align = "h", ncol=2, rel_widths = c(0.25,0.75))


## Phylotype IIB megaplasmid positions

phyloIIB_mp_IS <- IS_info[(IS_info$Phylotype=="IIB" & IS_info$Region=="MP"),]

phyloIIB_mp_IS_order <- phyloIIB_mp_IS %>%
  mutate(Isolate =  factor(Isolate, levels = ordered_tipsdf_list)) %>%
  arrange(Isolate) 

phyloIIB_ispos_mp_fam <- ggplot(data=phyloIIB_mp_IS_order, aes(x=x, y=Isolate, colour = IS_family))+
  geom_point(size=1)+
  theme_bw()+
  theme(axis.ticks.y = element_blank()) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title=element_text(size=12,face="bold"),axis.text.y=element_blank(),axis.title.y=element_blank()) + 
  scale_colour_manual(values=IS_fam_cols,name="IS family")+
  xlab("Megaplasmid position (bp)") +
  theme(legend.title = element_text(colour="black", size=12,face="bold"),legend.text = element_text(size=12)) +
  theme(plot.margin=margin(0,0,0,10))+
  guides(color = guide_legend(
    override.aes=list(shape = 15, size=2)))

plot_grid(IIB_tree, phyloIIB_ispos_mp_fam, align = "h", ncol=2, rel_widths = c(0.25,0.75))


## Phylotype IIB IS position congruence analysis

matrix <- read.csv("IS_location_presabs_phylotype_IIB.csv")
rownames(matrix) <- matrix[, 1];
matrix2 <- matrix[, -1];
mat <- as.matrix(matrix2)

x <- vegdist(mat, method="jaccard")
PhyloI_BC <- as.matrix(x)
PhyloI_BC[is.na(PhyloI_BC)] <- 0
PhyloI_UPGMA <- upgma(PhyloI_BC)

# Start running PACo
host.D <- cophenetic(IIB)
BC.D <- cophenetic(PhyloI_UPGMA)

HP <- read.csv("IIB_PACo_binarymatrix.csv")
row.names(HP) <- HP$X
HP <- HP[, -1];
HP2 <- as.matrix(HP)
host.D <- host.D[rownames(HP2),colnames(HP2)]
BC.D <- BC.D[colnames(HP2),colnames(HP2)] 

D <- prepare_paco_data(H=host.D,P=BC.D,HP=HP2)
D <- add_pcoord(D,correction="cailliez")
D <- PACo(D,nperm=10,seed=12,method="r0")

D <- paco_links(D)
res <- residuals_paco(D$proc)

# Congruence result
D$gof

assoc <- data.frame(pol=rownames(HP2)[which(HP2==1, arr.ind=TRUE)[,'row']], pla=colnames(HP2)[which(HP2==1, arr.ind=TRUE)[,'col']])
cophyloplot(IIB, PhyloI_UPGMA, assoc, show.tip.label=FALSE, use.edge.length=FALSE,
            lwd=1, col='steelblue', length.line=0, gap=0, space=70,rotate=TRUE)


## Figure 3C + D (Prevalence of IS positions within phylotypes (measure of IS activity))

IS_fam_Activity <- read.csv("IS_fam_activity.csv",fileEncoding="UTF-8-BOM")

phyloI_IS <- IS_fam_Activity[(IS_fam_Activity$Phylotype=="I"),]
phyloI_IS$IS_family <- as.factor(phyloI_IS$IS_family)
phyloI_IS$Number_isolates_with_IS_in_pos <- as.numeric(phyloI_IS$Number_isolates_with_IS_in_pos)
phyloI_IS$percent <- as.numeric(phyloI_IS$percent)

phyloI_ISactivity <- ggplot(phyloI_IS, aes(x = percent,fill = IS_family)) + 
  #call geom_histogram with position="dodge" to offset the bars and manual binwidth of 2
  geom_histogram(position = "stack", binwidth=0.02)+
  ylim(0,350)+
  xlab("Proportion of isolates with IS insertion position") +
  ylab("Number of IS insertion positions") +
  scale_fill_manual(values=IS_fam_cols,name="IS family")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.title=element_text(size=12,face="bold"))


phyloIIB_IS <- IS_fam_Activity[(IS_fam_Activity$Phylotype=="IIB"),]
phyloIIB_IS$IS_family <- as.factor(phyloIIB_IS$IS_family)
phyloIIB_IS$Number_isolates_with_IS_in_pos <- as.numeric(phyloIIB_IS$Number_isolates_with_IS_in_pos)
phyloIIB_IS$percent <- as.numeric(phyloIIB_IS)

phyloIIB_ISactivity <- ggplot(phyloIIB_IS, aes(x = percent,fill = IS_family)) + 
  #call geom_histogram with position="dodge" to offset the bars and manual binwidth of 2
  geom_histogram(position = "stack", binwidth=0.02)+
  ylim(0,350)+
  xlab("Proportion of isolates with IS insertion position") +
  scale_fill_manual(values=IS_fam_cols,name="IS family")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.title.y=element_blank())+
  theme(axis.title=element_text(size=12,face="bold"))

phyloI_ISactivity + phyloIIB_ISactivity  + plot_layout(guides = "collect",widths = c(1, 1))


## Figure S6 (Phylotype IIA IS positions)

IIA_tree <- ggtree(IIA,ladderize=TRUE)+ theme_tree2()+theme(plot.margin=margin(0,0,30,10))

x <- ladderize(IIA, right=FALSE)
is_tip <- x$edge[,2] <= length(x$tip.label)

# Extract tip order to order IS position matrix
ordered_tips <- x$edge[is_tip, 2]
ordered_tips_df <- data.frame(x$tip.label[ordered_tips])

ordered_tips_df$x.tip.label.ordered_tips. <- substr(ordered_tips_df$x.tip.label.ordered_tips., 0, 5)
ordered_tipsdf_list <- ordered_tips_df$x.tip.label.ordered_tips.
ordered_tipsdf_listrev <- rev(ordered_tipsdf_list)


## Phylotype IIA chromosome positions

phyloIIA_chr_IS <- IS_info[(IS_info$Phylotype=="IIA" & IS_info$Region=="Chr"),]

phyloIIA_chr_IS_order <- phyloIIA_chr_IS %>%
  mutate(Isolate =  factor(Isolate, levels = ordered_tipsdf_list)) %>%
  arrange(Isolate) 

phyloIIA_ispos_chr_fam <- ggplot(data=phyloIIA_chr_IS_order, aes(x=x, y=Isolate, colour = IS_family))+
  geom_point(size=1)+
  theme_bw()+
  theme(axis.ticks.y = element_blank()) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title=element_text(size=12,face="bold"),axis.text.y=element_blank(),axis.title.y=element_blank()) + 
  scale_colour_manual(values=IS_fam_cols,name="IS family")+
  xlab("Chromosome position (bp)") +
  theme(legend.title = element_text(colour="black", size=12,face="bold"),legend.text = element_text(size=12)) +
  theme(plot.margin=margin(0,0,0,10))+
  guides(color = guide_legend(
    override.aes=list(shape = 15, size=2)))

plot_grid(IIA_tree, phyloIIA_ispos_chr_fam, align = "h", ncol=2, rel_widths = c(0.25,0.75))

## Phylotype IIA megaplasmid positions

phyloIIA_mp_IS <- IS_info[(IS_info$Phylotype=="IIA" & IS_info$Region=="MP"),]

phyloIIA_mp_IS_order <- phyloIIA_mp_IS %>%
  mutate(Isolate =  factor(Isolate, levels = ordered_tipsdf_list)) %>%
  arrange(Isolate) 

phyloIIA_ispos_mp_fam <- ggplot(data=phyloIIA_mp_IS_order, aes(x=x, y=Isolate, colour = IS_family))+
  geom_point(size=1)+
  theme_bw()+
  theme(axis.ticks.y = element_blank()) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title=element_text(size=12,face="bold"),axis.text.y=element_blank(),axis.title.y=element_blank()) + 
  scale_colour_manual(values=IS_fam_cols,name="IS family")+
  xlab("Megaplasmid position (bp)") +
  theme(legend.title = element_text(colour="black", size=12,face="bold"),legend.text = element_text(size=12)) +
  theme(plot.margin=margin(0,0,0,10))+
  guides(color = guide_legend(
    override.aes=list(shape = 15, size=2)))

plot_grid(IIA_tree, phyloIIA_ispos_mp_fam, align = "h", ncol=2, rel_widths = c(0.25,0.75))


# Figure S7 (IS distance to neighbouring genes)

hist <- read.csv("IS_distancetogene_histogram.csv",fileEncoding="UTF-8-BOM")

ggplot(hist, aes(x=Shortest_distance_to_gene,fill = Region,color=Region)) + 
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


## Figure S8 (Difference in IS prox between chromosome and megaplasmid)

ggplot(data=IStotal_long,aes(x=IS_prox_region,y=IS_prox_freq,fill=IS_prox_region)) + geom_violin(alpha=0.5, position = position_dodge(width = .75),size=1,color=NA)+
  geom_boxplot(notch = FALSE,  outlier.size = -1, color="black",lwd=1, alpha = 0.3,show.legend = F, varwidth=TRUE,position = position_dodge(width = .75))+
  ggbeeswarm::geom_quasirandom(shape = 21,size=2, dodge.width = .75, color="black",alpha=.3,show.legend = F)+
  ylab("Number of gene-proximate ISs per genome") +
  xlab("Genomic element") +
  theme_bw()+
  scale_fill_manual(values=c("#88CCEE","#0033FF"),labels = c("Chromosome", "Megaplasmid"))+
  geom_smooth(method = "lm")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title=element_text(size=11,face="bold")) +
  labs(fill="Genomic element") +
  scale_x_discrete(labels=c("Chromosome", "Megaplasmid"))

wilcox.test(IStotal$MP_prox_100bp , IStotal$Chr_prox_100bp, paired = TRUE)



## Figure 4 (gene disruptions in RSSC isolates)

# Figure 4A (Phylotype I gene disruptions)

# Phylotype I chromosome disruptions

phyloI_chr_IS_genedis <- IS_info[(IS_info$Phylotype=="I" & IS_info$Region=="Chr" & IS_info$gene_interruption=="TRUE"),]

# Add in isolates that do not contain IS disruptions
genedis <- unique(phyloI_chr_IS_genedis$Isolate)
unique_isolates <- unique(phyloI_chr_IS$Isolate)

missing_vals <- tibble(setdiff(unique_isolates, genedis))
names(missing_vals)[1] <- "Isolate"

phyloI_chr_IS_genedis_edit <- bind_rows(phyloI_chr_IS_genedis,missing_vals)

# Extract tip order to order IS position matrix
x <- ladderize(I, right=FALSE)
is_tip <- x$edge[,2] <= length(x$tip.label)
ordered_tips <- x$edge[is_tip, 2]
ordered_tips_df <- data.frame(x$tip.label[ordered_tips])

ordered_tips_df$x.tip.label.ordered_tips. <- substr(ordered_tips_df$x.tip.label.ordered_tips., 0, 5)
ordered_tipsdf_list <- ordered_tips_df$x.tip.label.ordered_tips.
ordered_tipsdf_listrev <- rev(ordered_tipsdf_list)

phyloI_chr_IS_genedis_order <- phyloI_chr_IS_genedis_edit %>%
  mutate(Isolate = factor(Isolate, levels = ordered_tipsdf_list)) %>%
  arrange(Isolate) 

phyloI_ispos_chr_genedis_fam <- ggplot(data=phyloI_chr_IS_genedis_order, aes(x=x, y=Isolate,colour=IS_family))+
  geom_point(size=1)+
  theme_bw()+
  theme(axis.ticks.y = element_blank()) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title=element_text(size=12,face="bold"),axis.text.y=element_blank(),axis.title.y=element_blank()) + 
  scale_colour_manual(values=IS_fam_cols,name="IS family")+
  xlab("Chromosome position (bp)") +
  theme(legend.title = element_text(colour="black", size=12,face="bold"),legend.text = element_text(size=12)) +
  theme(plot.margin=margin(0,0,0,10))+
  guides(color = guide_legend(
    override.aes=list(shape = 15, size=2)))

plot_grid(I_tree, phyloI_ispos_chr_genedis_fam, align = "h", ncol=2, rel_widths = c(0.25,0.75))


## Phylotype I megaplasmid disruptions

phyloI_mp_IS <- IS_info[(IS_info$Phylotype=="I" & IS_info$Region=="MP"),]

phyloI_mp_IS_order <- phyloI_mp_IS %>%
  mutate(Isolate =  factor(Isolate, levels = ordered_tipsdf_list)) %>%
  arrange(Isolate) 

phyloI_ispos_mp_fam <- ggplot(data=phyloI_mp_IS_order, aes(x=x, y=Isolate, colour = IS_family))+
  geom_point(size=1)+
  theme_bw()+
  theme(axis.ticks.y = element_blank()) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title=element_text(size=12,face="bold"),axis.text.y=element_blank(),axis.title.y=element_blank()) + 
  scale_colour_manual(values=IS_fam_cols,name="IS family")+
  xlab("Megaplasmid position (bp)") +
  theme(legend.title = element_text(colour="black", size=12,face="bold"),legend.text = element_text(size=12)) +
  theme(plot.margin=margin(0,0,0,10))+
  guides(color = guide_legend(
    override.aes=list(shape = 15, size=2)))

plot_grid(I_tree, phyloI_ispos_mp_fam, align = "h", ncol=2, rel_widths = c(0.25,0.75))


# Figure 4C (Phylotype I IS disruptions > 5 isolates coloured by IS family)

phyloI_IS_genedis <- IS_info[(IS_info$Phylotype=="I" & IS_info$gene_interruption=="TRUE"),]

phyloI_IS_genedis_filter <- phyloI_IS_genedis %>% 
  group_by(Round) %>% 
  mutate(freq = n()) %>% 
  ungroup() %>% 
  filter(freq > 5) %>%
  select(-freq)

phyloI_IS_genedis_num <- phyloI_IS_genedis_filter %>% group_by(Round,Region,IS_family) %>% summarize(n=n()) %>% mutate(freq=n/sum(n)) 
newtab <- merge(phyloI_IS_genedis_num, phyloI_IS_genedis[, c("Round", "right_description")], by="Round")

region_names <- c('Chr' = "Chromosome", 'MP' = "Megaplasmid")

ggplot(data=newtab, aes(x=reorder(Round, -n),fill=IS_family))+
  geom_bar(aes(y=n), stat="identity", position = "dodge")+
  theme_bw()+
  ylim(0,59)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title=element_text(size=12,face="bold")) + 
  scale_fill_manual(values=IS_fam_cols,name="IS family")+
  scale_x_discrete(breaks = newtab$Round, labels = newtab$right_description)+
  xlab("Gene interrupted") +
  ylab("IS insertion frequency")+
  theme(legend.title = element_text(colour="black", size=12,face="bold"),legend.text = element_text(size=12)) +
  theme(plot.margin=margin(0,0,0,90))+
  guides(color = guide_legend(
    override.aes=list(shape = 15, size=2)))+
  theme(axis.text.x = element_text(angle = 55, hjust = 1))+
  facet_grid(~ Region, scales = "free", space = "free",labeller = as_labeller(region_names))


## Figure S9 (Phylotype I IS disruptions coloured by IS subgroup)

phyloI_IS_genedis_num <- phyloI_IS_genedis_filter %>% group_by(Round,Region,IS) %>% summarize(n=n()) %>% mutate(freq=n/sum(n)) 
newtab <- merge(phyloI_IS_genedis_num, phyloI_IS_genedis[, c("Round", "right_description")], by="Round")

newtab$IS <- factor(newtab$IS, levels = c("ISBma3", "ISRso7", "IS401","ISButh1","IS1021","IS1405","IS1420","IS1421","ISRso1","ISRso18","Unknown_IS5","ISRso17"), labels = c("ISBma3 (IS110)", "ISRso7 (IS256)", "IS401 (IS3)","ISButh1 (IS3)","IS1021 (IS5)","IS1405 (IS5)","IS1420 (IS5)","IS1421 (IS5)","ISRso1 (IS5)","ISRso18 (IS5)","Unknown (IS5)","ISRso17 (IS701)"))

nb.cols <- 12
mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols)

ggplot(data=newtab, aes(x=reorder(Round, -n),fill=IS))+
  geom_bar(aes(y=n), stat="identity", position = "dodge")+
  theme_bw()+
  ylim(0,59)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title=element_text(size=12,face="bold")) + 
  scale_x_discrete(breaks = newtab$Round, labels = newtab$right_description)+
  scale_fill_manual(values=mycolors)+
  xlab("Gene interrupted") +
  ylab("IS insertion frequency")+
  labs(fill = "IS subgroup") +
  theme(legend.title = element_text(colour="black", size=11,face="bold"),legend.text = element_text(size=10)) +
  theme(plot.margin=margin(0,0,0,90))+
  theme(axis.text.x = element_text(angle = 55, hjust = 1))+
  facet_grid(~ Region, scales = "free", space = "free",labeller = as_labeller(region_names))


## Figure 4B (Phylotype IIB gene disruptions)

x <- ladderize(IIB, right=FALSE)
is_tip <- x$edge[,2] <= length(x$tip.label)

ordered_tips <- x$edge[is_tip, 2]
ordered_tips_df <- data.frame(x$tip.label[ordered_tips])
ordered_tips_df$x.tip.label.ordered_tips. <- substr(ordered_tips_df$x.tip.label.ordered_tips., 0, 5)
ordered_tipsdf_list <- ordered_tips_df$x.tip.label.ordered_tips.
ordered_tipsdf_listrev <- rev(ordered_tipsdf_list)


## Phylotype IIB chromosome disruptions

phyloIIB_chr_IS_genedis <- IS_info[(IS_info$Phylotype=="IIB" & IS_info$Region=="Chr" & IS_info$gene_interruption=="TRUE"),]

genedis <- unique(phyloIIB_chr_IS_genedis$Isolate)
unique_isolates <- unique(phyloIIB_chr_IS$Isolate)

missing_vals <- tibble(setdiff(unique_isolates, genedis))
names(missing_vals)[1] <- "Isolate"

phyloIIB_chr_IS_genedis_edit <- bind_rows(phyloIIB_chr_IS_genedis,missing_vals)

phyloIIB_chr_IS_genedis_order <- phyloIIB_chr_IS_genedis_edit %>%
  mutate(Isolate =  factor(Isolate, levels = ordered_tipsdf_list)) %>%
  arrange(Isolate) 

phyloIIB_ispos_chr_genedis_fam <- ggplot(data=phyloIIB_chr_IS_genedis_order, aes(x=x, y=Isolate,colour=IS_family))+
  geom_point(size=1)+
  theme_bw()+
  theme(axis.ticks.y = element_blank()) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title=element_text(size=12,face="bold"),axis.text.y=element_blank(),axis.title.y=element_blank()) + 
  scale_colour_manual(values=IS_fam_cols,name="IS family")+
  xlab("Chromosome position (bp)") +
  theme(legend.title = element_text(colour="black", size=12,face="bold"),legend.text = element_text(size=12)) +
  theme(plot.margin=margin(0,0,0,10))+
  guides(color = guide_legend(
    override.aes=list(shape = 15, size=2)))

plot_grid(IIB_tree, phyloIIB_ispos_chr_genedis_fam, align = "h", ncol=2, rel_widths = c(0.25,0.75))


## Phylotype IIB megaplasmid disruptions

phyloIIB_mp_IS_genedis <- IS_info[(IS_info$Phylotype=="IIB" & IS_info$Region=="MP" & IS_info$gene_interruption=="TRUE"),]

genedis <- unique(phyloIIB_mp_IS_genedis$Isolate)
unique_isolates <- unique(phyloIIB_mp_IS$Isolate)

missing_vals <- tibble(setdiff(unique_isolates, genedis))
names(missing_vals)[1] <- "Isolate"

phyloIIB_mp_IS_genedis_edit <- bind_rows(phyloIIB_mp_IS_genedis,missing_vals)

phyloIIB_mp_IS_genedis_order <- phyloIIB_mp_IS_genedis_edit %>%
  mutate(Isolate =  factor(Isolate, levels = ordered_tipsdf_list)) %>%
  arrange(Isolate) 

phyloIIB_ispos_mp_genedis_fam <- ggplot(data=phyloIIB_mp_IS_genedis_order, aes(x=x, y=Isolate,colour=IS_family))+
  geom_point(size=1)+
  theme_bw()+
  theme(axis.ticks.y = element_blank()) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title=element_text(size=12,face="bold"),axis.text.y=element_blank(),axis.title.y=element_blank()) + 
  scale_colour_manual(values=IS_fam_cols,name="IS family")+
  xlab("Megaplasmid position (bp)") +
  theme(legend.title = element_text(colour="black", size=12,face="bold"),legend.text = element_text(size=12)) +
  theme(plot.margin=margin(0,0,0,10))+
  guides(color = guide_legend(
    override.aes=list(shape = 15, size=2)))

plot_grid(IIB_tree, phyloIIB_ispos_mp_genedis_fam, align = "h", ncol=2, rel_widths = c(0.25,0.75))


# Figure 4D (Phylotype IIB IS disruptions > 5 isolates coloured by IS family)

phyloIIB_IS_genedis <- IS_info[(IS_info$Phylotype=="IIB" & IS_info$gene_interruption=="TRUE"),]

phyloIIB_IS_genedis_filter <- phyloIIB_IS_genedis %>% 
  group_by(Round) %>% 
  mutate(freq = n()) %>% 
  ungroup() %>% 
  filter(freq > 5) %>%
  select(-freq)

phyloIIB_IS_genedis_num <- phyloIIB_IS_genedis_filter %>% group_by(Round,Region,IS) %>% summarize(n=n()) %>% mutate(freq=n/sum(n)) 
newtab <- merge(phyloIIB_IS_genedis_num, phyloIIB_IS_genedis[, c("Round", "right_description")], by="Round")

ggplot(data=newtab, aes(x=reorder(Round, -n),fill=IS))+
  geom_bar(aes(y=n), stat="identity", position = "dodge")+
  theme_bw()+
  ylim(0,269)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title=element_text(size=12,face="bold")) + 
  #scale_fill_manual(values=IS_fam_cols,name="IS family")+
  scale_x_discrete(breaks = newtab$Round, labels = newtab$right_description)+
  xlab("Gene interrupted") +
  ylab("IS insertion frequency")+
  theme(legend.title = element_text(colour="black", size=12,face="bold"),legend.text = element_text(size=12)) +
  theme(plot.margin=margin(0,0,0,90))+
  guides(color = guide_legend(
    override.aes=list(shape = 15, size=2)))+
  theme(axis.text.x = element_text(angle = 55, hjust = 1))+
  facet_grid(~ Region, scales = "free", space = "free",labeller = as_labeller(region_names))

