library(patchwork)

# 1. Add CHELSA to physico-chemical variables
meta <- read.table("data/metadata.tsv", header = TRUE, row.names = 1, sep = "\t")
coords <- data.frame(lon = meta$longitude, lat = meta$latitude, row.names = meta$sample)
chelsa <- list.files(path = "data/chelsa", pattern = '.tif$', full.names = T)
names(chelsa) <- c("bio1", "bio12")
chelsa_raster <- raster::stack(chelsa)
points <- sp::SpatialPoints(coords, proj4string = chelsa_raster@crs)
chelsa_val <- raster::extract(chelsa_raster, points)
df <- cbind.data.frame(sp::coordinates(points), chelsa_val)
df$sample <- rownames(df)
meta$id <- rownames(meta)
meta.new <- meta
meta.new$id <- rownames(meta)
meta.new <- meta.new |> dplyr::left_join(df)

# 2. Perform Principal Component Analysis
pca.data <- meta.new |> 
  dplyr::select(id, carbon, nitrogen, c_n, bio1, bio12) |> 
  dplyr::select(-c(id)) |>
  dplyr::rename(MAT = bio1, MAP = bio12, N = nitrogen, OrgC = carbon,
                `C:N`=c_n)

rownames(pca.data) <- meta.new$id
res.pca <- prcomp(pca.data, scale. = T)

# Extract PCA results
pca_ind <- as.data.frame(res.pca$x) # Individual coordinates
pca_var <- as.data.frame(res.pca$rotation) # Variable loadings

# Add grouping information from metadata
pca_ind$phyto <- as.factor(meta.new$phyto)

# Plot using ggplot2
p1 <- ggplot2::ggplot(pca_ind, ggplot2::aes(x = PC1, y = PC2, color = phyto)) +
  ggplot2::geom_point(size = 2.5) +  # Individuals
  ggplot2::geom_segment(data = pca_var, ggplot2::aes(x = 0, y = 0, xend = PC1 * 5, yend = PC2 * 5), 
               arrow = ggplot2::arrow(length = ggplot2::unit(0.2, "cm")), color = "grey") +  # Variables
  ggrepel::geom_text_repel(data = pca_var, ggplot2::aes(x = PC1 * 5, y = PC2 * 5, label = rownames(pca_var)), 
                  size = 4, color = "grey") +  # Variable name
  ggplot2::labs(x = "PC1 (55.18 %)", y = "PC2 (24.25 %)", title = NULL,  color = "Phytogeographic regions") +
  theme_Publication()

export_plot(p1, "fig/Figure2.png", type="png", width = 8, height = 7)

# 3. Pariwise ADONIS analysis
kdata <- pca.data

# Separate environmental data
kdata.env <- meta.new |>
  dplyr::select(sample, age_category, phyto, forest)
kdata.env$age_category <- as.factor(kdata.env$age_category)
kdata.env$phyto <- as.factor(kdata.env$phyto)

ado.phyto <- metagMisc::adonis_pairwise(
  x = kdata.env, 
  dd = vegan::vegdist(kdata), 
  permut = 9999, 
  group.var = "phyto"
  )
ado.phyto$Adonis.tab
ado.phyto$Betadisper.tab

# 4. Carbon and nitrogen plots
cn.data <- data.frame(carbon = meta$carbon,
                      nitrogen = meta$nitrogen, 
                      cn = meta$c_n,
                      age_cat = meta$age_category, 
                      phyto = meta$phyto, 
                      crop_type=meta$crop_type, 
                      row.names = rownames(meta)
)

desired_order <- list("Evergreen", "Semi-deciduous", "Dry")

cn1 <- ggplot2::ggplot(cn.data, ggplot2::aes(x = phyto, y = carbon)) +
  ggplot2::geom_boxplot() +
  ggplot2::geom_jitter() +
  ggpubr::stat_compare_means(label.y = 5) +
  theme_Publication() +
  ggplot2::labs(x = NULL, y = "Organic Carbon")

cn1$data$phyto <- factor(cn1$data$phyto, levels=desired_order)

cn2 <- ggplot2::ggplot(cn.data, aes(x = phyto, y = nitrogen)) +
  ggplot2::geom_boxplot() +
  ggplot2::geom_jitter() +
  ggpubr::stat_compare_means(method = "wilcox.test", comparisons = list(c("Evergreen", "Semi-deciduous"), c("Evergreen", "Dry"), c("Semi-deciduous", "Dry"))) +
  ggpubr::stat_compare_means(label.y = 0.37) +
  theme_Publication() +
  ggplot2::labs(x = NULL, y = "Nitrogen")

cn2$data$phyto <- factor(cn2$data$phyto, levels = desired_order)

cn3 <- ggplot2::ggplot(cn.data, aes(x = phyto, y = cn)) +
  ggplot2::geom_boxplot() +
  ggplot2::geom_jitter() +
  ggpubr::stat_compare_means(method = "wilcox.test", comparisons = list(c("Evergreen", "Semi-deciduous"), c("Evergreen", "Dry"), c("Semi-deciduous", "Dry"))) +
  ggpubr::stat_compare_means(label.y = 20) +
  theme_Publication() +
  ggplot2::labs(x = NULL, y = "Carbon nitrogen ratio")

cn3$data$phyto <- factor(cn3$data$phyto, levels = desired_order)

cn <- cn1 + cn2 + cn3
export_plot(cn, "fig/Figure3.png", type = "png", width = 13, height = 8)

# 5. Phyloseq import
vtable <- read.csv("data/amf_table.csv", header=T, row.names = 1)
vtax <- read.table("data/18s_tax.tsv", header=T, row.names=1, sep="\t")
met <- read.table("data/metadata.tsv", header=T, row.names = 1, sep="\t")
vps <- phyloseq::phyloseq(phyloseq::otu_table(vtable, taxa_are_rows = T),
                          phyloseq::tax_table(as.matrix(vtax)),
                          phyloseq::sample_data(met))
btable <- read.csv("data/asv_table.csv", row.names = 1)
btax <- read.table("data/v3v4_tax.tsv", header=T, row.names=1, sep="\t")
bps <- phyloseq::phyloseq(phyloseq::otu_table(btable, taxa_are_rows = T),
                          phyloseq::tax_table(as.matrix(btax)),
                          phyloseq::sample_data(met))

# 6. iNEXT plot
bac.prepare.inext <- metagMisc::prepare_inext(bps)
amf.prepare.inext <- metagMisc::prepare_inext(vps)
bac.inext <- iNEXT::iNEXT(bac.prepare.inext, q=c(0,1,2))
amf.inext <- iNEXT::iNEXT(amf.prepare.inext, q=c(0,1,2))

pi1 <- iNEXT::ggiNEXT(bac.inext, type=1, facet.var="Order.q") + 
  ggplot2::facet_wrap(~Order.q, scales="free")
pi2 <- iNEXT::ggiNEXT(amf.inext, type=1, facet.var="Order.q") + 
  ggplot2::facet_wrap(~Order.q, scales="free")

export_plot(pi1, "fig/FigureS1.png", type="png", width=21, height=13)
export_plot(pi2, "fig/FigureS2.png", type="png", width=21, height=13)

# 7. Share venn diagram
get_venn_data <- function(phys, sample_factor) {
  # Convert to presence-absence data
  ps.pa <- phyloseq::transform_sample_counts(phys, function(x) ifelse(x > 0, 1, 0))
  
  # Extract metadata
  metadata <- phyloseq::sample_data(ps.pa)
  
  # Get unique levels of the sample factor
  groups <- unique(metadata[[sample_factor]])
  
  # Create a list for Venn diagram
  venn.list <- list()
  
  for (group in groups) {
    # Subset samples for each group
    sample_names <- rownames(metadata[metadata[[sample_factor]] == group, ])
    ps_subset <- phyloseq::prune_samples(sample_names, ps.pa)
    
    # Get ASVs/OTUs present in the group
    present_ASVs <- phyloseq::taxa_names(phyloseq::prune_taxa(phyloseq::taxa_sums(ps_subset) > 0, ps_subset))
    
    venn.list[[group]] <- present_ASVs
  }
  return(venn.list)
}

venn.data.vps <- get_venn_data(vps, sample_factor="phyto")
g1 <- ggvenn::ggvenn(venn.data.vps, fill_color = c("#00A087FF","#4DBBD5FF", "#E64B35FF"))

venn.data.bps <- get_venn_data(bps, sample_factor="phyto")
g2 <- ggvenn::ggvenn(venn.data.bps, fill_color = c("#E64B35FF","#4DBBD5FF", "#00A087FF"))

g <- g1 + g2 + patchwork::plot_annotation(tag_levels='A')
export_plot(g, "fig/Figure4.png", type = "png", width = 13, height = 8)

# 8. Representativeness plots
compute_prevalence <- function(physeq, rank) {
  # Aggregate at the specified taxonomic rank
  physeq_ranked <- phyloseq::tax_glom(physeq, taxrank = rank, NArm = TRUE)
  
  # Extract OTU table as a matrix
  otu_mat <- as(phyloseq::otu_table(physeq_ranked), "matrix")
  
  # Compute prevalence for each taxon
  prevalence_df <- data.frame(
    Taxon = phyloseq::tax_table(physeq_ranked)[, rank],
    Prevalence = rowSums(otu_mat > 0) / ncol(otu_mat)  # Fraction of non-zero samples
  )
  
  return(prevalence_df)
}

# Function to compute average relative abundance
compute_avg_rel_abundance <- function(physeq, rank) {
  # Aggregate at taxonomic rank
  physeq_ranked <- phyloseq::tax_glom(physeq, taxrank = rank, NArm = TRUE)
  
  # Transform to relative abundance
  physeq_rel <- phyloseq::transform_sample_counts(physeq_ranked, function(x) x / sum(x))
  
  # Extract OTU table
  otu_mat_rel <- as(phyloseq::otu_table(physeq_rel), "matrix")
  
  # Compute mean relative abundance per taxon
  rel_abundance_df <- data.frame(
    Taxon = phyloseq::tax_table(physeq_ranked)[, rank],
    Avg_Rel_Abundance = rowMeans(otu_mat_rel, na.rm = TRUE)
  )
  
  return(rel_abundance_df)
}

amf.prev.abun <- merge(
  compute_prevalence(vps, "Genus"), 
  compute_avg_rel_abundance(vps, "Genus"), 
  by = "Genus"
  )
amf.prev.abun <- amf.prev.abun |> 
  dplyr::left_join(as.data.frame(phyloseq::tax_table(vps))) |> 
  dplyr::select(Genus, Family, Prevalence, Avg_Rel_Abundance) |> 
  unique()

pa1 <- ggplot2::ggplot(amf.prev.abun, ggplot2::aes(x = -Prevalence, y = reorder(Genus, Avg_Rel_Abundance), fill = Family)) + 
  ggplot2::geom_bar(stat = "identity") + 
  ggsci::scale_fill_npg() + 
  ggplot2::scale_y_discrete(position = "right") +
  ggplot2::scale_x_continuous(labels = c(100, 75, 50, 25, 0), position = "top") + 
  ggplot2::labs(x = "Proportion of sampling sites (%)", y = NULL) + 
  theme_Publication() + 
  ggplot2::theme(axis.text.y = element_text(face = "italic"))

pa2 <- ggplot2::ggplot(amf.prev.abun, ggplot2::aes(x = Avg_Rel_Abundance * 100, y = reorder(Genus, Avg_Rel_Abundance), fill = Family)) + 
  ggplot2::geom_bar(stat = "identity") + 
  ggsci::scale_fill_npg() + 
  ggplot2::scale_x_continuous(position = "top") + 
  ggplot2::labs(x = "Average relative abundance of VTs (%)", y = NULL) + 
  theme_Publication() + 
  ggplot2::guides(fill = "none") +
  ggplot2::theme(axis.text.y = element_blank())

pa3 <- pa1 + pa2 + patchwork::plot_layout(guides = 'collect') & 
  ggplot2::theme(legend.position='bottom')

export_plot(pa3, "fig/Figure5.png", type = "png", width = 13, height = 8)

bac.prev <- merge(
  compute_prevalence(bps, "Genus"), 
  compute_avg_rel_abundance(bps, "Genus"), 
  by = "Genus"
  )
bac.prev <- bac.prev |> 
  dplyr::left_join(as.data.frame(phyloseq::tax_table(bps))) |> 
  dplyr::select(Genus, Phylum, Prevalence, Avg_Rel_Abundance) |> 
  unique()

# selecting the 15 top
bac.prev.abun <- bac.prev |> dplyr::arrange(Avg_Rel_Abundance) |> dplyr::top_n(20)

phyla_levels <- levels(factor(bac.prev.abun$Phylum))

# Get the color palette used in pa4
phyla_colors <- ggsci::pal_npg("nrc")(length(phyla_levels))
names(phyla_colors) <- phyla_levels 

pa4 <- ggplot2::ggplot(bac.prev.abun, ggplot2::aes(x = -Prevalence, y = reorder(Genus, Avg_Rel_Abundance), fill = Phylum, order = Phylum)) + 
  ggplot2::geom_bar(stat = "identity") + 
  ggsci::scale_fill_npg() + 
  ggplot2::scale_y_discrete(position = "right") +
  ggplot2::scale_x_continuous(labels = c(100, 75, 50, 25, 0), position = "top") + 
  ggplot2::labs(x = "Proportion of sampling sites (%)", y = NULL) + 
  theme_Publication() + 
  ggplot2::theme(axis.text.y = element_blank())

pa5 <- ggplot2::ggplot(bac.prev.abun, aes(x = Avg_Rel_Abundance * 100, y = reorder(Genus, Avg_Rel_Abundance), fill = Phylum)) + 
  ggplot2::geom_bar(stat = "identity") + 
  ggplot2::scale_fill_manual(values = phyla_colors) +  # Use the extracted colors
  ggplot2::scale_x_continuous(position = "top") + 
  ggplot2::labs(x = "Average relative abundance of ASVs (%)", y = NULL) + 
  theme_Publication() + 
  ggplot2::guides(fill = "none") +
  ggplot2::theme(axis.text.y = element_text(face = "italic"))

pa6 <- pa4 + pa5 + patchwork::plot_layout(guides = 'collect') & 
  ggplot2::theme(legend.position = 'bottom')

export_plot(pa6, "fig/Figure6.png", type = "png", width = 13, height = 8)

# 9. Core members
vt.abun <- phyloseq::otu_table(vps)
vt.abun.rel <- sweep(vt.abun, 2, colSums(vt.abun), "/")
vtt <- vps
phyloseq::otu_table(vtt) <- phyloseq::otu_table(vt.abun.rel, taxa_are_rows = TRUE)
core.amf <- phyloseq::prune_taxa(microbiome::core_members(vps, prevalence=0.8), vtt)
core.amf.abund.sum.all <- apply(phyloseq::otu_table(core.amf), 1, mean)
core.amf.abund.df.all <- data.frame(taxon=names(core.amf.abund.sum.all), mean.rel.abun = core.amf.abund.sum.all)
vt.tax.df <- as.data.frame(phyloseq::tax_table(vtt))
vt.tax.df$taxon <- rownames(vt.tax.df)
core.amf.abund.df.all <- merge(core.amf.abund.df.all, vt.tax.df, by="taxon", all.x=T)
core.amf.abund.df.all <- core.amf.abund.df.all[, c("taxon", "Genus", "Species", "mean.rel.abun")]
core.amf.abund.df.all <- core.amf.abund.df.all[order(-core.amf.abund.df.all$mean.rel.abun),]
core.amf.abund.df.all |> dplyr::group_by(Species) |> dplyr::summarise(mm = sum(mean.rel.abun))

bac.abun <- phyloseq::otu_table(bps)
bac.abun.rel <- sweep(bac.abun, 2, colSums(bac.abun), "/")
btt <- bps
phyloseq::otu_table(btt) <- phyloseq::otu_table(bac.abun.rel, taxa_are_rows = TRUE)
core.bac <- phyloseq::prune_taxa(microbiome::core_members(bps, prevalence=0.8), btt)
core.bac.abund.sum.all <- apply(phyloseq::otu_table(core.bac), 1, mean)
core.bac.abund.df.all <- data.frame(taxon=names(core.bac.abund.sum.all), mean.rel.abun = core.bac.abund.sum.all)
bac.tax.df <- as.data.frame(phyloseq::tax_table(btt))
bac.tax.df$taxon <- rownames(bac.tax.df)
core.bac.abund.df.all <- merge(core.bac.abund.df.all, bac.tax.df, by="taxon", all.x=T)
core.bac.abund.df.all <- core.bac.abund.df.all[, c("taxon", "Genus", "mean.rel.abun")]
core.bac.abund.df.all <- core.bac.abund.df.all[order(-core.bac.abund.df.all$mean.rel.abun),]
core.bac.abund.df.all |> dplyr::group_by(Genus) |> dplyr::summarise(mm = sum(mean.rel.abun))


# 10. Alpha diversity
## Bacteria
bcom <- metagMisc::phyloseq_to_MetaCommunity(bps)
richness.bac <- entropart::DivPart(q = 0, MC = bcom)$CommunityAlphaDiversities
shannon.bac <- entropart::DivPart(q = 1, MC = bcom)$CommunityAlphaDiversities
simpson.bac <- entropart::DivPart(q = 2, MC = bcom)$CommunityAlphaDiversities

## AMF
vcom <- metagMisc::phyloseq_to_MetaCommunity(vps)
richness.amf <- entropart::DivPart(q = 0, MC = vcom)$CommunityAlphaDiversities
shannon.amf <- entropart::DivPart(q = 1, MC = vcom)$CommunityAlphaDiversities
simpson.amf <- entropart::DivPart(q = 2, MC = vcom)$CommunityAlphaDiversities

alpha.div <- cbind.data.frame(richness.bac, shannon.bac, simpson.bac, 
                              richness.amf, shannon.amf, simpson.amf)
alpha.div$sample <- rownames(alpha.div)

alpha.div <- alpha.div |> 
  tidyr::pivot_longer(c(richness.bac:simpson.amf),
                      names_to = "index",
                      values_to = "estimate") |>
  dplyr::mutate(organism = dplyr::case_when(stringr::str_detect(index, "bac") ~ "Bacteria", .default = "AMF"),
                index = dplyr::recode(index, richness.bac = "Richness", richness.amf = "Richness", shannon.bac = "Shannon", shannon.amf = "Shannon", simpson.bac = "Simpson", simpson.amf = "Simpson"))

met$xid <- rownames(met)
alpha.div <- alpha.div |> 
  dplyr::left_join(met, by = c("sample"="xid")) |>
  dplyr::select(sample, index, estimate, organism, age, age_category, phyto, forest)

summary_stats <- alpha.div |>
  dplyr::group_by(phyto, organism, index) |>
  dplyr::summarise(
    mean_estimate = mean(estimate, na.rm = TRUE),
    sd_estimate = sd(estimate, na.rm = TRUE)
  )

# 11. Bray-Curtis dissimilarity
c1 <- list(
  c("1-10", "11-20"),
  c("1-10", "21-30"),
  c("1-10", "30-sup"),
  c("1-10", "OF"),
  c("11-20", "21-30"),
  c("11-20", "30-sup"),
  c("11-20", "OF"),
  c("21-30", "30-sup"),
  c("21-30", "OF"),
  c("30-sup", "OF")
)

c2 <- list(
  c("Evergreen", "Dry"),
  c("Evergreen", "Semi-deciduous"),
  c("Dry", "Semi-deciduous")
)

generate_beta_diversity_plot <- function(ps,  groups_column, comp_clusters, desired_order, label, ref, col) {
  require(phyloseq)
  require(ggplot2)
  require(reshape2)
  require(ggpubr)
  require(ggsci)
  
  # Compute dissimilarity matrix using mult_dissim_func
  dis_matrix <- metagMisc::mult_dissim(ps, average=T)
  dis_matrix <- as.matrix(dis_matrix)
  
  # Group samples based on the provided column
  sub_dist <- list()
  groups_all <- sample_data(ps[[1]])[[groups_column]]
  groups_all <- as.factor(groups_all)
  
  for (group in levels(groups_all)) {
    row_group <- which(groups_all == group)
    sample_group <- sample_names(ps[[1]])[row_group]
    sub_dist[[group]] <- dis_matrix[sample_group, sample_group]
    sub_dist[[group]][!lower.tri(sub_dist[[group]])] <- NA
  }
  
  # Melt the list of matrices into a data frame
  braygroups <- melt(sub_dist)
  df.bray <- braygroups[complete.cases(braygroups), ]
  df.bray$L1 <- factor(df.bray$L1, levels = names(sub_dist))
  
  # Create the plot
  p <- ggplot2::ggplot(df.bray, aes(x = L1, y = value, colour = L1)) +
    ggplot2::geom_jitter() +
    ggplot2::geom_boxplot(alpha = 0.6) +
    theme_Publication() +
    ggplot2::scale_colour_manual(values=col) +
    ggpubr::stat_compare_means(ref.group = ref, method = "wilcox.test",label = "p.signif") +
    ggplot2::labs(x = NULL, y = paste("Beta Diversity by ", label), color = label)
  
  # Reorder levels and sort data
  p$data$L1 <- factor(p$data$L1, levels = desired_order)
  p$data <- p$data[order(p$data$L1), ]
  
  return(p)
}

d1 <- list("Evergreen", "Semi-deciduous", "Dry")
d2 <- list("1-10", "11-20", "21-30", "30-sup", "OF")
bbp <- metagMisc::phyloseq_mult_raref(bps, multithread = T, seeds=1:1000)
vvp <- metagMisc::phyloseq_mult_raref(vps, multithread = T, seeds=1:1000)

bet1 <- generate_beta_diversity_plot(bbp, "phyto", c2, d1, "Phytogeographic area", "Evergreen", c("#00A087FF","#4DBBD5FF", "#E64B35FF"))
bet2 <- generate_beta_diversity_plot(bbp, "age_category", c1, d2, "Age category", "OF", col=c("#E64B35FF", "#4DBBD5FF", "#3C5488FF", "#F39B7FFF", "#00A087FF"))

bet3 <- generate_beta_diversity_plot(vvp, "phyto", c2, d1, "Phytogeographic area", "Evergreen", col=c("#00A087FF","#4DBBD5FF", "#E64B35FF"))
bet4 <- generate_beta_diversity_plot(vvp, "age_category", c1, d2, "Age category", "OF", col=c("#E64B35FF", "#4DBBD5FF", "#3C5488FF", "#F39B7FFF", "#00A087FF"))

bet.bac <- bet1 + bet2 + plot_annotation(tag_levels = 'A')
bet.amf <- bet3 + bet4 + plot_annotation(tag_levels = 'A')

export_plot(bet.amf, "fig/Figure7.png", type = "png", width = 13, height = 8)
export_plot(bet.bac, "fig/Figure8.png", type = "png", width = 13, height = 8)

# 12. Indicator species
# Create nitrogen and carbon groups
km.nitro <- kmeans(vps@sam_data$nitrogen, centers = 3)
km.carbon <- kmeans(vps@sam_data$carbon, centers = 3)

cn.env <- data.frame(nitrogen = vps@sam_data$nitrogen, carbon = vps@sam_data$carbon)
rownames(cn.env) <- rownames(vps@sam_data)

cn.env$clusn <- km.nitro$cluster
cn.env$clusc <- km.carbon$cluster

cluster_names <- cn.env |>
  dplyr::group_by(clusn) |>
  dplyr::summarize(Mean_Nitrogen = mean(nitrogen)) |>
  dplyr::arrange(Mean_Nitrogen) |>
  dplyr::mutate(nitrogen_group = c("Low", "Moderate", "High"))

cn.env$nitrogen_group <- ifelse(cn.env$clusn == 1, "Low", ifelse(cn.env$clusn == 2, "Moderate", "High"))

cluster_names <- cn.env |>
  dplyr::group_by(clusc) |>
  dplyr::summarize(Mean_carbon = mean(carbon)) |>
  dplyr::arrange(Mean_carbon) |>
  dplyr::mutate(carbon_group = c("Low", "Moderate", "High"))

cn.env$carbon_group <- ifelse(cn.env$clusc == 2, "Low", ifelse(cn.env$clusc == 3, "Moderate", "High"))

# indicator species analysis
indval.amf.nitrogen <- indicspecies::multipatt(
  as.data.frame(t(vps@otu_table)),
  cn.env$nitrogen_group,
  control = permute::how(nperm=9999),
  func="r.g",
  allow.negative = T
)

indval.amf.carbon <- indicspecies::multipatt(
  as.data.frame(t(vps@otu_table)),
  cn.env$carbon_group,
  func="r.g",
  allow.negative = T,
  control = permute::how(nperm=9999)
)

indval.bac.nitrogen <- indicspecies::multipatt(
  as.data.frame(t(core.bac@otu_table)),
  cn.env$nitrogen_group,
  func="r.g",
  allow.negative = T,
  control = permute::how(nperm=9999)
)

indval.bac.carbon <- indicspecies::multipatt(
  as.data.frame(t(core.bac@otu_table)),
  cn.env$carbon_group,
  func="r.g",
  allow.negative = T,
  control = permute::how(nperm=9999)
)

# 13. Variance Partition
# Data preparation
varbac <- t(cbind.data.frame(richness.bac, shannon.bac, simpson.bac))
varbac <- varbac[,order(colnames(varbac))]
rownames(varbac) <- c("Richness", "Shannon", "Simpson")

varamf <- t(cbind.data.frame(richness.amf, shannon.amf, simpson.amf))
varamf <- varamf[,order(colnames(varamf))]
rownames(varamf) <- c("Richness", "Shannon", "Simpson")

env <- meta.new |>
  dplyr::rename(MAT=bio1, MAP=bio12, C=carbon, N=nitrogen, `C:N`=c_n) |>
  dplyr::select(C, N, `C:N`, MAT, MAP)
rownames(env) <- rownames(meta)
env <- as.data.frame(scale(env))
env <- env[order(row.names(env)),]

# Core variance partitioning
form <- ~ C + N + C:N + MAT + MAP
vp1 <- variancePartition::fitExtractVarPartModel(varamf, form, env)
vp2 <- variancePartition::fitExtractVarPartModel(varbac, form, env)

# Visualisation
pvp1 <- variancePartition::plotPercentBars(
  variancePartition::sortCols(vp1), 
  col = c("#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF", "#F39B7FFF", "grey"), 
  show.legend = FALSE
  ) +
  variancePartition::plotVarPart(
    variancePartition::sortCols(vp1), 
    col = c("#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF", "#F39B7FFF", "grey"), 
    show.legend = FALSE
  ) + 
  ggplot2::labs(x = "Attributes") + 
  ggplot2::theme(legend.position = "none")

pvp2 <- variancePartition::plotPercentBars(
  variancePartition::sortCols(vp2), 
  col = c("#4DBBD5FF", "#E64B35FF", "#00A087FF", "#F39B7FFF", "#3C5488FF", "grey"), 
  show.legend = FALSE
  ) +
  variancePartition::plotVarPart(
    variancePartition::sortCols(vp2), 
    col = c("#4DBBD5FF", "#E64B35FF", "#00A087FF", "#F39B7FFF", "#3C5488FF", "grey"), 
    show.legend = FALSE
  ) + 
  ggplot2::labs(x = "Attributes") +
  ggplot2::theme(legend.position = "none")

# Combine the plots and keep only the fill legend
pvp <- pvp1 / pvp2 + 
  patchwork::plot_annotation(tag_levels = "A") + 
  patchwork::plot_layout(guides = "collect")

export_plot(pvp, "fig/Figure9.png", type="png", width=13, height=10)
