library(readxl)
library(plotly)
library(gapminder)
library(ggplot2)
library(ggfortify)
library(reshape2)
library(PCAtools)
library(factoextra)

dat <- read_excel("data/3rd_data_shortlisted_new.xlsx")
dat <- as.data.frame(dat)
datTraits <- read_excel("data/3rd_data_shortlisted_trait_new.xlsx")
var_names <- variable.names(dat)

label <- c(1:11, 13:20)
rownames(dat) <- label

# PCA on infected mice
pca_res <- prcomp(dat[-c(1:5),], center = T, scale. = TRUE)
summary(pca_res)
str(pca_res)

basis_new <- pca_res[["rotation"]][,1:2]
dat_new <- pca_res[["x"]][,1:2]

# visualize PCA proportion
fviz_eig(pca_res) # num of pc = 6
group <- c(0,0,1,0,0,0,0,0,0,1,0,1,1,1)
group_color <- ifelse(group=='0', 'green', 'red')

fviz_pca_ind(pca_res,
             colour = group_color,
             habillage=group,
             addEllipses=T,
             palette = c('black', 'red'),
             repel = TRUE     # Avoid text overlapping
) +
  geom_point(shape = 21, fill = group_color, size = 6.5)


fviz_pca_var(pca_res,
             geom = c("arrow", "text"),
             repel = TRUE     # Avoid text overlapping
)

# Project the clinical traits to PCA space
quanti.sup <- datTraits
quanti.coord <- cor(quanti.sup[-c(1:5),], pca_res$x)
p <- fviz_pca_var(pca_res,repel = TRUE, alpha.var = 0.6, label = F, labelsize=3)
fviz_add(p, quanti.coord, label = F, color ="blue", geom="arrow",
         linetype = 'solid', repel = TRUE, labelsize=7)

# hierarchical plot of variables and clinical traits

var_and_weight <- cbind(dat, datTraits)
var_and_weight0 <- t(var_and_weight[-c(1:5),])
dd <- get_dist(var_and_weight0, method = "pearson")
hc <- hclust(dd, method ='complete')
plot(hc)

# divided into three clusters
rect.hclust(hc , k = 3, border = 2:6)
abline(h = 3, col = 'red')

