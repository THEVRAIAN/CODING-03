library(tidyverse)
library(snpR)
library(rMVP)
library(bigmemory)
library(parallel)

# Configuración del entorno de trabajo
vcf_path <- "D:/2026/MAESTRIA 2026 II/DIVERSIDAD_GENÉTICA/GIT/CODING-03/chickens_small.vcf.gz"
out_dir <- "dataset/rmvp/chickens_gwas"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cores <- detectCores() - 1

# 1. Carga de datos genotípicos y definición de fenotipos
geno <- snpR::read_vcf(vcf_path)

sample_ids <- sample.meta(geno)$sample

# Estructura fenotípica base adaptada a las muestras del VCF
pheno <- data.frame(
  Taxa = sample_ids,
  Trait1 = rnorm(length(sample_ids), mean = 10, sd = 2) # Reemplazar con datos fenotípicos reales
)

geno
sample.meta(geno)
snp.meta(geno)

# 2. Análisis exploratorio inicial (PCA)
sample.meta(geno)$population <- c("PopA", "PopB")
pca <- plot_clusters(geno
                     , facets = "population"  
                     , check_duplicates = TRUE
                     , plot_type = "pca"
)


pca$plots$pca <- pca$plots$pca +
  scale_color_manual(values = c("PopA" = "green", "PopB" = "yellow", "PopC" = "green", "PopD" = "orange")) +
  scale_fill_manual(values = c("PopA" = "green", "PopB" = "yellow", "PopC" = "green", "PopD" = "orange"))
pca$plots$pca

# 3. Preparación de insumos para rMVP
#AJUSTES EN LA DATA 
# 1. Extraer los IDs reales usando la columna 'sampID' de los metadatos de snpR
ids_geno <- sample.meta(geno)$sampID

# 2. Verificar que ahora sí muestre los nombres de los pollos
print(head(ids_geno))

# 3. Construir el dataframe de fenotipos con los IDs correctos
pheno <- data.frame(
  Taxa = ids_geno,
  Trait1 = rnorm(length(ids_geno), mean = 10, sd = 2),
  Trait2 = rnorm(length(ids_geno), mean = 50, sd = 5),
  stringsAsFactors = FALSE
)

# 4. Guardar el archivo tabulado limpio
pheno_file <- file.path(out_dir, "chickens_pheno.txt")
write.table(pheno, pheno_file, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# 5. Ejecutar MVP.Data con los datos sincronizados
MVP.Data(fileVCF = uncompressed_vcf
         , filePhe = pheno_file
         , sep.phe = "\t"
         , header = TRUE
         , priority = "speed"
         , fileKin = FALSE
         , filePC = FALSE
         , out = file.path(out_dir, "gwas")
# 4. Cálculo de matriz de parentesco (Kinship)
MVP.Data.Kin(fileKin = TRUE
             , mvp_prefix = file.path(out_dir, 'gwas')
             , out = file.path(out_dir, 'gwas'))

# 5. Carga de archivos binarios transformados
phenotype <- read.table(paste0(out_dir, "/gwas.phe"), header = TRUE)
genotype <- attach.big.matrix(paste0(out_dir, "/gwas.geno.desc"))
mapping <- read.table(paste0(out_dir, "/gwas.geno.map"), header = TRUE)
kinship <- attach.big.matrix(paste0(out_dir, "/gwas.kin.desc"))

# 6. Ejecución de modelos GWAS (GLM, MLM y FarmCPU)
phenotype %>% names()

trait <- phenotype %>% 
  select(1, "Trait1")

imMVP <- MVP(
  phe = trait
  , geno = genotype
  , map = mapping
  , K = kinship
  , method = c("GLM", "MLM", "FarmCPU")
  , outpath = out_dir
  , file.output = TRUE
  , file.type = "csv"
  , ncpus = cores
)

gc()

# 7. Generación de reportes y visualizaciones
# Densidad de SNPs
MVP.Report(MVP = imMVP
           , plot.type = "d"
           , col = c("darkgreen", "yellow", "red")
           , file.type = "jpg", dpi = 300
           , outpath = out_dir
)

# Histograma de distribución fenotípica
MVP.Hist(phe = trait
         , outpath = out_dir
         , file.type = "jpg", breakNum = 18, dpi = 300)

# Gráfico de Manhattan
MVP.Report(imMVP
           , plot.type = "m"
           , multracks = TRUE
           , threshold = c(1e-6, 1e-4)
           , threshold.lty = c(1, 2)
           , threshold.lwd = c(1, 1)
           , threshold.col = c("black", "grey")
           , amplify = TRUE
           , bin.size = 1e6
           , chr.den.col = c("darkgreen", "yellow", "red")
           , signal.col = c("red", "green")
           , signal.cex = c(1, 1)
           , file.type = "jpg", memo = "", dpi = 300
           , outpath = out_dir
)

# Gráfico Q-Q
MVP.Report(imMVP
           , plot.type = "q"
           , col = c("dodgerblue1", "olivedrab3", "darkgoldenrod1")
           , threshold = 1e-6
           , signal.pch = 19
           , signal.cex = 1.5
           , signal.col = "red"
           , box = FALSE 
           , multracks = TRUE
           , file.type = "jpg", memo = "", dpi = 300
           , outpath = out_dir
)
