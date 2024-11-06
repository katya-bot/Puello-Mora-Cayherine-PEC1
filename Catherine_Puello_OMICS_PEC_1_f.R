# ANÁLISIS DE DATOS ÓMICOS
# CATHERINE PUELLO MORA
# DATASET ESCOJIDA: 2023-UGrX-4MetaboAnalystTutorial

###########################  1. CARGA DE DATOS  ################################

# primero cargo mis archivos y creo My_data para que guarde el contenido
# la data está sucia, en comentarios del archivo sale que necesita limpieza

My_data <- read.csv("C:/Users/Usuario/Downloads/D_omics/my_dataset_git/ST000002_AN000002_clean.csv", 
                    header = TRUE, sep = "\t", skip = 70)

# ahora exploro el contenido de mis archivos para ver que tipo de variables hay

head(My_data)
str(My_data)

# preparo y limpio mis datos según description.md del dataset 

# Renombro la primera columna como 'Samples' y las demás columnas 
# como 'Sample_' para mayor claridad
colnames(My_data)[1] <- "Samples"
colnames(My_data)[-1] <- paste0("Sample_", 1:(ncol(My_data) - 1))

# Verifico los nombres de las columnas y la estructura

colnames(My_data)
head(My_data)


is.na(My_data) 
sum(is.na(My_data))
summary(My_data)


## Defino grupos 

grupos <- c(rep("After", 6), rep("Before", 6))
print(grupos)

################ 2. Creación del contenedor SummarizedExperiment ###############

if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
  BiocManager::install("SummarizedExperiment")
}
library(SummarizedExperiment)

# Separo los datos de expresión y asigno nombres de fila a las muestras

expresion_data <- as.matrix(My_data[, -1])
rownames(expresion_data) <- My_data$Samples

# Creo un DataFrame con los datos de los grupos

grupo_info <- DataFrame(Samples = colnames(expresion_data), Groups = grupos)

# Creo el objeto SummarizedExperiment

se <- SummarizedExperiment(assays = list(counts = expresion_data), colData = grupo_info)

se



####################### 3. Análisis Exploratorio ###############################

# Resumen estadístico:

summary_stats <- summary(assay(se))
print(summary_stats)

# Histograma de los valores de expresión

hist(assay(se), main = "Distribución de Valores de Expresión", xlab = "Valores de Expresión", 
     breaks = 30, col = "lightblue")

# Transposición de los datos de expresión

expresion_t <- t(assay(se))

# Boxplot log-transformado por grupos (mis datos tienen muchos outliers)

boxplot(log10(expresion_t) ~ grupos, main = " Distribución de Expresión por Grupo (Trans Logarítmica)", 
        ylab = "Log10(Valores de Expresión)", xlab = "Grupo", col = c("lightblue", "lightgreen"))

# Gráfico de densidad

if (!requireNamespace("reshape2", quietly = TRUE)) install.packages("reshape2")
library(reshape2)

expresion_log <- log10(assay(se) + 1)
expresion_long <- melt(expresion_log)
colnames(expresion_long) <- c("Metabolite", "Sample", "Expression")
expresion_long$Group <- rep(grupo_info$Groups, each = nrow(expresion_log))

if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
library(ggplot2)

ggplot(expresion_long, aes(x = Expression, fill = Group)) +
  geom_density(alpha = 0.5) +
  labs(title = "Densidad de Expresión por Grupo (Transformación Logarítmica)", 
       x = "Log10(Valores de Expresión)", y = "Densidad") +
  theme_minimal() +
  scale_fill_manual(values = c("tomato", "lightgreen"))


################## 3.2 Análisis Estadístico ####################################

# Prueba de Shapiro-Wilk para cada muestra en cada grupo

shapiro_test_results <- apply(assay(se), 1, function(x) {
  before <- x[grupo_info$Groups == "Before"]
  after <- x[grupo_info$Groups == "After"]
  
  c(before = shapiro.test(before)$p.value, after = shapiro.test(after)$p.value)
})

shapiro_df <- as.data.frame(t(shapiro_test_results))
colnames(shapiro_df) <- c("Shapiro_Before", "Shapiro_After")
head(shapiro_df)

# como tengo multitud de datos y cada metabolito puede tener valors muy diferentes,
# Creo una Función para aplicar la prueba de significancia adecuada

# https://explainedstatistics.com/what-is-error-rate-in-statistics-a-comprehensive-guide/
# https://statisticsbyjim.com/hypothesis-testing/mann-whitney-u-test/

#  si ambos grupos tienen distribución normal: se usa prueba t

analizar_metabolito <- function(expresion, grupo, shapiro_before, shapiro_after) {
  if (shapiro_before > 0.05 && shapiro_after > 0.05) {
    test_result <- t.test(expresion[grupo == "Before"], expresion[grupo == "After"])
    return(c(P_value = test_result$p.value, Test_Type = "t-test"))
  } else {
    
    #  si al menos un grupo no tiene distribución normal: se usa prueba de Mann-Whitney
    test_result <- wilcox.test(expresion[grupo == "Before"], expresion[grupo == "After"])
    return(c(P_value = test_result$p.value, Test_Type = "Mann-Whitney"))
  }
}

# Aplico análisis estadístico:

analisis_estadistico <- function(se, grupos) {
  shapiro_test_results <- apply(assay(se), 1, function(x) {
    before <- x[grupos == "Before"]
    after <- x[grupos == "After"]
    c(Shapiro_Before = shapiro.test(before)$p.value, Shapiro_After = shapiro.test(after)$p.value)
  })
  
  shapiro_df <- as.data.frame(t(shapiro_test_results))
  stat_test_results <- apply(assay(se), 1, function(x, i) {
    analizar_metabolito(x, grupos, shapiro_df[i, "Shapiro_Before"], shapiro_df[i, "Shapiro_After"])
  }, i = rownames(shapiro_df))
  
  stat_test_df <- as.data.frame(t(stat_test_results))
  colnames(stat_test_df) <- c("P_value", "Test_Type")
  rownames(stat_test_df) <- rownames(shapiro_df)
  
  return(list(shapiro_df = shapiro_df, stat_test_df = stat_test_df))
}

# Ejecuto el análisis completo y lo visualizo:

resultado <- analisis_estadistico(se, grupo_info$Groups)
head(resultado$shapiro_df)
head(resultado$stat_test_df)

# Ajusto los valores de p

resultado$stat_test_df$Adjusted_P_value <- p.adjust(resultado$stat_test_df$P_value, method = "BH")

# filtrado con valor 0.1

metabolitos_significativos <- resultado$stat_test_df[resultado$stat_test_df$Adjusted_P_value < 0.1, ]
print(metabolitos_significativos)

############### PCA  ##########################################################3

# Filtro los datos significativos:

significant_metabolites <- rownames(metabolitos_significativos)
expresion_significativa <- assay(se)[significant_metabolites, ]


pca_result <- prcomp(t(expresion_significativa), scale = TRUE)

# gráfico: 

pca_df <- data.frame(PC1 = pca_result$x[, 1], PC2 = pca_result$x[, 2], Group = grupo_info$Groups)
ggplot(pca_df, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 4, alpha = 0.7) +
  labs(title = "PCA de Metabolitos Significativos") +
  scale_color_manual(values = c("tomato", "lightblue")) +
  theme_minimal()


#################### Análisis de Correlación y Redes Metabólicas ###############

# matriz de correlación y mapa de calor

correlation_matrix <- cor(t(expresion_significativa), method = "spearman")
head(correlation_matrix)

library(pheatmap)
pheatmap(correlation_matrix, color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
         main = "Mapa de Calor de la Correlación de Spearman")



# Filtrado de correlaciones

cor_matrix <- correlation_matrix
cor_matrix[abs(cor_matrix) < 0.7] <- 0 # así Solo conserva correlaciones fuertes

#  igraph y construcción de la red 

if (!requireNamespace("igraph", quietly = TRUE)) install.packages("igraph")
library(igraph)

cor_graph <- graph_from_adjacency_matrix(cor_matrix, mode = "undirected", weighted = TRUE, diag = FALSE)
cor_graph <- simplify(cor_graph)

# Visualización de la red 

plot(cor_graph, vertex.label = V(cor_graph)$name, vertex.size = 5, vertex.label.cex = 0.7, 
     main = "Red Metabólica Basada en Correlación Significativa")


# se aprecian posibles biomarcadores¡

###################### Validación de Biomarcadores con PLS-DA ##################

# me aseguro que los paquetes se han cargado

if (!requireNamespace("ropls", quietly = TRUE)) BiocManager::install("ropls")
library(ropls)

# hago el escalado: 

expresion_significativa <- scale(t(expresion_significativa))

# creo el modelo: 

plsda_model <- opls(expresion_significativa, grupo_info$Groups, predI = 2)
summary(plsda_model)

#visualizo los resultados:

#resumen del modelo
plot(plsda_model, typeVc = "Resumen del Modelo")
#visualización de las cargas
plot(plsda_model, typeVc = "Cargas")


###############################################################################
###############################################################################
##################### DESCARGA DEL CONTENEDOR ##################################

# Simplemente guardo el archivo con extensión .Rda y lo pongo en la ruta donde 
# tengo mi repositorio local vinculado a GitHub

#save(se, file = "C:/Users/Usuario/Desktop/Puello-Mora-Catherine-PEC1/
     #SE_CatherinePM_PEC1.Rda")

################################################################################
################################################################################
################## DATOS EN FORMATO TEXTO #####################################

# guardo mi matriz de expresión en un formato de texto, en este caso he elejido csv

#expresion_data <- as.data.frame(assay(se))
#write.csv(expresion_data, file = "C:/Users/Usuario/Desktop/Puello-Mora-Catherine-PEC1/Se_DatosE_en texto.csv", row.names = TRUE)

