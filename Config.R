# config.R
# Глобальные параметры для анализа

# ПУТИ К ФАЙЛАМ

DATA_DIR <- "data_all"

PHENOTYPE_FILE <- file.path(DATA_DIR, "phenotable.tsv")
COUNTS_FILE <- file.path(DATA_DIR, "miR.Counts.csv")

OUTPUT_DIR <- "results"
FIGURES_DIR <- "figures"

# Создание директорий если их нет
dir.create(OUTPUT_DIR, showWarnings = FALSE)
dir.create(FIGURES_DIR, showWarnings = FALSE)


# ПАРАМЕТРЫ ФИЛЬТРАЦИИ
# Минимальный уровень экспрессии
MIN_COUNTS <- 10

# Минимальное количество образцов
MIN_SAMPLES <- 3

# Порог для значимости
PADJ_THRESHOLD <- 0.05
LOG2FC_THRESHOLD <- 1.0

# Порог для "высокой" экспрессии (для Venn диаграммы)
HIGH_EXPRESSION_THRESHOLD <- 150

# ПАРАМЕТРЫ НОРМАЛИЗАЦИИ

# Метод трансформации: "rlog" или "vst"
TRANSFORM_METHOD <- "rlog"

# ПАРАМЕТРЫ ВИЗУАЛИЗАЦИИ

# Количество топ генов для heatmap
TOP_GENES_HEATMAP <- 50

# Размер текста на графиках
BASE_TEXT_SIZE <- 12

# Цветовая палитра для групп
GROUP_COLORS <- c(
  "no_complications" = "#4DAF4A",
  "humoral" = "#E41A1C",
  "cellular" = "#377EB8",
  "TCAD" = "#FF7F00"
)

# ПАРАМЕТРЫ ENRICHMENT АНАЛИЗА

# P-value cutoff для обогащения
ENRICH_PVAL_CUTOFF <- 0.05

# Минимальный размер генного набора
MIN_GENESET_SIZE <- 10

# Максимальный размер генного набора
MAX_GENESET_SIZE <- 500


CONTRASTS <- list(
  humoral = c("condition", "AMR", "NR"),
  cellular = c("condition", "ACR", "NR"),
  TCAD = c("condition", "CAV", "NR")
)

set.seed(42)

# НЕОБХОДИМЫЕ БИБЛИОТЕКИ

required_packages <- c(
  "DESeq2", "ggplot2", "pheatmap", "RColorBrewer", 
  "vsn", "hexbin", "tidyverse", "ggrepel", 
  "EnhancedVolcano", "VennDiagram", "clusterProfiler",
  "org.Hs.eg.db", "multiMiR", "miRBaseConverter", "ggVennDiagram"
)

# Функция для установки недостающих пакетов
install_if_missing <- function(packages) {
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      if (pkg %in% c("DESeq2", "clusterProfiler", "org.Hs.eg.db", 
                     "EnhancedVolcano", "vsn", "hexbin")) {
        BiocManager::install(pkg)
      } else {
        install.packages(pkg)
      }
    }
  }
}

# Загрузка всех библиотек
suppressPackageStartupMessages({
  lapply(required_packages, library, character.only = TRUE)
})

cat("✓ Configuration loaded successfully\n")
