library(tidyverse)
library(readxl)

path = "~/Documents/Bioinformatics Projects/Evidence for YhfX/2 Expression/"

#RNA-seq microarray data set from GEO, series GSE20305 (Jozefczuk et al. 2010)
array_raw <- as.data.frame(read_tsv(paste(path, "GSE20305_series_matrix.txt", sep = ""), skip = 60, col_names = T))

#Some quality adjustments made in Excel, obtaining sample names to replace "GSE____" colnames
samples <- read_xlsx(paste(path,"samples.xlsx", sep=""), col_names = F)
sampleID <- unlist(samples[1])
sampleName <- unlist(samples[2])

#gene names obtained from Series Matrix File on GEO Accession Display
geneNames <- read_xlsx(paste(path, "GEO_Array_geneNames.xlsx", sep = "")) %>%
  select(ID, ORF)

#aggregated file with all isomerases, decarboxylases, and racemases of E. coli
#obtained from Uniprot
all_ases <- read_csv(paste(path, "all_ases.csv", sep = ""))
all_ases <- unlist(all_ases[["Gene"]])

#joining array to all_ases data to match gene name
array <- as_tibble(array_raw) %>%
  #rename array ID references to gene name
  left_join(geneNames, by = c("ID_REF" = "ID")) %>%
  mutate(ID_REF = ORF) %>%
  select(-ORF) %>%
  rename(GENE = ID_REF) %>%
  subset(GENE %in% c("yhfX", "yhfW", "php", "yhfT", "yhfS", all_ases))

#rename columns to sample names
colnames(array)[-1] <- sampleName

#f(x) to summarize values by timepoint and stress condition
calc_mean <- function(ntimepoints, stress = c("heatstress", "coldstress", "oxidativestress","lactose","control")) {
  avg <- NULL
    for(i in 1:ntimepoints){
        avg[i] <- list(rowMeans(array[str_detect(colnames(array), pattern = paste(stress, "_timepoint", i, sep = ""))]))
        avg
    }
  avg
}


#application of above averaging f(x) and renaming list items
oxidative <- calc_mean(15, "oxidativestress")
names(oxidative) <- paste("oxidative", 1:length(oxidative))

lactose <- calc_mean(15, "lactose")
names(lactose) <- paste("lactose", 1:length(lactose))

control <- calc_mean(15, "control")
names(control) <- paste("control", 1:length(control))

cold <- calc_mean(15, "coldstress")
names(cold) <- paste("cold", 1:length(cold))

heat <- calc_mean(15, "heatstress")
names(heat) <- paste("heat", 1:length(heat))

not_any_na <- function(x) all(!is.na(x))


#concat values, remove NA cols
array_average <- data.frame(array[1], control, cold, heat, oxidative, lactose) %>%
                    select(where(not_any_na))


#file to use in HemI 2.0 for heatmap visualization
write_excel_csv(array_average, file = paste(path, "HemI_allases.csv", sep = ""), col_names = T)

#range check for consider transformation in HemI
max(array_average[-1]) #max = 205540.7
min(array_average[-1]) #min = 0.462717

#wide range, log10 transformation will be applied
