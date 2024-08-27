suppressPackageStartupMessages(library(gplots))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(tidyverse))
options(scipen=999)
options(warn=-1)

#Receiving the List with IDs and the output prefix
sample_names_all <- c("sample_names_braziliensis", "Lourador_Lbraziliensis")

sample_names_all <- commandArgs(trailingOnly=TRUE)

#Checking if the exact number of files were submitted to the script:
if (length(sample_names_all)!=2) {
  stop("You need exactly 2 arguments: The name of your file with sample IDs and the output prefix, as in Rscript CI_Estimation_server.v4.R sample_names_braziliensis Lourador_Lbraziliensis", call.=FALSE)
} 

#Receiving the list of samples and output prefix:
sample_names <- read.table(sample_names_all[1])
complexity_table <- data.frame()
temp_name <- sample_names_all[2]

#Setting the image sizes based on the number os samples:
number_of_samples <- length(sample_names$V1)

print("Reading and plotting the complexity results!")

#Loop to read each "table.simulation_replicates_final.csv" file and combine in a single table
i=1
for (i in 1:nrow(sample_names)) {
  sample_id_tem <- sample_names[i,"V1"]
  temp_input_table <- read.table(paste(sample_id_tem, "table.simulation_replicates_final.csv", sep = "."), sep = ";", header = TRUE)
  complexity_table <- rbind(complexity_table, temp_input_table)
}

complexity_table_all <- complexity_table

#Filtering by coverage and Het_SNP_Count:
complexity_table <- complexity_table[complexity_table$ct>=100,]
complexity_table <- complexity_table[complexity_table$coverage>=25,]

#Setting the cutoff based on population data:
complexity_table_sim <- complexity_table[complexity_table$Pop=="Simulated",]
complexity_table_sim_mean <- mean(complexity_table_sim$md)
complexity_table_sim_SD <- sd(complexity_table_sim$md)
cutoff <- complexity_table_sim_mean+3*complexity_table_sim_SD

#Identifying the complex samples
complexity_table$species <- temp_name
complexity_table$Complex <- "No"
complexity_table$Complex[complexity_table$md >= 0.1 & complexity_table$Q_value_BH < 10^-10 & complexity_table$Prop_Chr_CM_higher_0.1 >= 50] <- "Complex"
complexity_table$Complex[complexity_table$md <= 0.1 & complexity_table$md >= cutoff & complexity_table$Q_value_BH < 10^-10 & complexity_table$Prop_Chr_CM_higher_0.1 >= 50] <- "Pot_Complex"

#Ploting the Complexity Mean and Complexity SD
png(paste (temp_name,"complexity_all_merged.png", sep = "."), width = 1500, height = 1000, units = "px", pointsize = 30, res = 200)
print(ggplot(complexity_table, aes(x=md,y=sd, color=Complex, shape= Pop)) + geom_point(alpha=0.5, size = 3) + theme_bw() + 
  scale_color_manual(values = c( "gray", "red", "orange") , breaks = c("No","Complex", "Pot_Complex")) + scale_shape_manual(values = c(16, 17, 18)) +
  geom_vline(xintercept = 0.1, color = "red", linetype='dashed', size = 0.5) +
  geom_vline(xintercept = cutoff, color = "orange", linetype='dashed', size = 0.5) +
  theme(legend.key.size = unit(1, 'cm')) + xlab("Complexity Mean") + ylab("Complexity SD"))
invisible(dev.off())

#Ploting the Complexity Mean and Proportion of complex chromosomes
png(paste(temp_name,"complexity_all_merged_chromosome.png", sep = "."), width = 1500, height = 1000, units = "px", pointsize = 30, res = 200)
print(ggplot(complexity_table, aes(x=md,y=Prop_Chr_CM_higher_0.1, color=species, shape= Complex)) + geom_point(alpha=0.5, size = 3) + theme_bw() + 
  scale_color_manual(values = c( "gray", "red") , breaks = c("No","Complex")) + scale_shape_manual(values = c(16, 17, 18)) +
  geom_vline(xintercept = 0.1, color = "red", linetype='dashed', size = 0.5) +
  theme(legend.key.size = unit(1, 'cm')) + xlab("Complexity Mean") + ylab("Proportion of complex chromosomes"))
invisible(dev.off())

#Ploting the Complexity and Proportion of complex chromosomes with names and cutoffs
png(paste(temp_name,"complexity_all_separated_chr_new.png",sep = "."), width = 2500, height = 2500, units = "px", pointsize = 30, res = 300)
print(ggplot(complexity_table, aes(x=md,y=Prop_Chr_CM_higher_0.1, color=Chr_evaluated, shape= Complex)) + geom_point(alpha=0.5, size = 3) + theme_bw() + 
  #scale_color_manual(values = c("#DC143C", "dark golden rod", "cyan", "green", "darkgreen", "darkblue" )) +
   scale_shape_manual(values = c(16, 17, 18, 15)) +
  geom_hline(yintercept = 50, color = "red", linetype='dotted', size = 1) + ylab("Prop.Complex.Chr") + xlab("Complexity") +
  # geom_vline(xintercept = EA_frassen_cutof, color = "#DC143C", linetype='dotted', size = 1) +
  # geom_vline(xintercept = EA_HIV_cutoff, color = "orange", linetype='dotted', size = 1) +
  geom_vline(xintercept = 0.1, color = "red", linetype='dotted', size = 1) +
  geom_text_repel(aes(label=samp_name),  size=2, max.overlaps=Inf,  show.legend=FALSE, segment.linetype="dotted") + theme(legend.key.size = unit(1, 'cm')))
invisible(dev.off())

#Same as before, but more compact
png(paste(temp_name,"complexity_all_separated_chr_new_short.png",sep = "."), width = 1500, height = 1500, units = "px", pointsize = 30, res = 300)
print(ggplot(complexity_table, aes(x=md,y=Prop_Chr_CM_higher_0.1, color=Chr_evaluated, shape= Complex)) + geom_point(alpha=0.5, size = 3) + theme_bw() + 
  #scale_color_manual(values = c("#DC143C", "dark golden rod", "cyan", "green", "darkgreen", "darkblue" )) +
  scale_shape_manual(values = c(16, 17, 18, 15)) +
  geom_hline(yintercept = 50, color = "red", linetype='dotted', size = 1) + ylab("Prop.Complex.Chr") + xlab("Complexity") +
  # geom_vline(xintercept = EA_frassen_cutof, color = "#DC143C", linetype='dotted', size = 1) +
  # geom_vline(xintercept = EA_HIV_cutoff, color = "orange", linetype='dotted', size = 1) +
  geom_vline(xintercept = 0.1, color = "red", linetype='dotted', size = 1) +
  geom_vline(xintercept = cutoff, color = "orange", linetype='dashed', size = 0.5) +
  geom_text_repel(aes(label=samp_name),  size=2, max.overlaps=Inf,  show.legend=FALSE, segment.linetype="dotted") + theme(legend.key.size = unit(1, 'cm')))
invisible(dev.off())

#Writting the complexity table with all samples
write.table(complexity_table, paste(temp_name, "Table_with_updated_complex.csv",sep = "."), sep = ",", col.names = TRUE, row.names = FALSE, quote = FALSE )

############################################################
#Generating the AARD density plots for all sets:

#Removing samples without 100 SNPs an < 25 coverage
sample_names <- sample_names[sample_names$V1 %in% complexity_table$samp_name, , drop=FALSE]
print("Reading and plotting the AARD distributions!")

#Reading all the AARD files: "real_data_to_MH.csv" and "simulated_data_to_MH.csv"
a=1
All_samples_df <- data.frame()
SNP_counts_dataframe <- data.frame()
for (a in 1:nrow(sample_names)) {
  temp_name2 <- sample_names[a,1]
  
  isc2_real <- read.table(paste(temp_name2, "real_data_to_MH.csv", sep = "."), sep = ";", header = TRUE)
  isc2_sim <- read.table(paste(temp_name2,"simulated_data_to_MH.csv", sep = "."), sep = ";", header = TRUE)
  isc2_real2 <- isc2_real %>% select(Sample, coverage, median_coverage, ALTF, RDR, RDA, Complexity, MAF)
  isc2_real2$type = "Real"
  isc2_sim$type = "Sim"
  isc2_sim$Sample <- isc2_real2$Sample
  isc2_all_data <- rbind(isc2_real2 , isc2_sim)
  
  #Getting the SNP counts
  snp_count <- nrow(isc2_real)
  temp_snp_count <- data.frame(Sample = temp_name2, SNPs = snp_count)
  SNP_counts_dataframe <- rbind(SNP_counts_dataframe, temp_snp_count)
  
  All_samples_df <- rbind(All_samples_df, isc2_all_data)
}


#Setting the complex isolates to a different colour:
All_samples_df2 <- All_samples_df
All_samples_df$Complex <- "No"

complex_isolates <- complexity_table$samp_name[complexity_table$Complex=="Complex"]
Pot_Complex_isolates <- complexity_table$samp_name[complexity_table$Complex=="Pot_Complex"]

All_samples_df$Complex[All_samples_df$type=="Real"] <- "Not_complex"
All_samples_df$Complex[All_samples_df$type=="Sim"] <- "Simulated_clone"
All_samples_df$Complex[All_samples_df$Sample %in% complex_isolates & All_samples_df$type =="Real" ] <- "Complex"
All_samples_df$Complex[All_samples_df$Sample %in% Pot_Complex_isolates & All_samples_df$type =="Real" ] <- "Pot_Complex"

#Plotting the AARD density graphs
png(paste(temp_name, "simulated_data_distribution_complex", "png", sep = "."), width = (number_of_samples*150)+1000, height =(number_of_samples*150)+500, units = "px", pointsize = 30, res = 300)
print(ggplot(All_samples_df, aes(y=ALTF, fill = Complex)) + geom_density(alpha = 0.4) + theme_bw() + theme(axis.text.x = element_text(angle = 90), strip.text.x = element_text(size = 10)) +
        ylab("AARD") +  facet_wrap(~Sample, scales = "free_y") + coord_flip() + geom_hline(yintercept=0.5, color = "red") + geom_hline(yintercept=0.25, color = "blue") +
        geom_hline(yintercept=0.75, color = "blue") + 
        scale_fill_manual(values = c( "#48D1CC","white", "#DC143C", "orange"), breaks = c("Not_complex","Simulated_clone","Complex", "Pot_Complex"))) + ylim(0,1)
invisible(dev.off())

############################################################

#Ploting the SNP counts:
print("Plotting the SNPcount data!")

complexity_table_all_nocont <- complexity_table_all[!grepl("Control.", complexity_table_all$samp_name),]
complexity_table_all_nocont$Class <- "Not usable"
complexity_table_all_nocont$Class[complexity_table_all_nocont$coverage >=20 & complexity_table_all_nocont$ct >= 100] <- "OK"

png(paste(temp_name, "usable_samples", "png", sep = "."), width = 1600, height =1300, units = "px", pointsize = 30, res = 300)
ggplot(complexity_table_all_nocont, aes(x=ct,y=coverage, label=samp_name, color=Class)) + geom_point() +
  geom_hline(yintercept = 25, color = "red", linetype='dotted', size = 1) + ylab("Genome Coverage") + xlab("Het SNP count") +
  # geom_vline(xintercept = EA_frassen_cutof, color = "#DC143C", linetype='dotted', size = 1) +
  # geom_vline(xintercept = EA_HIV_cutoff, color = "orange", linetype='dotted', size = 1) +
  geom_vline(xintercept = 100, color = "red", linetype='dotted', size = 1) + theme_bw() + scale_x_log10() + scale_y_log10() +
  geom_text_repel(aes(label=samp_name),  size=2, max.overlaps=Inf,  show.legend=FALSE, segment.linetype="dotted") + theme(legend.key.size = unit(1, 'cm')) +
  scale_color_manual(values = c("#DC143C", "#1E90FF"), breaks = c("Not usable","OK") )
invisible(dev.off())

print("All done!")

