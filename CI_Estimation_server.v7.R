#CI Estimation V2, model Server
#Author: Joao Luis Reis Cunha, University of York, Jeffares group
#This script receives a VCF file, genome coverage and somy estimations of chromosomes and estimates the isolate complexity

#Note: Chromosome names have to be only numerals, as in 1, 2, 3, etc..
#We strongly suggest that the VCF should be filtered to retain only single-copy genes prior to be used here

#How to run:

#For a single sample:
#Rscript CI_Estimation_server.v4.R ERR205724 ERR205724.rec2.recode.vcf
#where ERR205724 is your sample name and ERR205724.rec2.recode.vcf is your vcf file

#For a list of samples in the file "sample_names":
#for i in $(cat sample_names); do Rscript CI_Estimation_server.v4.R ${i} ${i}.rec2.recode.vcf ; done
#Where "sample_names" is a test file containing the IDs of your samples, one at each line, as in:
#ERR205724
#ERR205726
#ERR205730
#ERR205731
#ERR205732
#ERR205741

#The script automatically removes chromosomes and genes with extra copies 
#To that end the script estimate the genome coverage, and chromosome copies directly from the VCF file
#This is NOT recommended. It is best to provide files containing the genome coverage of each isolate as well as its chromosomal copies

#The genome coverage file should be a text file named "Genome_Cov.Table", have no header, and two columns. A TAB separated file.
#The first column is the sequence ID and the second the genome coverage, as in:

# ERR3956078	25
# ERR3956080	26.00
# ERR3956084	28
# ERR3956086	25
# ERR3956093	25.00
# ERR3956094	26
# ERR3956095	27.00

#The Chromosome copies file has to be named "CCNV.Table.ordered", a TAB separatted file
#The first column correspond to the chromosome number, similar to what was present in the VCF file, and the header must be "Chromosome"
#Chromosome names have to be only numerals, as in 1, 2, 3...
#The headers from columns 2 to "n" must be the sample IDs, identical to the ones in the "sample_names" file.
#The lines correspond to each chromosome copy number
#And example of the "CCNV.Table.ordered" file for the isolates ERR3956070 and ERR3956071, from L. donovani, that has 36 chromosomes

#Chromosome	ERR3956070	ERR3956071
#1	1	0.9474
#2	1.048	1
#3	1.048	1
#4	1	1
#5	1	1
#6	1.048	1
#7	1.048	1
#8	1	0.9474
#9	1	1
#10	1	1
#11	1	1
#12	1.048	1
#13	1	0.9474
#14	1	1
#15	1	1
#16	1	1
#17	1	1
#18	1.048	1
#19	1	0.9474
#20	1	1
#21	1	1
#22	1	1
#23	1	0.9474
#24	1.048	0.9474
#25	1	0.9474
#26	1.048	1
#27	1	0.9474
#28	1	0.9474
#29	1	0.9474
#30	1	0.9474
#31	1.952	1.895
#32	1	0.9474
#33	1	1
#34	1	0.9474
#35	1	0.9474
#36	1	0.9474

#The files "sample_names", "Genome_Cov.Table" and "CCNV.Table.ordered" have to be in the same folder that you are running the script

#The SNP and coverage cutoffs can be seen and altered in the section:
#Selecting the cutoffs

print("Loading the required libraries")
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(viridis))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(vcfR))
suppressPackageStartupMessages(library(ggrepel))

###################################################################################################################################
#Loading input files and Data:

#This is the version of the script to run with batch samples
sample_names <- commandArgs(trailingOnly=TRUE)

#Checking if the exact number of files were submitted to the script:
if (length(sample_names)!=2) {
  stop("You need exactly 2 arguments: Your sample ID and the vcf file name, as in Rscript CI_Estimation_server.v4.R ERR205724 ERR205724.rec2.recode.vcf", call.=FALSE)
} 

#Saving the sample name and the vcf file
sample_name <- sample_names[1]
vcf_file <- sample_names[2]

print(paste("Processing sample ",sample_name, "and the vcf file", vcf_file, sep = " "))

#mporting VCF
#Change the suffix "rec2.recode.vcf" of this to the 
vcf2 <- read.vcfR(vcf_file, verbose = FALSE )# VCF file
#Note: Chr names have to be only numerals, as in 1, 2, 3, etc..
#We strongly suggest that the VCF should be filtered to retain only single-copy genes prior to be used here

#You can provide CCNV predictions, by having a TSV file named "CCNV.Table.ordered" in the folder, with the chromosome copy number values from your sample, 
#You can also provide a TSV with the genome coverage, by having a file named "Genome_Cov.Table" in the folder

#If they are not provided, these values will be estimated directly from the VCF file (not recommended)

if (file.exists("CCNV.Table.ordered")) {
input_ccnv <- read.csv("CCNV.Table.ordered", sep = "\t", header = TRUE, check.names = FALSE)
input_ccnv <- input_ccnv[,c("Chromosome",sample_name)] #CSV with the CCNV data
print("CCNV file found!!")
}

if (file.exists("Genome_Cov.Table")) {
genome_cov_temp <- read.csv("Genome_Cov.Table", sep = "\t", header = FALSE, check.names = FALSE)
input_cov <- genome_cov_temp[genome_cov_temp$V1==sample_name, "V2"] #Mean genome coverage value
print("Genome coverage file found!!")
}


###################################################################################################################################
#Selecting the cutoffs:
#OBS: These are the recommended values. Usually you would not need to change those.

#Read depth in SNPs
cutoff_read_depth_alternate = 5 # Minimum number of reads in the alternate allele for the SNP position to be kept
dubious_SNP_minor_cutoff = 1 #  Minimum number of reads in the reference allele for the SNP to be classified as "dubious"
dubious_SNP_major_cutoff = 5 #  Maximum number of reads in the reference allele for the SNP to be classified as "dubious"
SNP_minimun_read_depth = 15 #  Minimum number of reads in the position for the SNP to be kept

#Genome coverage:
minimun_coverage = 25 # Minimum genome coverage for the sample to be kept

#CCNV cutoff
#Somy data values - Important to remove chromosomes with evidence of extra copies or loss of copies
CCNV_low = 0.85  # Minimum Chromosome copy number
CCNV_high = 1.15 # Maximum Chromosome copy number

#SNP coverage cutoff
SNP_cov_low = 0.85 #Minimum normalized coverage when compared to the chromosome coverage in which the SNP is
SNP_cov_high = 1.15 #Maximum normalized coverage when compared to the chromosome coverage in which the SNP is


###################################################################################################################################

#Processing the analysis
print("Importing and processing the VCF data!!")

############# SNP data ############# 
#Reading the VCF file and transforming it into a dataframe
sample_vcf <- vcfR2tidy(
  vcf2,
  info_only = FALSE,
  single_frame = FALSE,
  toss_INFO_column = TRUE
)

#Obtaining the genotype and read counts in each allele
sample_vcf_gt <- as.data.frame(sample_vcf$gt)

#Correcting the chromosome IDs and separating the columns with genotypes and read depths of alleles
sample_vcf_chr_names <-  unique(as.data.frame(sample_vcf$fix[,1:2]))
sample_vcf_gt2 <- merge(sample_vcf_gt, sample_vcf_chr_names, by = "ChromKey")

#Selecting the sample name:
sample_name = sample_vcf_gt2[1,colnames(sample_vcf_gt2)=="Indiv"]

#Parte para os dados da Laila:
# sample_vcf_gt2$CHROM <- gsub("-S","", sample_vcf_gt2$CHROM)
# sample_vcf_gt2$CHROM <- gsub("TcChr","", sample_vcf_gt2$CHROM)

sample_vcf_gt3 <- sample_vcf_gt2[!is.na(sample_vcf_gt2$gt_GT),]
sample_vcf_gt4 <- sample_vcf_gt3 %>% separate(gt_AD, c("RDR", "RDA"), ",")
sample_vcf_gt5 <- sample_vcf_gt4 %>% separate(gt_GT_alleles, c("All1", "All2"), '/')
sample_vcf_gt5 <- sample_vcf_gt5 %>% separate(gt_GT, c("GEN1", "GEN2"), '/')

#Setting the Unique SNP ID for the table:
sample_vcf_gt5$SNP_Uniq <- paste(sample_vcf_gt5$CHROM, sample_vcf_gt5$POS, sep = "_")

#Setting the genotype read depth as numeric
sample_vcf_gt5$RDR <- as.numeric(as.character(sample_vcf_gt5$RDR))
sample_vcf_gt5$RDA <- as.numeric(as.character(sample_vcf_gt5$RDA))
sample_vcf_gt5$gt_DP <- as.numeric(as.character(sample_vcf_gt5$gt_DP))

#Cleaning temporary data
rm(sample_vcf_gt, sample_vcf_gt2, sample_vcf_gt3, sample_vcf_gt4, vcf2)

#############################################################################################
#Obtaining the genome coverage and chromosome copy number variation directly from VCF files##
#This will only be done if the coverage and CCNV were not provided in the initial parts    ##
#############################################################################################

#Estimating potential coverage and somy:
#Estimating coverage
#Finding chr with extra copies:
#Getting the raw mean and median coverage
#I will use the median coverage
temp_coverage_mean <- mean(sample_vcf_gt5$gt_DP)
temp_coverage_median <- median(sample_vcf_gt5$gt_DP)

temp_chr_cov <- data.frame(sample_vcf_gt5 %>% group_by(CHROM) %>% summarise(mean_coverage = mean(gt_DP/temp_coverage_mean)))
#Selecting the most common somy as a second normalizer:
#This is the median value of the chromosome somies for all chromosomes. This corrects the bias of inflating counts from trisomic/tetrasomic chromosomes
chr_normalizer2 <- median(temp_chr_cov$mean_coverage)
temp_chr_cov$mean_coverage <- temp_chr_cov$mean_coverage/chr_normalizer2

#Selecting the euploid chromosomes to be in the analysis
chr_between_0.85_1.15 <- temp_chr_cov$CHROM[temp_chr_cov$mean_coverage > 0.85 & temp_chr_cov$mean_coverage <=1.15]

#Getting the mean genome coverage only using chrs in the interval of global coverage of 0.85-1.15
mean_cov_temp <- as.numeric(as.character(sample_vcf_gt5 %>% filter(CHROM %in% chr_between_0.85_1.15) %>% #group_by(CHROM) %>%
  summarise(mean_coverage = mean(gt_DP))))

#Getting the normalized chromosome coverage for each chromosome
chr_coverage_temp <- data.frame(sample_vcf_gt5 %>%  group_by(CHROM) %>%
  summarise(mean_coverage = mean(gt_DP/mean_cov_temp)))

colnames(chr_coverage_temp) <- c("Chromosome", sample_name)

#Checking if the coverage and CCNV files already exist or need to be generated
if (!exists("input_ccnv")) { 
input_ccnv <- chr_coverage_temp
print("Estimating manually the Chromosome copy number!!")

}
if (!exists("input_cov")) {
input_cov <- mean_cov_temp
print("Estimating manually the Genome Coverage!!")

}

print("Processing the Data!!")

############# CCNV data ############# 
#Processing the CCNV data:
input_ccnv_melt <- melt(input_ccnv, id=c("Chromosome"))
colnames(input_ccnv_melt) <- c("CHROM", "Indiv", "CCN")
input_ccnv_melt$CHROM <- as.character(input_ccnv_melt$CHROM)

#Merging the CCNV data with the VCF data
input_ccnv_melt$CHROM <- as.character(input_ccnv_melt$CHROM)
sample_vcf_gt5_temp <- sample_vcf_gt5 %>% left_join(input_ccnv_melt, by=c("CHROM","Indiv"))
#Removing chromosomes with extra or loss of copies
sample_vcf_gt5_temp <- sample_vcf_gt5_temp[sample_vcf_gt5_temp$CCN > CCNV_low & sample_vcf_gt5_temp$CCN < CCNV_high,]

############# Coverage data ############# 
#Getting the coverage data:
sample_vcf_gt5_temp2 <- sample_vcf_gt5_temp
sample_vcf_gt5_temp2$median_coverage <- input_cov
sample_vcf_gt5_temp2$depth <- as.numeric(as.character(sample_vcf_gt5_temp2$RDR)) + as.numeric(as.character(sample_vcf_gt5_temp2$RDA))

############# Filtering SNPs by SNP coverage normalized by chromosome coverage ############# 
#Filtering by coverage taking into account CCN:
sample_vcf_gt5_temp2$cutff_up <- SNP_cov_low*sample_vcf_gt5_temp2$CCN*sample_vcf_gt5_temp2$median_coverage
sample_vcf_gt5_temp2$cutff_down <- SNP_cov_high*sample_vcf_gt5_temp2$CCN*sample_vcf_gt5_temp2$median_coverage
sample_vcf_gt5_temp3 <- sample_vcf_gt5_temp2[sample_vcf_gt5_temp2$depth > sample_vcf_gt5_temp2$cutff_up & sample_vcf_gt5_temp2$depth < sample_vcf_gt5_temp2$cutff_down,] 

#Adequating the data Columns to the script:
vcf_gt6 <- sample_vcf_gt5_temp3 %>% select(CHROM, POS, All1, All2, RDR, RDA, Indiv, gt_DP, CCN, median_coverage, cutff_up, cutff_down)
colnames(vcf_gt6) <- c("Chr", "Pos", "Ref", "Alt", "RDR", "RDA", "Sample", "Coverage", "CCN", "median_coverage", "cutff_up", "cutff_down")

######Filtering SNPs ############# 
#Getting the MAF, REFF and Total Coverage:
vcf_gt6$RDR <- as.numeric(as.character(vcf_gt6$RDR))
vcf_gt6$RDA <- as.numeric(as.character(vcf_gt6$RDA))

#Estimating MAF and ALTF
vcf_gt6$MAF = pmin(vcf_gt6$RDR, vcf_gt6$RDA)/(vcf_gt6$RDR + vcf_gt6$RDA)
vcf_gt6$REFF = (vcf_gt6$RDR/(vcf_gt6$RDR + vcf_gt6$RDA))
vcf_gt6$coverage = vcf_gt6$RDR + vcf_gt6$RDA
vcf_gt6$ALTF = (vcf_gt6$RDA/(vcf_gt6$RDR + vcf_gt6$RDA))

#Removing positions with coverage < 15
vcf_gt6 <- vcf_gt6[vcf_gt6$coverage >=SNP_minimun_read_depth,]

#Removing samples with mean coverage below 15
vcf_gt6 <- vcf_gt6[vcf_gt6$median_coverage >=minimun_coverage,]

#Classifying the SNPs in Homo, Hetero or Dubious
vcf_gt6 <- vcf_gt6[vcf_gt6$RDA >=cutoff_read_depth_alternate,] #At least 5 reads in the alternate allele
vcf_gt6$Snp_Type <- "hom"
vcf_gt6$Snp_Type[vcf_gt6$RDA >cutoff_read_depth_alternate & vcf_gt6$RDR <dubious_SNP_major_cutoff & vcf_gt6$RDR >=dubious_SNP_minor_cutoff] <- "Dub"
vcf_gt6$Snp_Type[vcf_gt6$RDR >cutoff_read_depth_alternate & vcf_gt6$RDA <dubious_SNP_major_cutoff & vcf_gt6$RDA >=dubious_SNP_minor_cutoff] <- "Dub"
vcf_gt6$Snp_Type[vcf_gt6$RDA >=cutoff_read_depth_alternate & vcf_gt6$RDR >=dubious_SNP_major_cutoff] <- "het"
count_table_SNPs <- vcf_gt6 %>% count(Sample, Snp_Type)
colnames(count_table_SNPs) <- c("Sample","Snp_Type","Filt_mq30")

#Getting the snp counts
Snp_count_overall <- vcf_gt6 %>% count(Sample, Snp_Type)
Snp_count_chr <- vcf_gt6 %>% count(Sample, Chr, Snp_Type)
write.table(Snp_count_overall, file =paste(sample_name,"SNP_counts_table.csv", sep = "."), quote = FALSE, sep = ";", row.names = FALSE)

#Estimating the normalized coverage:
vcf_gt6$cov_norm <- vcf_gt6$Coverage/vcf_gt6$median_coverage

#Plot with the SNP counts that was classified in Homozygous, Heterozygous or Dubious
png(paste(sample_name, "Plot_SNPs_overall_pre.png", sep = "."), width = 1000, height = 1000, units = "px", pointsize = 30, res = 200)
print(ggplot(Snp_count_overall, aes(y=n, x=Sample,  fill=Snp_Type)) + geom_col(position = "dodge") + theme_bw() +
  ylab("SNP_count") + xlab("Isolate") + scale_fill_manual(values = c("#C0C0C0","#40E0D0", "#008B8B")))
invisible(dev.off())

vcf_gt6 <- vcf_gt6[!is.na(vcf_gt6$Sample),]

#Recoveing only the Heterozygous SNPs
input_table_het3 <- vcf_gt6[vcf_gt6$Snp_Type == "het",]
input_table_het3$Chr <- factor(input_table_het3$Chr, levels = unique(input_table_het3$Chr))

#Estimating the normalized coverage for each SNP position
input_table_het3$cov_norm <- input_table_het3$Coverage/input_table_het3$median_coverage

#Plots tge AARD in eacg SNP position in each non super or under-numerary chromosome
png(paste(sample_name, "Coverage_across_chr_SNPs_scaled_cov_SMAF_inv.png", sep = "."), width = 6000, height = 3000, units = "px", pointsize = 30, res = 200)
print(ggplot(input_table_het3, aes(x = Pos, y=ALTF, color = cov_norm, shape = Snp_Type)) + geom_point() + theme_bw() +
 theme(legend.key.size = unit(2, "cm"),legend.text = element_text(size = 15),
       legend.title = element_text(size = 15), legend.key.height= unit(2, 'cm'), legend.key.width= unit(2, 'cm'),
       axis.title.x = element_text(size = 30),axis.title.y = element_text(size = 30)) + labs(y = "AARD", x = "Chr pos") +
 facet_wrap(~Chr, scales="free_x") + ylim(0,1) + scale_x_continuous(labels = scales::comma) +
 scale_colour_viridis_c(direction=-1) + geom_hline(yintercept=0.5, color = "red") + geom_hline(yintercept=0.25, color = "blue") +
 geom_hline(yintercept=0.75, color = "blue"))
invisible(dev.off())

#Plot the Cov in home and het snps
png(paste(sample_name,"Coverage_across_Snps_het_hom_filt.png", sep = "."), width = 600, height = 600, units = "px", pointsize = 30, res = 200)
  print(ggplot(vcf_gt6, aes(x = Snp_Type, y=cov_norm, fill = Snp_Type)) + geom_boxplot() + scale_fill_manual(values = c("#C0C0C0","#40E0D0","#008B8B")) +
  theme_bw() + theme(axis.text.x = element_text(angle = 90)))
invisible(dev.off())

# plotting the distribution based on the REFF
png(paste(sample_name,"Ref.allele.freq.Density.png", sep = "."), width = 1000, height = 500, units = "px", pointsize = 30, res = 200)
  print(ggplot(input_table_het3, aes(y=REFF, fill = "")) + geom_density() + theme_bw() + theme(axis.text.x = element_text(angle = 90), strip.text.x = element_text(size = 20)) +
  ylab("AARD") +  facet_wrap(~Sample, ncol = 7, scale = "free_y") + coord_flip() + geom_hline(yintercept=0.5, color = "red") + geom_hline(yintercept=0.25, color = "blue") +
  geom_hline(yintercept=0.75, color = "blue") + scale_fill_manual(values = "#9400D3"))
invisible(dev.off())

#Getting the Complexity Index (CI), the absolute value of the Alternate allele freq subtracted by 0.5. This represents the variation from the "0.5" expected value
input_table_het3$Complexity <- abs(input_table_het3$ALTF - 0.5)
input_table_het3$Sample <- as.character(input_table_het3$Sample)

#Plotting the complexity frequency in each chromosome
png(paste(sample_name,"Compl.chr.png", sep = "."), width = 2000, height = 1000, units = "px", pointsize = 30, res = 200)
  print(ggplot(input_table_het3, aes(x=Chr, y=Complexity)) + geom_boxplot() + theme_bw() + xlab("Chromosome") +
        theme(axis.text.x = element_text(angle = 90)) + ylab("CI")) +  geom_hline(yintercept=0.1, linetype="dashed", color = "red")
invisible(dev.off())
 
#Saving the heterozygous input table:
input_table_het3$Pop <- sample_name
write.table(input_table_het3, file = paste(sample_name,"filtered.table",  "csv", sep = "."), quote = FALSE, sep = ";", row.names = FALSE)

#Generating the metrics across chromosomes
chromosomes_mean_median_complexity <- data.frame(input_table_het3 %>% group_by(Chr) %>% summarise(chr_mean_CM = mean(Complexity), chr_median_CM = median(Complexity), snps_evaluated=n()))

number_of_evaluated_chromosomes <- nrow(chromosomes_mean_median_complexity)
number_of_evaluated_chromosomes_mean_higher_0.1 <- nrow(chromosomes_mean_median_complexity[chromosomes_mean_median_complexity$chr_mean_CM >=0.1,])

###########################################################################################
#Generating the "Simulated data" with the same coverage and SNP count as the real isolate##
#And combining the simulated data and the real data                                      ##
###########################################################################################

print("Generating the simulated data!")
table_input <- input_table_het3
table_input_summary <- as.data.frame(summarise(group_by(table_input, Sample), coverage = median(median_coverage), soma=sum(Complexity), sd=sd(Complexity), 
                                               md=median(Complexity), ct =n()))
#Getting the vector of 1 and 0s:
vector.of.ones.and.zeros<-rep(c(0,1),10000)

#Runnig the replicates sampling from the population:
snp_prop_dataframe <- data.frame()

#Selecting the data from a ginven isolate and getting the simulation:
name_sample_mean_cov <- unique(table_input$median_coverage)
    
c=1
for (c in 1:nrow(table_input)) {
    tempcov = table_input[c,"coverage"]
    allele.count<-sum(sample(vector.of.ones.and.zeros,tempcov))  
    ref_allele_count <- tempcov - allele.count 
    proportion <-allele.count/tempcov 
    
    #add to vector
    temp_df <- data.frame(  Sample = c(paste("Control", sample_name, sep = ".")), coverage = c(tempcov), median_coverage =c(name_sample_mean_cov),
                            ALTF  = c(proportion),
                            RDR = c(ref_allele_count), RDA = c(allele.count))
    temp_df$Complexity <- abs(0.5-temp_df$ALTF)
    temp_df$MAF = pmin(temp_df$RDR, temp_df$RDA)/(temp_df$RDR + temp_df$RDA)
    snp_prop_dataframe <- rbind(snp_prop_dataframe, temp_df)
       
}
table_simmulation_summary_mean <- as.data.frame(summarise(group_by(snp_prop_dataframe, Sample), coverage = median(median_coverage), soma=sum(Complexity), sd=sd(Complexity), 
                                                            md=median(Complexity), ct =n()))

table_input_summary$Pop <- "Real"
table_simmulation_summary_mean$Pop <- "Simulated"

#Exporting the table with the replicates
write.table(table_input, file = paste(sample_name, "real_data_to_MH", "csv", sep = "."), sep = ";", quote = FALSE, row.names = FALSE)
write.table(snp_prop_dataframe, file = paste(sample_name, "simulated_data_to_MH", "csv", sep = "."), sep = ";", quote = FALSE, row.names = FALSE)
  
#Plot simulated data distribution_and_rea_togheter
data_frame_to_plot1 <- snp_prop_dataframe
data_frame_to_plot1$Sample <- table_input$Sample
data_frame_to_plot1$Type <- "Simulated"
data_frame_to_plot2 <- table_input
data_frame_to_plot2$Type <- "real"

data_frame_to_plot3 <- bind_rows(data_frame_to_plot1, data_frame_to_plot2) #Used dyplr to add NA in the simulated data

png(paste(sample_name, "simulated_data_distribution", "png", sep = "."), width = 1000, height = 500, units = "px", pointsize = 30, res = 200)
print(ggplot(data_frame_to_plot3, aes(y=ALTF, fill = Type)) + geom_density(alpha = 0.4) + theme_bw() +
        theme(axis.text.x = element_text(angle = 90), strip.text.x = element_text(size = 20), legend.key.size = unit(0.5, "cm")) +
        ylab("AARD") + coord_flip() + geom_hline(yintercept=0.5, color = "red") + geom_hline(yintercept=0.25, color = "blue") +
        geom_hline(yintercept=0.75, color = "blue") + scale_fill_manual(values = c("#9400D3", "#48D1CC")) + ylim(0,1))
invisible(dev.off())

###################################################################################################################################
##############################
#Doing the mantelhaen tests ##
##############################
print("Running the CMHtest!")

estimate_pop_mantelhaen_comb <- function(input_table1, input_table2) {
  temp_dataframe <- data.frame()
  
  input_table_vec <- input_table1
  input_table_sim <- input_table2
  sample_names <- as.vector(unique(input_table_vec$Sample))
  sample_sims <- as.vector(unique(input_table_sim$Sample))
  
  a = 1
  for (a in a:length(sample_names)) {
    samp_name <- sample_names[a]
    sample_sim <- sample_sims[a] 
    
    input_table_vec2 <- input_table_vec[input_table_vec$Sample == samp_name,]
    input_table_vec2 <- input_table_vec2[,c("RDR","RDA")]
    input_table_vec2$min <- apply(input_table_vec2,1, min) 
    input_table_vec2$max <- apply(input_table_vec2,1, max) 
    input_table_vec2 <- input_table_vec2[,c("min", "max")]
    
    
    input_table_sim2 <- input_table_sim[input_table_sim$Sample == sample_sim,]
    input_table_sim2 <- input_table_sim2[,c("RDR","RDA")]
    input_table_sim2$smin <- apply(input_table_sim2,1, min) 
    input_table_sim2$smax <- apply(input_table_sim2,1, max) 
    input_table_sim2 <- input_table_sim2[,c("smin", "smax")]
    
    input_table_combine1 <- cbind(input_table_vec2, input_table_sim2)
    
     input_table_vec3 <- input_table_combine1[,c("min","smin","max","smax")]
    unlisted_df <- unlist(as.data.frame(t(input_table_vec3)), use.names=FALSE)
    
    num_rows <- nrow(input_table_vec3)
    
    test_array_sample <- array(unlisted_df, dim = c(2,2,num_rows), 
                               dimnames = list(Alleles = c("real","simulated"), counts = c("1","2"),
                                               SNP_position = paste(input_table_vec2$Chr,input_table_vec2$Pos, sep = "_")))
    
    mant_test <- mantelhaen.test(test_array_sample)
    test_results_df <- data.frame(samp_name, num_rows, mant_test$statistic, mant_test$parameter, mant_test$p.value, mant_test$estimate, mant_test$null.value, 
                                  mant_test$alternative, mant_test$method, mant_test$data.name)
    
    temp_dataframe <- rbind(temp_dataframe, test_results_df)
    
  }  
  temp_dataframe$Q_value_BH <- p.adjust(temp_dataframe$mant_test.p.value, method = "BH")    
  return(temp_dataframe)
  
  
}
realsample_MT <- estimate_pop_mantelhaen_comb(table_input, snp_prop_dataframe)

print("Estimaing final CI results!")

#Getting all the data togheter:
colnames(table_input_summary)[1] <- "samp_name"
colnames(table_simmulation_summary_mean)[1] <- "samp_name"

table_input_summary <- merge(table_input_summary, realsample_MT, by = "samp_name")

table_input_summary$Q_value_BH <- as.numeric(as.character(table_input_summary$Q_value_BH))
table_input_summary$logQv <-  -log10(table_input_summary$Q_value_BH)

#This replace the "Inf" value for "1000" to allow downstream analysis
table_input_summary$logQv[table_input_summary$logQv =="Inf"] <- "1000"
table_input_summary$logQv <- as.numeric(as.character(table_input_summary$logQv))

#Including the number and proportion of chrs with Complexity higher than 0.1 in the output files chr
table_input_summary$Chr_evaluated <- number_of_evaluated_chromosomes
table_input_summary$Chr_CM_higher_0.1 <- number_of_evaluated_chromosomes_mean_higher_0.1
table_input_summary$Prop_Chr_CM_higher_0.1 <- (number_of_evaluated_chromosomes_mean_higher_0.1/number_of_evaluated_chromosomes)*100

table_simmulation_summary_mean$Chr_evaluated <-NA
table_simmulation_summary_mean$Chr_CM_higher_0.1 <- NA
table_simmulation_summary_mean$Prop_Chr_CM_higher_0.1 <- NA

table_real_simulation <- bind_rows(table_input_summary, table_simmulation_summary_mean) #Used dyplr to add NA in the simulated data
  
###################################################################################################################################
#Exporting table:
#Final table with the complexity estimations for the real and the simulated data
write.table(table_real_simulation, file = paste(sample_name, "table", "simulation_replicates_final", "csv", sep = "."), sep = ";", quote = FALSE, row.names = FALSE)

#Plotting the result:
table_real_simulation2 <- table_real_simulation
table_real_simulation2[table_real_simulation2$Pop=="Simulated", "Prop_Chr_CM_higher_0.1"] <- 
  table_real_simulation2[table_real_simulation2$Pop=="Real", "Prop_Chr_CM_higher_0.1"]
table_real_simulation2[table_real_simulation2$Pop=="Simulated", "Chr_evaluated"] <- 
  table_real_simulation2[table_real_simulation2$Pop=="Real", "Chr_evaluated"]

table_real_simulation2$Complex <- "No"
table_real_simulation2$Complex[table_real_simulation2$md>=0.1 & table_real_simulation2$Prop_Chr_CM_higher_0.1 >=50 & table_real_simulation2$Q_value_BH < 0.0000000001] <- "Yes"

png(paste(sample_name, "CI_final_result", "png", sep = "."), width = 1200, height = 1000, units = "px", pointsize = 30, res = 200)
print(ggplot(table_real_simulation2, aes(x=md,y=Prop_Chr_CM_higher_0.1,  shape= Complex, color=Complex,size=Chr_evaluated)) + geom_point(alpha=0.5) + theme_bw() +
   geom_vline(xintercept = 0.1, color = "red", linetype='dotted', linewidth = 1) +  geom_hline(yintercept = 50, color = "red", linetype='dotted', linewidth = 1) +
  ylab("Prop.Complex.Chr") + xlab("CI") + scale_colour_manual(values = c("blue", "red"), breaks = c("No","Yes")) +  
  geom_text_repel(aes(label=samp_name), max.overlaps=Inf,  show.legend=FALSE, segment.linetype="dotted"))
invisible(dev.off())

print("All done!")
