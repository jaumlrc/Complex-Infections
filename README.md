# Complex-Infections
Contains the scripts and test files for the manuscript "Detecting complex infections in Trypanosomatids using in whole genome sequencing reads"   
   
The script *CI_Estimation_server.v.7.R* receives a VCF file that contains the read depth information for each allele in a sample, estimates the **"Complexity Index" (CI)** and classify the isolate in complex (multiclonal or poliploid) or not.   
   
Examples of SNP callers that can generate these vcfs are [GATK](https://gatk.broadinstitute.org/hc/en-us), [freebayes](https://github.com/freebayes/freebayes) and [octopus](https://luntergroup.github.io/octopus/).  
   
We **Strongly** suggest that the VCF file be filtered to maintain only biallelic SNP sites, removing insertions/deletions.   


# Required R libraries:
* ggplot2
* reshape2
* viridis
* dplyr
* tidyr
* vcfR
* ggrepel

These can be installed in R with:
```
$  install.packages('<package_name')
as in
$  install.packages('vcfR')
```

# How to run the script:
The script receives a VCF file and do all the analysis and plots.

**Running for one sample**:
```
$  Rscript CI_Estimation_server.v7.R <sample_Id> <vcf.file>
Example:
$  Rscript CI_Estimation_server.v7.R ERR205724 ERR205724.rec2.recode.vcf
```

**Running for "n" samples listes in the ***sample_names*** file**
```
$  for i in $(cat sample_names); do Rscript CI_Estimation_server.v7.R ${i} ${i}.<vcf_pattern> ; done
Example
$  for i in $(cat sample_names); do Rscript CI_Estimation_server.v7.R ${i} ${i}.rec2.recode.vcf ; done
```

**Supporting data:**
The "sample_names" file is a text file containing the IDs of your samples, one at each line, as in:
```
ERR205724
ERR205789
ERR3956088
ERR4678145
```
The script can estimate the Chromosomal Copy Number Variation (CCNV) and genome coverage directly from the VCF file.<br /> 
However, we strongly suggest that the CCNV and genome coverage be estimated with other tools.
The script automatically detect the files ***Genome_Cov.Table*** (has the genome coverage information) and ***CCNV.Table.ordered*** (has the CCNV information) in the same folder that the script is run.

The genome coverage file should be a TAB sepparated file file named "Genome_Cov.Table", have no header, and two columns. 
The first column is the sequence ID and the second the genome coverage, as in:
```
ERR4678145	54
ERR205724	59
ERR205789	61
```

The Chromosome copies file has to be named "CCNV.Table.ordered", also TAB separatted file
The first column correspond to the chromosome number, similar to what was present in the VCF file, and the header must be "Chromosome"
Chromosome names should have only numerals, as in 1, 2, 3, so that the plots always have the chromosomes in order.
The headers from columns 2 to "n" must be the sample IDs, identical to the ones in the "sample_names" file.
The lines correspond to each chromosome copy number
And example of the "CCNV.Table.ordered" file for the isolates ERR205724, ERR205789, ERR3956088 and ERR4678145 from L. donovani, that has 36 chromosomes

```
Chromosome	ERR205724	ERR205789	ERR3956088	ERR4678145
1	1	0.9576	1.066	1	1.019
2	2	1.051	1.131	0.9091	1.046
3	3	0.9661	0.8033	1.045	1.056
4	4	1.339	1.131	1	1
5	5	0.9661	1.344	1.045	1.037
6	6	1.271	0.8033	1.045	1.037
7	7	0.9831	1.082	1.045	1.019
8	8	1.237	1.131	1	0.9907
9	9	1.051	0.8033	1.045	1.019
10	10	1.271	0.8033	1.045	1.019
11	11	1.271	1.131	1.045	1.037
12	12	0.9831	1.098	1.023	1.037
13	13	0.9831	1.115	1.045	1.019
14	14	0.9831	1.098	1	1
15	15	0.9661	1.115	1	1.009
16	16	0.9831	0.8197	1.045	1.019
17	17	0.9831	0.8197	1.045	1.037
18	18	0.9831	0.8525	1	1.019
19	19	0.9831	0.8033	1	1.019
20	20	0.9746	1.082	1	1
21	21	0.9831	0.8197	1.045	1
22	22	1.051	1.131	1	1.019
23	23	1.068	1.148	1.136	1.019
24	24	1.051	1.082	1.045	1.019
25	25	0.9661	0.8033	1	1
26	26	1.203	0.8197	1	1
27	27	0.9661	0.8197	1	1
28	28	0.9661	0.8197	1	1
29	29	1.322	0.8197	1	1
30	30	0.9661	0.8197	1	1
31	31	1.695	1.41	1.909	1.778
32	32	0.9661	1.082	1	1
33	33	1.271	1.361	1	1
34	34	0.9831	0.8197	1	0.963
35	35	1	1.148	1	0.9815
36	36	0.9661	1.082	0.9545	0.9815
```

# Test data and expected outpus
Example data for 4 VCFs can be obteined in the [Test_data](https://github.com/jaumlrc/Complex-Infections/tree/main/Test_data) folder in the repository.   
This folder also has:   
A example [CCNV.Table.ordered](https://github.com/jaumlrc/Complex-Infections/blob/main/Test_data/CCNV.Table.ordered) CCNV data file;   
A example [Genome_Cov.Table](https://github.com/jaumlrc/Complex-Infections/blob/main/Test_data/Genome_Cov.Table) Genome coverage file   
And a exemple [sample_names](https://github.com/jaumlrc/Complex-Infections/blob/main/Test_data/sample_names) file with sample IDs   

The output files generated for this four isolates can be seen in the folder (Test_outputs)[https://github.com/jaumlrc/Complex-Infections/tree/main/Test_outputs]

#Output files







