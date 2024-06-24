# Complex-Infections
Contains the scripts and test files for the manuscript "Detecting complex infections in Trypanosomatids using in whole genome sequencing reads"
The script *CI_Estimation_server.v.7.R* receives a VCF file that contains the read depth information for each allele in a sample, estimates the **"Complexity Index" (CI)** and classify the isolate in complex (multiclonal or poliploid) or not.


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

The genome coverage file should be a text file named "Genome_Cov.Table", have no header, and two columns. A TAB separated file.
The first column is the sequence ID and the second the genome coverage, as in:
```
ERR4678145	54
ERR205724	59
ERR205789	61
```






