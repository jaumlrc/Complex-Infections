# Complex-Infections
Contains the scripts and test files for the manuscript "Detecting complex infections in Trypanosomatids using in whole genome sequencing reads"
The script *CI_Estimation_server.v.7.R* receives a VCF file that contains the read depth information for each allele in a sample, estimates the **"Complexity Index" (CI)** and classify the isolate in complex (multiclonal or poliploid) or not.


# How to run the script:
The script receives a VCF file and do all the analysis and plots.

**Running for one sample**:
```
$ #Rscript CI_Estimation_server.v7.R <sample_Id> <vcf.file>
$example:
$ #Rscript CI_Estimation_server.v7.R ERR205724 ERR205724.rec2.recode.vcf
```

**Running for "n" samples listes in the ***sample_names*** file**
```
$ for i in $(cat sample_names); do Rscript CI_Estimation_server.v7.R ${i} ${i}.<vcf_pattern> ; done
$ for i in $(cat sample_names); do Rscript CI_Estimation_server.v7.R ${i} ${i}.rec2.recode.vcf ; done

```

The script can estimate the Chromosomal Copy Number Variation (CCNV) and genome coverage directly from the VCF file. However, we strongly suggest that the CCNV and genome coverage be estimated with other tools. The scipt can






