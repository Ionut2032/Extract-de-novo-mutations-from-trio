##########################################################################################################################
### INSTALL PACKAGES AND LOAD DATA
##########################################################################################################################

#install & load vcfR package (useful for reading the .vcf file)
install.packages("vcfR")
library("vcfR")

#read the vcf file as a vcfR object
multisample_vcf_trio <- read.vcfR("filtered_trios.vcf")

#extract the genotype as reference vs. al1ternative (e.g. 0/1)
extracted_genotype_ref_vs_alt <- extract.gt(multisample_vcf_trio, element="GT")


#read the matrix as a data frame so it can be easily accessed from the name of the column
data_frame <- as.data.frame(extracted_genotype_ref_vs_alt, row.names = NULL, stringsAsFactors=FALSE)

#extract the DNA bases i.e. alleles from reference vs. alternative allele for x3 samples
ref_vs_alt_x3samples <- extract.haps(multisample_vcf_trio, mask = FALSE, unphased_as_NA = FALSE)
#read the matrix as a data frame so it can be easily accessed from the name of the column
data_frame_genotype <- as.data.frame(ref_vs_alt_x3samples, row.names = NULL, stringsAsFactors=FALSE)

########################################################################################################################
### COUNT THE NUMBER OF INCOMPATIBLE VARIANTS (i.e. which are parents=highest no, which is the child=the other one?)
########################################################################################################################
NA12878_VS_NA12891 <- with(data_frame, data_frame[(data_frame$NA12878 == "0/0" & data_frame$NA12891 == "1/1") | (data_frame$NA12878 == "1/1" & data_frame$NA12891 == "0/0"), ])
nrow(NA12878_VS_NA12891)

NA12878_VS_NA12892 <- with(data_frame, data_frame[data_frame$NA12878 == "0/0" & data_frame$NA12892 == "1/1" | data_frame$NA12878 == "1/1" & data_frame$NA12892 == "0/0", ])
nrow(NA12878_VS_NA12892)

NA12891_VS_NA12892 <- with(data_frame, data_frame[data_frame$NA12891 == "0/0" & data_frame$NA12892 == "1/1" | data_frame$NA12891 == "1/1" & data_frame$NA12892 == "0/0", ])
nrow(NA12891_VS_NA12892)

########################################################################################################################
### FIND DE-NOVO VARIANTS
########################################################################################################################

#find de-novo 0/0 parents & 0/1 child
de_novo_1 <- with(data_frame, data_frame[(data_frame$NA12891 == "0/0" & data_frame$NA12892 == "0/0" 
                                          & data_frame$NA12878 == "0/1"), ])
#find de-novo 0/0 parents & 1/1 child
de_novo_2 <- with(data_frame, data_frame[(data_frame$NA12891 == "0/0" & data_frame$NA12892 == "0/0" 
                                          & data_frame$NA12878 == "1/1"), ])
#find total de-novos
de_novo <- with(data_frame, data_frame[(data_frame$NA12891 == "0/0" & data_frame$NA12892 == "0/0" & data_frame$NA12878 == "0/1" | data_frame$NA12891 == "0/0" & data_frame$NA12892 == "0/0" & data_frame$NA12878 == "1/1"), ])
nrow(de_novo)

#find the types of DNA variations present in the de-novo variations
#put the chromosome & position as different columns
library(data.table)
library(plyr)

#include the row names as a column
hh1 <- setDT(data_frame_genotype, keep.rownames = TRUE)
hh2 <- setDT(de_novo, keep.rownames = TRUE)
#merge the variants with the genotype (i.e. type of mutation present)
dd <- subset(hh1, rn %in% hh2$rn)

#further reduce the number of variants to the 40 validated by Sanger sequencing & their genotype
validated_de_novos <- dd[dd$rn %in% c("1_75884343", "1_110583335", "1_182974758", "2_39556621", "2_152899032",
                                      "2_182693277", "3_101454745", "3_118900031", "4_15227519", "4_104624818",
                                      "5_6466106", "5_126385924", "5_145107247", "6_52120843", "6_145808310",
                                      "6_160334960", "8_21568355", "8_74680107", "9_38096405", "9_89775629",
                                      "9_123350598", "10_56256294", "10_120325447", "11_40851625", "11_76051551",
                                      "11_85631211", "11_117264626", "14_56763353", "14_78809184", "15_50953965",
                                      "15_51248561", "15_58669774", "15_85715561", "15_87791919", "15_99175976",
                                      "16_29601086", "16_82540504", "17_53502801", "17_71375843", "17_80712655",
                                      "18_74117587", "19_6661912", "20_55356548", "21_20165354", "22_24262476",
                                      "X_8598739")]

