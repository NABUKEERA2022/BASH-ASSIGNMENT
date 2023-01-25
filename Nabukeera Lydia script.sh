#!/bin/bash

# Files to use can be downloaded from https://drive.google.com/drive/folders/11UD52i99CaCSBEJFNb8Y1afo9p3hL8cL?usp=sharing

# Unzip the vcf file using the "gunzip" code

gunzip sample.vcf.gz 

# Qn 3: Number of samples in the file
    bcftools query -l sample.vcf.gz | wc -l
    
# Qn 4. Number of variants in the file
    bcftools query -f '%ALT\n' sample.vcf | wc -l


# Qn 5. How would you extract the chromosome, position, QualByDepth and RMSMappingQuality fields? 
#Save the output to a tab-delimited file
    bcftools query -f '%CHROM\t%POS\t%INFO/QD\t%INFO/MQ\n' sample.vcf > QR.txt

# Qn 6. Extract data that belongs to chromosomes 2,4 and MT
    awk '$1 ~ /^(2|4|MT)$/ {print $0}' sample.vcf > chrom24MT.vcf

# Qn7. Print out variants that do not belong to chr20:1-30000000
    awk '$1 != "20" || ($1 == "chr20" && ($2 < 1 || $2 > 30000000)) \
        {print $1, $2, $4, $5}' sample.vcf > Chr20.vcf

# Qn8. Extract variants that belong to SRR13107019
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' -s SRR13107019 sample.vcf > SRR19variant.txt 

# Qn9. Filter out variants with a QualByDepth above 7
 vcftools --vcf sample.vcf --minDP 7 --recode --out QD7

#Qn 10 How many contigs are referred to in the file. Check the header section
grep -c "^##contig" sample.vcf

#Qn 12 Extract data on the read depth of called variants for sample SRR13107018
bcftools query -f '%DP\n' -s SRR13107018 sample.vcf > readSRR18

#Qn 13 Extract data on the allele frequency of alternate alleles. Combine this data with the
#chromosome and position of the alternate allele

bcftools query -f '%CHROM\t%POS\t%AF\n' sample.vcf > ALTallele.txt

# Manipulating SAM files

#Qn 3: How many samples are in the file

grep -E '^@RG' sample.sam | awk '{print $2}' | awk -F ':' '{print $2}' | sort | uniq | wc -l

# Qn 4: How many alignments are in the file

grep -c "^@" sample.sam #I got 35511 alignments.

#Qn 5: Get summary statistics for the alignments in the file
samtools flagstat sample.sam > samstat.txt

# Qn 6: Count the number of fields in the file
awk '{print NF}' sample.sam

#Qn 7: Print all lines in the file that have @SQ and sequence name tag beginning with NT_

grep '@SQ.*NT_' sample.sam 

#Qn 8: Print all lines in the file that have @RG and LB tag beginning with Solexa

grep "@RG.*LB:Solexa" sample.sam

# Qn 9: Extract primarily aligned sequences and save them in another file

awk '$1 !~ /^@/ && $2 == "99" || $2 == "83"' sample.sam > savedalignedseq.sam

#10. Extract alignments that map to chromosomes 1 and 3. Save the output in BAM
#format

awk '$1 !~ /^@/ && ($3 == "1" || $3 == "3")' sample.sam | samtools view -Sb - > alignmap.bam

#Qn: 11. How would you obtain unmapped reads from the file
samtools view -f 4 sample.sam > unmapedreads.sam

# Qn 12. How many reads are aligned to chromosome 4
grep -c "^4\t" sample.sam #Gave an output of 0 reads

#Qn 14. Extract all optional fields of the file and save them in “optional_fields.txt”

awk '{for(i=11;i<=NF;i++) print $i}' sample.sam > optional_fields.txt
