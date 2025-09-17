#1) Get dbSNP GRCh37 VCF (once)
# (On NIH Biowulf: module load bcftools)
wget -nc https://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/00-All.vcf.gz
wget -nc https://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/00-All.vcf.gz.tbi

#2) Extract rsIDs from your table
cut -f1 gwas_b38.tsv | sed '1d' | grep -E '^rs[0-9]+' | sort -u > rsids.txt

#3) Pull b37 coordinates for just those rsIDs
bcftools view -i 'ID=@rsids.txt' 00-All.vcf.gz -Ou \
| bcftools query -f '%ID\t%CHROM\t%POS\n' > rsid_to_b37.tsv
# columns: SNP  Chr_b37  Pos_b37

#4) Left-join back to your file (keeps your P)
awk 'BEGIN{FS=OFS="\t"} FNR==NR{m[$1]=$2"\t"$3; next}
     NR==1{print $0,"Chr_b37","Pos_b37"; next}
     {split(m[$1],a,"\t"); print $0,(m[$1]?a[1]: "NA"),(m[$1]?a[2]:"NA")}' \
     rsid_to_b37.tsv gwas_b38.tsv > gwas_b37.tsv


##
awk 'BEGIN{FS=OFS="\t"; print "SNP","Chr","Pos","P"} NR>1 {print $1,$5,$6,$4}' gwas_b37.tsv > GCST90027158_buildGRCh37.tsv.gz.filtered.txt
