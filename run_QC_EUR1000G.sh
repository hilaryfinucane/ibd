CHR=$1
POP="EUR"

outdir=/groups/price/hilary/ibd/data/

if [ "$CHR" == "" ]; then
    echo "Error: Need CHR as input argument"
    exit
fi

# wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr$CHR.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
# gunzip ALL.chr$CHR.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz

# cp /groups/price/ldsc/reference_files/1000G_EUR_Phase3/TGP2261.txt $outdir
# awk '{if ($2=="EUR") {print $1, $1} }' ${outdir}TGP2261.txt > ${outdir}list_id.txt

# /groups/price/steven/soft/plink2_v1.90b3w/plink \
#   --vcf ALL.chr$CHR.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf \
#   --keep ${outdir}list_id.txt \
#   --geno 0.01 \
#   --hwe 0.001 \
#   --biallelic-only \
#   --snps-only \
#   --make-bed \
#   --out ${outdir}1000G.$POP.temp.$CHR
  
# #remove duplicate SNPs
# perl remove_duplicate_snps.pl ${outdir}1000G.$POP.temp.$CHR
# /groups/price/steven/soft/plink2_v1.90b3w/plink \
#   --bfile ${outdir}1000G.$POP.temp.$CHR \
#   --exclude ${outdir}1000G.$POP.temp.$CHR.list \
#   --make-bed \
#   --out ${outdir}1000G.$POP.temp2.$CHR

# /groups/price/steven/soft/plink2_v1.90b3w/plink \
#   --bfile ${outdir}1000G.$POP.temp2.$CHR \
#   --cm-map /groups/price/steven/data/recombination_map/genetic_map_b37/genetic_map_chr${CHR}_combined_b37.txt $CHR \
#   --zero-cms \
#   --make-bed \
#   --out ${outdir}1000G.$POP.QC.$CHR

# /groups/price/steven/soft/plink2_v1.90b3w/plink \
#   --bfile ${outdir}1000G.$POP.QC.$CHR \
#   --freq \
#   --out ${outdir}1000G.$POP.QC.$CHR
  
rm ALL.chr$CHR.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf
rm ${outdir}1000G.$POP.temp*.$CHR.*

