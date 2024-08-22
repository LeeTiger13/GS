# population structure analysis

## filter LD 

#50 10 0.2   50个SNP 窗口  step 10个  r2 0.2  ; 50 5 0.4
plink --vcf  $workdir/00.filter/clean.vcf.gz  --indep-pairwise 50 10 0.2 --out ld   \
    --allow-extra-chr --set-missing-var-ids @:# 
plink --vcf  $workdir/00.filter/clean.vcf.gz  --make-bed --extract ld.prune.in  \
    --out LDfiltered --recode vcf-iid  --keep-allele-order  --allow-extra-chr --set-missing-var-ids @:#  

#转换成plink格式
vcftools --vcf LDfiltered.vcf --plink \
    --out plink
#转换成admixture要求的bed格式
plink --noweb --file plink  --recode12 --out admixture \
     --allow-extra-chr  --keep-allele-order

#admixture 群体结构分析
for k in {2..10};do
    admixture -j2 -C 0.01 --cv admixture.ped $k >admixture.log$k.out
done

#绘图展示 
structure_plot.r  -d ./ -s admixture.nosex 
structure_plot.r  -d ./ -s admixture.nosex -f $GROUP -g Group  #按照分组顺序显示

#确定最佳K，CV值最小时对应的K值为最佳K

grep "CV error" *out