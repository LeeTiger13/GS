########################################################
#PCA分析
#########################################################
cd $workdir  #回到工作目录
mkdir 02.PCA
cd 02.PCA
#方法1：
## plink分析PCA
plink --vcf  $workdir/00.filter/clean.vcf.gz --pca 10 --out  plink_pca   \
    --allow-extra-chr --set-missing-var-ids @:#    --vcf-half-call missing

#绘图

pca_plink_plot.r -i plink_pca.eigenvec -f $GROUP -g Group --name plink_pca