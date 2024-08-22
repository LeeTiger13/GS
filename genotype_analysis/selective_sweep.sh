######################################################
#选择消除分析  fst pi ROD  基于多样性降低
########################################################
cd $workdir  #回到工作目录
mkdir 05.select_sweep
cd 05.select_sweep

#不同群体的样本分开放到不同的文件中
cat $GROUP |grep wild |cut -f 1 >wild_popid.txt
cat $GROUP |grep cultivated|grep -v "Indian"|grep -v "Xishuangbanna" |cut -f 1 >cultivated_popid.txt

#fst  pi  tajimaD分析
mkdir fst_pi_ROD
cd fst_pi_ROD

#设置输入的vcf文件与计算的窗口和步长
gzvcf=$workdir/00.filter/clean.vcf.gz
window=100000
step=10000

#pi 多样性
vcftools  --gzvcf $gzvcf \
    --window-pi $window --window-pi-step  $step  \
    --keep ../wild_popid.txt   --out pi.wild
vcftools  --gzvcf $gzvcf \
    --window-pi $window --window-pi-step  $step  \
    --keep ../cultivated_popid.txt  --out pi.cultivated

#pi多样性绘图输出
pi_manhattan_plot.r -i pi.wild.windowed.pi -F $FAI -f 19226500 -n pi.wild
pi_smooth_line_plot.r -i pi.wild.windowed.pi -F $FAI -f 19226500  -n pi.wild.smoothline

#Fst  群体间多样性差异
vcftools  --gzvcf $gzvcf --fst-window-size $window --fst-window-step $step  \
    --weir-fst-pop  ../wild_popid.txt --weir-fst-pop ../cultivated_popid.txt --out  Fst.wild.cultivated

#fst绘图输出
fst_manhattan_plot.r -i Fst.wild.cultivated.windowed.weir.fst -F $FAI -f 19226500 -n Fst.wild.cultivated
fst_manhattan_plot.r -i Fst.wild.cultivated.windowed.weir.fst -F $FAI -f 19226500  -n Fst.wild.cultivated_vline --vline
fst_smooth_line_plot.r -i Fst.wild.cultivated.windowed.weir.fst -F $FAI -f 19226500  -n Fst.wild.cultivated_smoothline

#fst 与 pi 联合 筛选受选择区域
fst_pi_select_sweep.r --fst Fst.wild.cultivated.windowed.weir.fst \
    --pi1 pi.wild.windowed.pi --pi2 pi.cultivated.windowed.pi --zscore --log2 \
    -A wild -B cultivated -c 0.05 -n fst-pi.wild-vs-cultivated  -f pdf

fst_pi_select_sweep.r --fst Fst.wild.cultivated.windowed.weir.fst \
    --pi1 pi.wild.windowed.pi --pi2 pi.cultivated.windowed.pi --zscore --log2 \
    -A wild -B cultivated -c 0.05 -n fst-pi.wild-vs-cultivated  -f png

#ROD计算与展示
rod_calculate.r --wild pi.wild.windowed.pi --domesticated pi.cultivated.windowed.pi -p ROD.wild.cultivated

#绘图输出
rod_manhattan_plot.r -i ROD.wild.cultivated.txt -F $FAI -f 19226500  -n ROD.wild.cultivated 
rod_smooth_line_plot.r -i ROD.wild.cultivated.txt -F $FAI -f 19226500  -n ROD.wild.cultivated.smoothline