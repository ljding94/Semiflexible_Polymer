L=100
kappa=0.0
f=0.0
g=0.0
#for kappa in $(seq 0.0 1.0 20.0)
#do
#    echo "Running L=$L kappa=$kappa, f=$f, g=$g"
#    nohup ./semiflexible_polymer $L $kappa $f $g lo &
#done
kappa=10.0
for f in 0.00 0.10 0.30
do
    for gL in 0.00 0.50 1.50
    do
        echo "Running L=$L kappa=$kappa, f=$f, gL=$gL"
        nohup ./semiflexible_polymer $L $kappa $f $gL 0 ~/Work/Semiflexible_Polymer/data/scratch_local/data_pool &
    done
done