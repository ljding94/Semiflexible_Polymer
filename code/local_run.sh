L=100
kappa=0.0
f=0.0
g=0.0
for kappa in $(seq 0.0 1.0 20.0)
do
    echo "Running L=$L kappa=$kappa, f=$f, g=$g"
    nohup ./semiflexible_polymer $L $kappa $f $g lo &
done
kappa=10.0
for f in $(seq 0.0 0.1 2.0)
do
    echo "Running L=$L kappa=$kappa, f=$f, g=$g"
    nohup ./semiflexible_polymer $L $kappa $f $g lo &
done
kappa=10.0
f=0.0
for g in $(seq 0.0 0.1 2.0)
do
    echo "Running L=$L kappa=$kappa, f=$f, g=$g"
    nohup ./semiflexible_polymer $L $kappa $f $g lo &
done

