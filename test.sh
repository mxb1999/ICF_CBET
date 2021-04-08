echo "Testing CUDA Ray Tracing"
TRIALS=3
#for i in {51,101,151,201,251,301,351}
#do
#    for ((t=0;t<$TRIALS;t++));
#    do
#        ./implSim 0 12 1 $i
#    
#done
#echo "Testing CPU Ray Tracing"
#for j in {1,12}
#do
#    for i in {51,101,151,201,251,301,351}
#    do
#        for ((t=0;t<$TRIALS;t++));
#        do
#            ./implSim 0 $j 0 $i
#        done
#    done
#done
echo "Testing CUDA CBET"
for i in {5,10,20,30,40,50,60,70,80,90,110,120,130,140,150}
do
    for ((t=0;t<$TRIALS;t++));
    do
        ./implSim 0 12 1 $i
    done
done
echo "Testing CPU CBET"
for j in {1,12}
do
    for i in {5,10,20,30,40,50,60,70,80,90,110,120,130,140,150}
    do
        for ((t=0;t<$TRIALS;t++));
        do
            ./implSim 0 $j 0 $i
        done
    done
done