echo "performance" > /sys/devices/system/cpu/cpufreq/policy3/scaling_governor

#Ordem das matrizes que serão efetuados os testes para ambos os métodos:
elements=(32 50 64 100 128) # 200 256 300 400 512 1000 1024 2000 2048 3000 4000 4096 5000 10000)

#Escala de valores em x
x_range="[10:10000]"

#Resolução de saída dos gráficos gerados pelo gnuplot.
image_size="1600,900" #1600x900

group="time"
magnitude="milissegundos"

cd t1
make
cd ..

cd t2
make
cd ..

#--------- L2CACHE ------------------

touch l2cache.out
cp /dev/null l2cache.out

# calcula o cache miss
for i in ${elements[@]}; do
	val_old=$(likwid-perfctr -C 3 -g L2CACHE -m -O ./t1/pdeSolver -nx $i -ny $i -i 10 | tail -n1 | cut -d ',' -f2)
	val_new=$(likwid-perfctr -C 3 -g L2CACHE -m -O ./t2/pdeSolver -nx $i -ny $i -i 10 | tail -n1 | cut -d ',' -f2)
	echo "$i $val_old $val_new " >> l2cache.out
done

gnuplot <<- EOF
reset
set terminal png size ${image_size} enhanced font 'Verdana,12' 
set style data linespoints
set style fill solid 2.00 border 0

set style line 1 lc rgb 'orange' lt 1 lw 2 pt 13 ps 1.0
set style line 2 lc rgb 'red' lt 1 lw 2 pt 7 ps 1.0
set style line 3 lc rgb 'grey' lt 1 lw 2 pt 13 ps 1.0
set style line 4 lc rgb 'green' lt 1 lw 2 pt 7 ps 1.0

set key font ",10" 
set key horizontal  
set key spacing 3 
set key samplen 3 
set key width 2

set title 'Medição de desempenho' font ",16" 

set xrange ${x_range}
set logscale x

set xlabel "tamanho da matriz (doubles)" 
set xlabel font ",13" 
set ylabel "$magnitude" 
set ylabel font ",13"

set output "$group-multMatColVet.png"
plot '$group.tmp' using 1:2 with linespoints ls 1 title "old", \
'$group.tmp' using 1:3 with linespoints ls 2 title "new"
EOF