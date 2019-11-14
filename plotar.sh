echo "performance" > /sys/devices/system/cpu/cpufreq/policy3/scaling_governor

#Ordem das matrizes que serão efetuados os testes para ambos os métodos:
elements=(32, 50, 64, 100, 128, 200, 256, 300, 400, 512, 1000, 1024, 2000, 2048, 3000, 4000, 4096, 5000, 10000)

#Escala de valores em x
x_range="[32:10000]"

#Resolução de saída dos gráficos gerados pelo gnuplot.
image_size="1600,900" #1600x900

cd t1
make
cd ..

cd t2
make
cd ..

#--------- L2CACHE ------------------

echo "Calculando L2CACHE"

arq_tmp=l2cache.out

touch $arq_tmp
cp /dev/null $arq_tmp

# calcula o cache miss
for i in ${elements[@]}; do
	val_old=$(likwid-perfctr -C 3 -g L2CACHE -m -O ./t1/pdeSolver -nx $i -ny $i -i 10 | tail -n1 | cut -d ',' -f2)
	val_new=$(likwid-perfctr -C 3 -g L2CACHE -m -O ./t2/pdeSolver -nx $i -ny $i -i 10 | tail -n1 | cut -d ',' -f2)
	echo "$i $val_old $val_new " >> $arq_tmp
done

ytag="L2 miss ratio"
arq_saida="l2missratio.png"

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

set title 'Medição de cache miss L2' font ",16" 

set xrange ${x_range}

set xlabel "Tamanho de nx e ny" 
set xlabel font ",13" 
set ylabel "$ytag" 
set ylabel font ",13"

set output "$arq_saida"
plot '$arq_tmp' using 1:2 with linespoints ls 1 title "Trabalho 1", \
'$arq_tmp' using 1:3 with linespoints ls 2 title "Trabalho 2 (Otimizado)"
EOF
echo "Grafico plotado em $arq_saida"

# ---------- L3 BANDIWIDTH ----------------------------------
echo "Calculando L3 BANDIWIDTH"

arq_tmp=l3bandwidth.out

touch $arq_tmp
cp /dev/null $arq_tmp

# calcula o cache miss
for i in ${elements[@]}; do
	val_old=$(likwid-perfctr -C 3 -g L3 -m -O ./t1/pdeSolver -nx $i -ny $i -i 10 | tail -n2 | head -n1 | cut -d ',' -f2)
	val_new=$(likwid-perfctr -C 3 -g L3 -m -O ./t2/pdeSolver -nx $i -ny $i -i 10 | tail -n2 | head -n1 | cut -d ',' -f2)
	echo "$i $val_old $val_new " >> $arq_tmp
done

ytag="Memory bandwidth [MBytes/s]"
arq_saida="l3bandwidth.png"

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

set title 'Medição de bandwidth L3' font ",16" 

set xrange ${x_range}

set xlabel "Tamanho de nx e ny" 
set xlabel font ",13" 
set ylabel "$ytag" 
set ylabel font ",13"

set output "$arq_saida"
plot '$arq_tmp' using 1:2 with linespoints ls 1 title "Trabalho 1", \
'$arq_tmp' using 1:3 with linespoints ls 2 title "Trabalho 2 (Otimizado)"
EOF
echo "Grafico plotado em $arq_saida"

# ---------- FLOPS_DP e AVX ----------------------------------
echo "Calculando FLOPS DP e AVX"

arq_tmp=flops.out

touch $arq_tmp
cp /dev/null $arq_tmp

# calcula o cache miss
for i in ${elements[@]}; do
	val_old_dp=$(likwid-perfctr -C 3 -g L3 -m -O ./t1/pdeSolver -nx $i -ny $i -i 10 | tail -n5 | head -n1 | cut -d ',' -f2)
	val_new_dp=$(likwid-perfctr -C 3 -g L3 -m -O ./t2/pdeSolver -nx $i -ny $i -i 10 | tail -n5 | head -n1 | cut -d ',' -f2)
	val_old_avx=$(likwid-perfctr -C 3 -g L3 -m -O ./t1/pdeSolver -nx $i -ny $i -i 10 | tail -n4 | head -n1 | cut -d ',' -f2)
	val_new_avx=$(likwid-perfctr -C 3 -g L3 -m -O ./t2/pdeSolver -nx $i -ny $i -i 10 | tail -n4 | head -n1 | cut -d ',' -f2)

	echo "$i $val_old_dp $val_new_dp $val_old_avx $val_new_avx" >> $arq_tmp
done

ytag="MFLOP/s"
arq_saida="mflops.png"

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

set title 'Medição de FLOPS DP e FLOPS AVX ' font ",16" 

set xrange ${x_range}

set xlabel "Tamanho de nx e ny" 
set xlabel font ",13" 
set ylabel "$ytag" 
set ylabel font ",13"

set output "$arq_saida"
plot '$arq_tmp' using 1:2 with linespoints ls 1 title "Trabalho 1 - FLOPS DP", \
'$arq_tmp' using 1:3 with linespoints ls 2 title "Trabalho 2 - FLOPS DP", \
'$arq_tmp' using 1:4 with linespoints ls 3 title "Trabalho 1 - AVX DP", \
'$arq_tmp' using 1:5 with linespoints ls 4 title "Trabalho 2 - AVX DP"
EOF
echo "Grafico plotado em $arq_saida"

#--------- TEMPO ------------------

echo "Calculando tempo de execucao"

arq_tmp=time.out

touch $arq_tmp
cp /dev/null $arq_tmp

# calcula o cache miss
for i in ${elements[@]}; do
	val_old=$(likwid-perfctr -C 3 -g L3 -m -O ./t1/pdeSolver -nx $i -ny $i -i 10 | tail -n10 | head -n1 | cut -d ',' -f2)
	val_new=$(likwid-perfctr -C 3 -g L3 -m -O ./t2/pdeSolver -nx $i -ny $i -i 10 | tail -n10 | head -n1 | cut -d ',' -f2)
	echo "$i $val_old $val_new " >> $arq_tmp
done

ytag="Tempo de execução [ms]"
arq_saida="time.png"

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
set logscale x

set title 'Medição de tempo de execução' font ",16" 

set xrange ${x_range}

set xlabel "Tamanho de nx e ny" 
set xlabel font ",13" 
set ylabel "$ytag" 
set ylabel font ",13"

set output "$arq_saida"
plot '$arq_tmp' using 1:2 with linespoints ls 1 title "Trabalho 1", \
'$arq_tmp' using 1:3 with linespoints ls 2 title "Trabalho 2 (Otimizado)"
EOF
echo "Grafico plotado em $arq_saida"
