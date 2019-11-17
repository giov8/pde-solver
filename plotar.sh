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

#--------- TEMPO (Gauss) ------------------

echo "Calculando tempo de execucao (Gauss-Seidel)"

arq_tmp=time_gauss.out

touch $arq_tmp
cp /dev/null $arq_tmp
aux=10

for i in ${elements[@]}; do
	val_old=$(likwid-perfctr -C 3 -g L3 -m -O ./t1/pdeSolver -nx $i -ny $i -i 10 | tail -n36 | head -n1 | cut -d ',' -f2)
	val_new=$(likwid-perfctr -C 3 -g L3 -m -O ./t2/pdeSolver -nx $i -ny $i -i 10 | tail -n36 | head -n1 | cut -d ',' -f2)
	val_old=$(echo "scale=5; $val_old / $aux" | bc)
	val_new=$(echo "scale=5; $val_new / $aux" | bc)
	echo "$i $val_old $val_new " >> $arq_tmp
done

ytag="Tempo de execução [s]"
arq_saida="time_gauss.png"

gnuplot <<- EOF
reset
set terminal png size ${image_size} enhanced font 'Verdana,18' 
set style data linespoints
set style fill solid 2.00 border 0

set style line 1 lc rgb 'orange' lt 1 lw 2 pt 13 ps 1.0
set style line 2 lc rgb 'red' lt 1 lw 2 pt 7 ps 1.0
set style line 3 lc rgb 'grey' lt 1 lw 2 pt 13 ps 1.0
set style line 4 lc rgb 'green' lt 1 lw 2 pt 7 ps 1.0

set key font ",18" 
set key horizontal  
set key spacing 3 
set key samplen 3 
set key width 2
set logscale y

set title 'Medição de tempo de execução (Gauss-Seidel)' font ",24" 

set xrange ${x_range}

set xlabel "Tamanho de nx e ny" 
set xlabel font ",20" 
set ylabel "$ytag" 
set ylabel font ",20"

set output "$arq_saida"
plot '$arq_tmp' using 1:2 with linespoints ls 1 title "Trabalho 1", \
'$arq_tmp' using 1:3 with linespoints ls 2 title "Trabalho 2 (Otimizado)"
EOF
echo "Grafico plotado em $arq_saida"

#--------- TEMPO (Residuo) ------------------

echo "Calculando tempo de execucao (Residuo)"

arq_tmp=time_residuo.out

touch $arq_tmp
cp /dev/null $arq_tmp
aux=10

for i in ${elements[@]}; do
	val_old=$(likwid-perfctr -C 3 -g L3 -m -O ./t1/pdeSolver -nx $i -ny $i -i 10 | tail -n10 | head -n1 | cut -d ',' -f2)
	val_new=$(likwid-perfctr -C 3 -g L3 -m -O ./t2/pdeSolver -nx $i -ny $i -i 10 | tail -n10 | head -n1 | cut -d ',' -f2)
	val_old=$(echo "scale=5; $val_old / $aux" | bc)
	val_new=$(echo "scale=5; $val_new / $aux" | bc)
	echo "$i $val_old $val_new " >> $arq_tmp
done

ytag="Tempo de execução [s]"
arq_saida="time_residuo.png"

gnuplot <<- EOF
reset
set terminal png size ${image_size} enhanced font 'Verdana,18' 
set style data linespoints
set style fill solid 2.00 border 0

set style line 1 lc rgb 'orange' lt 1 lw 2 pt 13 ps 1.0
set style line 2 lc rgb 'red' lt 1 lw 2 pt 7 ps 1.0
set style line 3 lc rgb 'grey' lt 1 lw 2 pt 13 ps 1.0
set style line 4 lc rgb 'green' lt 1 lw 2 pt 7 ps 1.0

set key font ",18" 
set key horizontal  
set key spacing 3 
set key samplen 3 
set key width 2
set logscale y

set title 'Medição de tempo de execução (Residuo)' font ",24" 

set xrange ${x_range}

set xlabel "Tamanho de nx e ny" 
set xlabel font ",20" 
set ylabel "$ytag" 
set ylabel font ",20"

set output "$arq_saida"
plot '$arq_tmp' using 1:2 with linespoints ls 1 title "Trabalho 1", \
'$arq_tmp' using 1:3 with linespoints ls 2 title "Trabalho 2 (Otimizado)"
EOF
echo "Grafico plotado em $arq_saida"

#--------- L2CACHE (Gauss) ------------------

echo "Calculando L2CACHE (Gauss-Seidel)"

arq_tmp=l2cache_gauss.out

touch $arq_tmp
cp /dev/null $arq_tmp


for i in ${elements[@]}; do
	val_old=$(likwid-perfctr -C 3 -g L2CACHE -m -O ./t1/pdeSolver -nx $i -ny $i -i 10 | tail -n24 | head -n1 | cut -d ',' -f2)
	val_new=$(likwid-perfctr -C 3 -g L2CACHE -m -O ./t2/pdeSolver -nx $i -ny $i -i 10 | tail -n24 | head -n1 | cut -d ',' -f2)
	echo "$i $val_old $val_new " >> $arq_tmp
done

ytag="L2 miss ratio"
arq_saida="l2missratio_gauss.png"

gnuplot <<- EOF
reset
set terminal png size ${image_size} enhanced font 'Verdana,18' 
set style data linespoints
set style fill solid 2.00 border 0

set style line 1 lc rgb 'orange' lt 1 lw 2 pt 13 ps 1.0
set style line 2 lc rgb 'red' lt 1 lw 2 pt 7 ps 1.0
set style line 3 lc rgb 'grey' lt 1 lw 2 pt 13 ps 1.0
set style line 4 lc rgb 'green' lt 1 lw 2 pt 7 ps 1.0

set key font ",18" 
set key horizontal  
set key spacing 3 
set key samplen 3 
set key width 2

set title 'Medição de cache miss L2 (Gauss-Seidel)' font ",24" 

set xrange ${x_range}

set xlabel "Tamanho de nx e ny" 
set xlabel font ",20" 
set ylabel "$ytag" 
set ylabel font ",20"

set output "$arq_saida"
plot '$arq_tmp' using 1:2 with linespoints ls 1 title "Trabalho 1", \
'$arq_tmp' using 1:3 with linespoints ls 2 title "Trabalho 2 (Otimizado)"
EOF
echo "Grafico plotado em $arq_saida"

#--------- L2CACHE (Residuo) ------------------

echo "Calculando L2CACHE (Residuo)"

arq_tmp=l2cache_residuo.out

touch $arq_tmp
cp /dev/null $arq_tmp


for i in ${elements[@]}; do
	val_old=$(likwid-perfctr -C 3 -g L2CACHE -m -O ./t1/pdeSolver -nx $i -ny $i -i 10 | tail -n1 | cut -d ',' -f2)
	val_new=$(likwid-perfctr -C 3 -g L2CACHE -m -O ./t2/pdeSolver -nx $i -ny $i -i 10 | tail -n1 | cut -d ',' -f2)
	echo "$i $val_old $val_new " >> $arq_tmp
done

ytag="L2 miss ratio"
arq_saida="l2missratio_residuo.png"

gnuplot <<- EOF
reset
set terminal png size ${image_size} enhanced font 'Verdana,18' 
set style data linespoints
set style fill solid 2.00 border 0

set style line 1 lc rgb 'orange' lt 1 lw 2 pt 13 ps 1.0
set style line 2 lc rgb 'red' lt 1 lw 2 pt 7 ps 1.0
set style line 3 lc rgb 'grey' lt 1 lw 2 pt 13 ps 1.0
set style line 4 lc rgb 'green' lt 1 lw 2 pt 7 ps 1.0

set key font ",18" 
set key horizontal  
set key spacing 3 
set key samplen 3 
set key width 2

set title 'Medição de cache miss L2 (Resíduo)' font ",24" 

set xrange ${x_range}

set xlabel "Tamanho de nx e ny" 
set xlabel font ",20" 
set ylabel "$ytag" 
set ylabel font ",20"

set output "$arq_saida"
plot '$arq_tmp' using 1:2 with linespoints ls 1 title "Trabalho 1", \
'$arq_tmp' using 1:3 with linespoints ls 2 title "Trabalho 2 (Otimizado)"
EOF
echo "Grafico plotado em $arq_saida"


# ---------- L3 BANDWIDTH (Gauss) ----------------------------------
echo "Calculando L3 BANDIWIDTH (Gauss-Seidel)"

arq_tmp=l3bandwidth_gauss.out

touch $arq_tmp
cp /dev/null $arq_tmp


for i in ${elements[@]}; do
	val_old=$(likwid-perfctr -C 3 -g L3 -m -O ./t1/pdeSolver -nx $i -ny $i -i 10 | tail -n28 | head -n1 | cut -d ',' -f2)
	val_new=$(likwid-perfctr -C 3 -g L3 -m -O ./t2/pdeSolver -nx $i -ny $i -i 10 | tail -n28 | head -n1 | cut -d ',' -f2)
	echo "$i $val_old $val_new " >> $arq_tmp
done

ytag="Memory bandwidth [MBytes/s]"
arq_saida="l3bandwidth_gauss.png"

gnuplot <<- EOF
reset
set terminal png size ${image_size} enhanced font 'Verdana,18' 
set style data linespoints
set style fill solid 2.00 border 0

set style line 1 lc rgb 'orange' lt 1 lw 2 pt 13 ps 1.0
set style line 2 lc rgb 'red' lt 1 lw 2 pt 7 ps 1.0
set style line 3 lc rgb 'grey' lt 1 lw 2 pt 13 ps 1.0
set style line 4 lc rgb 'green' lt 1 lw 2 pt 7 ps 1.0

set key font ",18" 
set key horizontal  
set key spacing 3 
set key samplen 3 
set key width 2

set title 'Medição de bandwidth L3 (Gauss-Seidel)' font ",24" 

set xrange ${x_range}

set xlabel "Tamanho de nx e ny" 
set xlabel font ",20" 
set ylabel "$ytag" 
set ylabel font ",20"

set output "$arq_saida"
plot '$arq_tmp' using 1:2 with linespoints ls 1 title "Trabalho 1", \
'$arq_tmp' using 1:3 with linespoints ls 2 title "Trabalho 2 (Otimizado)"
EOF
echo "Grafico plotado em $arq_saida"

# ---------- L3 BANDWIDTH (Residuo) ----------------------------------
echo "Calculando L3 BANDIWIDTH (Residuo)"

arq_tmp=l3bandwidth_residuo.out

touch $arq_tmp
cp /dev/null $arq_tmp


for i in ${elements[@]}; do
	val_old=$(likwid-perfctr -C 3 -g L3 -m -O ./t1/pdeSolver -nx $i -ny $i -i 10 | tail -n2 | head -n1 | cut -d ',' -f2)
	val_new=$(likwid-perfctr -C 3 -g L3 -m -O ./t2/pdeSolver -nx $i -ny $i -i 10 | tail -n2 | head -n1 | cut -d ',' -f2)
	echo "$i $val_old $val_new " >> $arq_tmp
done

ytag="Memory bandwidth [MBytes/s]"
arq_saida="l3bandwidth_residuo.png"

gnuplot <<- EOF
reset
set terminal png size ${image_size} enhanced font 'Verdana,18' 
set style data linespoints
set style fill solid 2.00 border 0

set style line 1 lc rgb 'orange' lt 1 lw 2 pt 13 ps 1.0
set style line 2 lc rgb 'red' lt 1 lw 2 pt 7 ps 1.0
set style line 3 lc rgb 'grey' lt 1 lw 2 pt 13 ps 1.0
set style line 4 lc rgb 'green' lt 1 lw 2 pt 7 ps 1.0

set key font ",18" 
set key horizontal  
set key spacing 3 
set key samplen 3 
set key width 2

set title 'Medição de bandwidth L3 (Resíduo)' font ",24" 

set xrange ${x_range}

set xlabel "Tamanho de nx e ny" 
set xlabel font ",20" 
set ylabel "$ytag" 
set ylabel font ",20"

set output "$arq_saida"
plot '$arq_tmp' using 1:2 with linespoints ls 1 title "Trabalho 1", \
'$arq_tmp' using 1:3 with linespoints ls 2 title "Trabalho 2 (Otimizado)"
EOF
echo "Grafico plotado em $arq_saida"

# ---------- FLOPS_DP e AVX (Gauss) ----------------------------------
echo "Calculando FLOPS DP e AVX (Gauss-Seidel)"

arq_tmp=flops_gauss.out

touch $arq_tmp
cp /dev/null $arq_tmp


for i in ${elements[@]}; do
	val_old_dp=$(likwid-perfctr -C 3 -g FLOPS_DP -m -O ./t1/pdeSolver -nx $i -ny $i -i 10 | tail -n31 | head -n1 | cut -d ',' -f2)
	val_new_dp=$(likwid-perfctr -C 3 -g FLOPS_DP -m -O ./t2/pdeSolver -nx $i -ny $i -i 10 | tail -n31 | head -n1 | cut -d ',' -f2)
	val_old_avx=$(likwid-perfctr -C 3 -g FLOPS_DP -m -O ./t1/pdeSolver -nx $i -ny $i -i 10 | tail -n30 | head -n1 | cut -d ',' -f2)
	val_new_avx=$(likwid-perfctr -C 3 -g FLOPS_DP -m -O ./t2/pdeSolver -nx $i -ny $i -i 10 | tail -n30 | head -n1 | cut -d ',' -f2)

	echo "$i $val_old_dp $val_new_dp $val_old_avx $val_new_avx" >> $arq_tmp
done

ytag="MFLOP/s"
arq_saida="mflops_gauss.png"

gnuplot <<- EOF
reset
set terminal png size ${image_size} enhanced font 'Verdana,18' 
set style data linespoints
set style fill solid 2.00 border 0

set style line 1 lc rgb 'orange' lt 1 lw 2 pt 13 ps 1.0
set style line 2 lc rgb 'red' lt 1 lw 2 pt 7 ps 1.0
set style line 3 lc rgb 'grey' lt 1 lw 2 pt 13 ps 1.0
set style line 4 lc rgb 'green' lt 1 lw 2 pt 7 ps 1.0

set key font ",18" 
set key horizontal  
set key spacing 3 
set key samplen 3 
set key width 2

set title 'Medição de FLOPS DP e FLOPS AVX (Gauss-Seidel)' font ",24" 

set xrange ${x_range}

set xlabel "Tamanho de nx e ny" 
set xlabel font ",20" 
set ylabel "$ytag" 
set ylabel font ",20"

set output "$arq_saida"
plot '$arq_tmp' using 1:2 with linespoints ls 1 title "Trabalho 1 - FLOPS DP", \
'$arq_tmp' using 1:3 with linespoints ls 2 title "Trabalho 2 - FLOPS DP", \
'$arq_tmp' using 1:4 with linespoints ls 3 title "Trabalho 1 - FLOPS AVX", \
'$arq_tmp' using 1:5 with linespoints ls 4 title "Trabalho 2 - FLOPS AVX"
EOF
echo "Grafico plotado em $arq_saida"

# ---------- FLOPS_DP e AVX (Residuo) ----------------------------------
echo "Calculando FLOPS DP e AVX (Residuo)"

arq_tmp=flops_residuo.out

touch $arq_tmp
cp /dev/null $arq_tmp


for i in ${elements[@]}; do
	val_old_dp=$(likwid-perfctr -C 3 -g FLOPS_DP -m -O ./t1/pdeSolver -nx $i -ny $i -i 10 | tail -n31 | head -n1 | cut -d ',' -f2)
	val_new_dp=$(likwid-perfctr -C 3 -g FLOPS_DP -m -O ./t2/pdeSolver -nx $i -ny $i -i 10 | tail -n31 | head -n1 | cut -d ',' -f2)
	val_old_avx=$(likwid-perfctr -C 3 -g FLOPS_DP -m -O ./t1/pdeSolver -nx $i -ny $i -i 10 | tail -n30 | head -n1 | cut -d ',' -f2)
	val_new_avx=$(likwid-perfctr -C 3 -g FLOPS_DP -m -O ./t2/pdeSolver -nx $i -ny $i -i 10 | tail -n30 | head -n1 | cut -d ',' -f2)

	echo "$i $val_old_dp $val_new_dp $val_old_avx $val_new_avx" >> $arq_tmp
done

ytag="MFLOP/s"
arq_saida="mflops_residuo.png"

gnuplot <<- EOF
reset
set terminal png size ${image_size} enhanced font 'Verdana,18' 
set style data linespoints
set style fill solid 2.00 border 0

set style line 1 lc rgb 'orange' lt 1 lw 2 pt 13 ps 1.0
set style line 2 lc rgb 'red' lt 1 lw 2 pt 7 ps 1.0
set style line 3 lc rgb 'grey' lt 1 lw 2 pt 13 ps 1.0
set style line 4 lc rgb 'green' lt 1 lw 2 pt 7 ps 1.0

set key font ",18" 
set key horizontal  
set key spacing 3 
set key samplen 3 
set key width 2

set title 'Medição de FLOPS DP e FLOPS AVX (Residuo)' font ",24" 

set xrange ${x_range}

set xlabel "Tamanho de nx e ny" 
set xlabel font ",20" 
set ylabel "$ytag" 
set ylabel font ",20"

set output "$arq_saida"
plot '$arq_tmp' using 1:2 with linespoints ls 1 title "Trabalho 1 - FLOPS DP", \
'$arq_tmp' using 1:3 with linespoints ls 2 title "Trabalho 2 - FLOPS DP", \
'$arq_tmp' using 1:4 with linespoints ls 3 title "Trabalho 1 - FLOPS AVX", \
'$arq_tmp' using 1:5 with linespoints ls 4 title "Trabalho 2 - FLOPS AVX"
EOF
echo "Grafico plotado em $arq_saida"


