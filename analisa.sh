echo "performance" > /sys/devices/system/cpu/cpufreq/policy3/scaling_governor

# get L2 miss ratio
likwid-perfctr -C 3 -O -g L2CACHE ./pdeSolver -nx 100 -ny 100 -i 10