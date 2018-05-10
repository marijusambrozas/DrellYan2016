#!bin/bash

qsub -V v20180128_M15to64_DiRapBin00.sh
qsub -V v20180128_M15to64_DiRapBin01.sh
qsub -V v20180128_M15to64_DiRapBin02.sh
qsub -V v20180128_M15to64_DiRapBin03.sh
qsub -V v20180128_M15to64_DiRapBin04.sh
qsub -V v20180128_M64to106_DiRapBin00.sh
qsub -V v20180128_M64to106_DiRapBin01.sh
qsub -V v20180128_M64to106_DiRapBin02.sh
qsub -V v20180128_M106to3000_DiRapBin00.sh
qsub -V v20180128_M106to3000_DiRapBin01.sh

echo "finished"
