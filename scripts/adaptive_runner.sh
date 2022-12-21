#!usr/bin/env bash

version="v5"

qbins="b10"
errmode="ABS"
prediction="lorenzo"

appidx=5
file_name="velocity_z"
fname=${file_name}
#fname=""
prefix1="/home/alpoulo/datasets"
prefix2="/zfs/fthpc/alpoulo/datasets"

apps=("HACC" "CESM2D" "exaalt-copper" "ISABEL" "CESM3D" "NYX" "MIRANDA")
dim=("-1 1073726487" "-2 3600 1800" "-2 5243 3137" "-3 500 500 100" "-3 3600 1800 26" "-3 512 512 512" "-3 384 384 256")
exts=(".f32" "_1_1800_3600.f32" ".f32.dat" ".bin.f32" "_1_26_1800_3600.f32" ".dat" ".d64")
dtypes=("-f" "-f" "-f" "-f" "-f" "-f" "-d")
prefixes=("${prefix2}" "${prefix1}" "${prefix1}" "${prefix1}" "${prefix2}" "${prefix2}" "${prefix1}")

app=${apps[$appidx]}
dims=${dim[$appidx]}
dtype=${dtypes[$appidx]}
ext=${exts[$appidx]}
prefix=${prefixes[$appidx]}
infile="${prefix}/${app}/${file_name}${ext}"

adaptivefile="${fname}_adaptive"
outfile="${adaptivefile}_${errmode}_${prediction}_${qbins}_${version}.out"
baselinefile="${fname}_adaptive"

scratch="/scratch1/alpoulo/adaptive_quantization/output/${app}/tests/${errmode}/${prediction}/${qbins}"
mkdir -p $scratch

entropydir="${scratch}/entropy/${version}"
mkdir -p $entropydir

histogramdir="${scratch}/histogram/${version}"
mkdir -p $histogramdir

cmpdir="${scratch}/cmp/${version}"
mkdir -p $cmpdir


for j in {1..1}
do
    for i in {3..3}
    do
        let adapt_eb=2**${i}
        quantentropyfile="${entropydir}/${baselinefile}_q${i}_1e-${j}_quant.out"
        adaptentropyfile="${entropydir}/${baselinefile}_q${i}_1e-${j}_adapt.out"
        quanthistfile="${histogramdir}/${baselinefile}_q${i}_1e-${j}_quanthist.out"
        adapthistfile="${histogramdir}/${baselinefile}_q${i}_1e-${j}_adapthist.out"
        export QUANTENTROPY=${quantentropyfile}
        export ADAPTENTROPY=${adaptentropyfile}
        export QUANTHISTOGRAM=${quanthistfile}
        export ADAPTHISTOGRAM=${adapthistfile}
        cmpfile="${cmpdir}/${adaptivefile}_${version}_q${i}_1e-${j}.dat.sz.out"
        ./sz3_$version $dtype -i $infile -o $cmpfile $dims -c sz3.config -M $errmode ${adapt_eb}e-$j -a -q $i >> $outfile
        printProperty -f $quantentropyfile 3 >> $outfile
        printProperty -f $adaptentropyfile 3 >> $outfile
        #rm ${cmpfile}
        rm ${quantentropyfile}
        rm ${adaptentropyfile}
        printf "\n" >> $outfile
    done
    printf "\n" >> $outfile
done

#app="HACC"
#dims="-1 1073726487"
#fname="vx"
#infile="/zfs/fthpc/alpoulo/datasets/${app}/${fname}.f32"

#app="CESM2D"
#dims="-2 3600 1800"
#fname="FSDTOA"
#infile="/home/alpoulo/datasets/${app}/${fname}_1_1800_3600.f32"

#app="exaalt-copper"
#dims="-2 5243 3137"
#infname="dataset1-5423x3137.z.f32.dat"
#fname="dataset1-z"
#infile="/home/alpoulo/datasets/${app}/${infname}"

#app="ISABEL"
#dims="-3 500 500 100"
#fname="TCf48"
#infile="/home/alpoulo/datasets/${app}/${fname}.bin.f32"

#app="CESM3D"
#dims="-3 3600 1800 26"
#fname="OMEGA"
#infile="/zfs/fthpc/alpoulo/datasets/${app}/${fname}_1_26_1800_3600.f32"

#app="NYX"
#dims="-3 512 512 512"
#fname="dark_matter_density"
#infile="/zfs/fthpc/alpoulo/datasets/${app}/${fname}.dat"

#dtype="-d"
#app="Miranda"
#dims="-3 384 384 256"
#fname="pressure"
#infile="/home/alpoulo/datasets/${app}/${fname}.d64"


