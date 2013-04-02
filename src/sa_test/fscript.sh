
rm -rf 0 1 
mkdir 0
mv feff-0.inp ./0/feff.inp
cp chi.chi3 ./0
cd 0
feff6 > /dev/null
bash ../process_ifeffit.sh
ifeffit -q process.iff
cd ..

mkdir 1
mv feff-1.inp ./1/feff.inp
cp chi.chi3 ./1
cd 1
feff6 > /dev/null
bash ../process_ifeffit.sh
ifeffit -q process.iff
cd ..

folder0=($(grep -v "#" ./0/my_chi.chi3))
folder1=($(grep -v "#" ./1/my_chi.chi3))

rm -rf dod.dat

tLen=${#folder0[@]}

for (( i=0; i<${tLen}; i=i+2 ));
do
	total=$( echo "${folder0[$i+1]} + ${folder1[$i+1]}" | awk -F "+" '{print $1 + $2 }')
	avg=$( echo $total / 2 | awk -F "/" '{print $1 / $2}')
	$(echo "${folder0[$i]} $avg" >> dod.dat)
done

rm -rf 0 1 
