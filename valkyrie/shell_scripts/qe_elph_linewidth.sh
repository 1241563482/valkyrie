#!/bin/bash
# merge ph eigv and a2f into x-y-size format

eigv="ph.freq.gp"
gam="elph.gamma.7"
elph="linewidth_phonon.dat"

# nbnd=  x, nks= y 
head_line=$(sed -n '1'p $gam)
nbnd=$(echo $head_line | sed "s/.*nbnd=\s*//;s/,.*//")
nks=$(echo $head_line | sed "s/.*nks=\s*//;s/\s*\///")

tot=$(wc -l $gam | awk '{print $1}')
num=$(( ( $tot - 1 ) / $nks ))

a2f_tmp=a2f.tmp
 > $a2f_tmp

# arrange a2F data
sed -n '2,$'p $gam | \
for ((i=1; i<=$tot; ++i))
do
	read line
	echo -n "$line " >> $a2f_tmp
	if [[ $(( $i % $num )) = 0 ]]; then
		echo ""  >> $a2f_tmp
	fi
done

> $elph
for (( i=0; i<$nbnd; ++i ))
do
	ph_col=$(($i + 2))
	a2f_col=$(($i + 4))

	awk '{print $1 "\t" $'$ph_col'}' $eigv > ph_eigv.$i.tmp
	awk '{print $'$a2f_col' }' $a2f_tmp > a2f.$i.tmp

	#paste ph_eigv.$i.tmp a2f.$i.tmp > $elph.$i.dat
	paste ph_eigv.$i.tmp a2f.$i.tmp >> $elph
	rm ph_eigv.$i.tmp a2f.$i.tmp
done
