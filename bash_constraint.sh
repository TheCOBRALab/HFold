input=$1
file="/Users/mateo2/Documents/Code/output/constraint/fasta/$input.txt";
# file2="/Users/mateo2/Documents/Code/output/structures/$input.txt";
TIMEFMT=$'%*E\t%M'
exec 5<$file
name="s"
while read line1 <&5; do

if ((i % 2 == 1))  

then
# echo $line1 > out1.txt
# echo $line2 >> out1.txt
./build/HFold -i $name -n 20 -d1 $line1
# { time ./build/Spark -d1 -r $line2 $line1 > "out1.txt"; } 2>"out.txt"
# echo "$i";
# { time ./build/Spark -P "params/rna_Turner04.par" -d1 -r $line2 $line1 > "out1.txt"; } 2>"out.txt"
# cat "out1.txt" >> "/Users/mateo2/Documents/Code/output/candidates/Turnerpk/$input.txt"
# cat "out.txt" >> "/Users/mateo2/Documents/Code/output/time/Spark/Turnerpk/$input.txt"
# cat "out1.txt" >> "/Users/mateo2/Documents/Code/output/proof/Spark/$input.txt"
else
name=$line1
fi
i=$((i+1));
done
echo "first is $i";
exec 5</dev/null
