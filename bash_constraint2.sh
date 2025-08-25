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
./build/HFold -n 200 $line1 > out.txt
echo $name >> "/Users/mateo2/Documents/Code/output/constraint/structures/$input.txt"
cat "out.txt" >> "/Users/mateo2/Documents/Code/output/constraint/structures/$input.txt"

else
name=$line1
fi
i=$((i+1));
done
echo "first is $i";
exec 5</dev/null