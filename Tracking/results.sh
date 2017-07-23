
cp /dev/null fitresults.txt 
mkdir -p results

HISTO_LIST=(deltaM_5 deltaM_6)
MACRO_LIST=(fit2bDATA.C fit4bDATA.C fit2bMC.C fit4bMC.C)

for file in ${HISTO_LIST[*]}; do

for macro in ${MACRO_LIST[*]}; do

sed -i -e "s/input_file/${file}/" $macro
sed -i -e "s/fitresults/${file}/" $macro
root -l $macro 
sed -i -e "s/${file}/input_file/" $macro
sed -i -e "s/${file}/fitresults/" $macro

done

done
python results.py


