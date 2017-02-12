#edit base dir-----path to folder with lc files

declare -a sdssID
let i=0
while IFS=$'\n' read -r line_data; do
    # Parse “${line_data}” to produce content 
    sdssID[i]="${line_data}" # Populate array.
    ((++i))
done < sdssList.dat
echo i


let x=0
while IFS=$'\n' read -r object; do
python combineLCPlotting.py -id ${object} -sdssid ${sdssID[x]}


echo "ran routine"
((++x))
echo x ${object} ${sdssID[x]}
done < k2List.dat
