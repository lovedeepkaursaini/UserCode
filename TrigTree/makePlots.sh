#!/bin/bash

# usage: 
# ./makePlots.sh 

if [[ $# != 0 ]]
    then
    echo "Usage: ./makePlots.sh"
    exit
fi

folderList="${CMSSW_BASE}/src/UserCode/TrigTree/SAVE"

root -b -q "${CMSSW_BASE}/src/UserCode/TrigTree/divide.C"

# generate webpages in each dir
for folder in ${folderList}
do
  
  cd $folder

  echo "Generating index.html for ${folder}"

  echo "<html><head><title>" >> index.html
  echo "Trigger Efficiency: ${folder}" >> index.html
  echo "</title></head>" >> index.html
  echo "<body>" >> index.html
  echo "<h1>Trigger Efficiency: ${folder} ${info}</h1>" >> index.html
  echo ""

  echo "Linking pictures for "
  pwd
  for picture in `ls *.gif`
  do
  # add stuff to filename and stick into index.html
    echo "<img src=\"${picture}\"> <br>" >> index.html
  done
  
  echo "<br><br>" >> index.html
  echo "</body></html>" >> index.html
  cd -
   
done

cd ..

exit
