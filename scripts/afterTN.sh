#!/bin/bash
prefix=$(echo "${snakemake_input[0]}" |grep -o "/.*/")
for couple in ${snakemake_input} "UU"
do
   if [[ ${couple} != "UU" ]]
   then
      tumor=$(echo "${couple}" |grep -o "[^/]*$" |cut -d "." -f1)
      couple=($(ls ${prefix}|grep -o "${tumor}_vs_[^.]*" |sort -u |sed "s/_vs_/ /" ))
      echo "completed!!" >${prefix}/${couple[0]}.final.completed
      [ ! -f ${prefix}/${couple[1]}.final.completed  ] && ln -s ${prefix}/${couple[0]}.final.completed ${prefix}/${couple[1]}.final.completed
   else
      echo "随便凑了个UU，因为如果样本少，就循环bu下去了~~~"
   fi
done

