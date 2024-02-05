
configfile: "null.yaml"
GATK = config["Application"]["GATK"]
java = config["Application"]["java"]
genomeversion = config["genomeversion"]
bundlehg = config["bundlehg"]
indelsvcf_1 = config["indelsvcf_1"]
indelsvcf_2 = config["indelsvcf_2"]
dbsnpvcf = config["dbsnpvcf"]
reference = config["reference"]
gatkinterval = config["gatkinterval"]
workp = config["workp2"]
tumorlist = config["tumorlist"]
normallist = config["normallist"]
datatype_gatk = config["datatype_gatk"]
scnv_Allelic_need = config["scnv_Allelic_need"]
chrom = config["chrom"]
workp = config["workp2"]


rule gatkCNVcallerignore:
   input:
      normal = "",
      tumor = "",
   output:
      ignore = workp+"/10_cnvcaller/{tumor}.part2.2.ignore"
   shell:
      """
      echo \"因为设置了sCNVcaller false 这部分就不跑了~~~\" >{output.ignore}
      """

rule gatkCNVcallerpairpart1:
   input:
      "null.bam",
   output:
      "part1.completed",
   resources:
      mem_mb=10000 # 就在这里报错了内存溢出，堆栈错误的情况，所以就特意分配了内存空间
   params:
      status = "null"
   shell:
      """
      normalline=$(echo \"{normallist}\" |sed -E \"s/\s+/;/g\")
      echo \"${{normalline}}\"
      outputp=$(echo \"{output[0]}\"|grep -o \"/.*/\")
      inputp=$(echo \"{input[0]}\"|grep -o \"/.*/\")
      if [ ! -d ${{outputp}} ] ;then mkdir -p ${{outputp}} ;fi
      interval_list=${{outputp}}/version.{genomeversion}.preprocessed.interval_list
      # if [[ {datatype_gatk} == \"WES\" ]] && [ -f {gatkinterval} ]
      echo "{gatkinterval}"
      if [ -f {gatkinterval} ]
      then
         time {java} -jar {GATK} PreprocessIntervals \\
            -L {gatkinterval} --reference {reference} \\
            --output ${{interval_list}} --padding 250 --bin-length 0 \\
            --interval-merging-rule OVERLAPPING_ONLY
      fi

      # time {java} -jar {GATK} AnnotateIntervals \\
      #   -L Covered.bed \\
      #   --reference {reference} \\
      #   --output ${{outputp}}/annotated.tsv \\
      #   --interval-merging-rule OVERLAPPING_ONLY

      # time {java} -jar {GATK} FilterIntervals \\
      #   -L ${{outputp}}/preprocessed.interval_list \\
      #   --annotated-intervals ${{outputp}}/annotated.tsv \\
      #   -O ${{outputp}}/filtered.interval_list \\
      #   --interval-merging-rule OVERLAPPING_ONLY
      # interval_list=\"${{outputp}}/$filtered.interval_list\"
      if [ -f ${{interval_list}} ]
      then
         for tumorone in {tumorlist}
         do
            tumor=${{tumorone}}
            normal=$(echo "${{normalline}}" | grep -o -E \"^[^;]+\")
            normalline=$(echo "${{normalline}}" | sed -E \"s/^[^;]+;//\")
            sname=${{tumor}}
            for couple in \"tumor\" \"normal\"
            do
               if [[ ${{couple}} == \"normal\"  ]]; then sname=${{normal}} ;fi
               time {java} -Xmx5000m -jar {GATK} CollectReadCounts \\
                  -L ${{interval_list}} \\
                  --reference {reference} \\
                  --input ${{inputp}}/${{sname}}.sorted.markdup.BQSR.bam \\
                  --format HDF5 --interval-merging-rule OVERLAPPING_ONLY \\
                  --output ${{outputp}}/${{sname}}.count.HDF5
            done
         done
         normalhdf5param=\"\"
         for normalone in {normallist}
         do 
            normalhdf5param=\"${{normalhdf5param}} --input ${{outputp}}/${{normalone}}.count.HDF5 \"
         done
         time {java} -Xmx5000m -jar {GATK} CreateReadCountPanelOfNormals \\
            ${{normalhdf5param}} \\
            --minimum-interval-median-percentile 5.0  \\
            --output ${{outputp}}/wholenormal.cnvpon.hdf5
      fi
      normalline=$(echo \"{normallist}\" |sed -E \"s/\s+/;/g\")
      for tumorone in {tumorlist}
      do
         tumor=${{tumorone}}
         normal=$(echo "${{normalline}}" | grep -o -E \"^[^;]+\")
         normalline=$(echo "${{normalline}}" | sed -E \"s/^[^;]+;//\")
         prefix=${{outputp}}/${{tumor}}_vs_${{normal}}
         if [ -f ${{interval_list}} ]
         then
            time {java} -Xmx5000m -jar {GATK} DenoiseReadCounts \\
               --input ${{outputp}}/${{tumor}}.count.HDF5 \\
               --count-panel-of-normals ${{outputp}}/wholenormal.cnvpon.hdf5 \\
               --standardized-copy-ratios ${{prefix}}.standardizedCR.tsv \\
               --denoised-copy-ratios ${{prefix}}.denoisedCR.tsv
         fi
         echo \"completed!!\" >${{outputp}}/${{tumor}}.part1.completed
      done
      """


rule gatkCNVcallerpairpart2bychrom:
   input:
      tumor = "",
      completed = "",
   output:
      completed = "",
   resources:
      mem_mb = 10000,
      tmpdir = "null",
   params:
      normal = lambda wildcards: "null",
      status = config["status"],
   shell:
      """
      tumor={wildcards.tumor}
      sname={params.normal}
      inputbam=$(echo \"{input.tumor}\" | sed \"s/\/${{tumor}}\./\/{params.normal}\./\" )
      outputp=$(echo \"{output.completed}\" |grep -o \"/.*/\")
      chromno={wildcards.chromno}
      interval_list=${{outputp%/*/}}/version.{genomeversion}.preprocessed.interval_list
      intervalpara=\"\"
      if [ -f ${{interval_list}} ] && {scnv_Allelic_need}
      then
         if [[ ${{chromno}} != \"chromall\" ]]
         then
            intervalpara=\"-L ${{chromno}} \"
            intervalpara=\"${{intervalpara}} -L ${{interval_list}} --interval-set-rule INTERSECTION\"
         fi
         for couple in \"normal\" \"tumor\"
         do
            if [[ ${{couple}} == \"tumor\" ]]; then inputbam={input.tumor} ;sname=${{tumor}} ;fi
            time {java} -Xmx10G -Djava.io.tmpdir={resources.tmpdir} -jar {GATK} CollectAllelicCounts \\
               ${{intervalpara}} --input ${{inputbam}} \\
               --reference {reference} --minimum-base-quality 20 \\
               --output ${{outputp}}/${{sname}}.${{chromno}}.part.allelic_counts
         done
      fi
      echo \"completed!!\" >{output.completed}
      """


rule gatkCNVcallerpairpart2merge:
   input:
      "null.cnv",
   output:
      completed = workp+"/10_cnvcaller/{tumor}.part2.2.completed",
   params:
      tumor_vs_normal = lambda wildcards: "null_vs_"+ config["tumor_vs_normal"][wildcards.tumor],
      GATKsCNVchrom = "chromall",
   resources:
      mem_mb = 20000,
   shell:
      """
      outputp=$(echo \"{output.completed}\" |grep -o \"/.*/\")
      outputcnv=${{outputp}}/{params.tumor_vs_normal}.called.seg.cnv
      inputp=$(echo \"{input[0]}\" |grep -o \"/.*/\")
      prefix=${{inputp}}/{params.tumor_vs_normal}
      inputp2=${{inputp%/*/}}
      prefix2=${{inputp2}}/{params.tumor_vs_normal}
      normal=$(echo \"{params.tumor_vs_normal}\" |sed \"s/^.*_//\")
      tumor_alleclic=${{inputp2}}/{wildcards.tumor}.chromall.allelic_counts
      normal_alleclic=${{inputp2}}/${{normal}}.chromall.allelic_counts
      interval_list=${{outputp}}/version.{genomeversion}.preprocessed.interval_list
      if [ -f ${{interval_list}} ]
      then
         allelicparams=\"\"
         if {scnv_Allelic_need} # ModelSegments的时候有两个策略，是否利用CollectAllelicCounts的结果
         then
            if [[ {params.GATKsCNVchrom} != \"chromall\" ]]
            then
               # CNV不知道可不可以按照染色体做，所以谨慎点还是不分染色体做好~
               for couple in \"${{normal}}\" \"{wildcards.tumor}\"
               do
                  cat ${{inputp}}/${{couple}}.chr1.part.allelic_counts \
                     >${{inputp2}}/${{couple}}.chromall.allelic_counts
                  pppp=$(ls ${{inputp}}/${{couple}}.chr*.part.allelic_counts |grep -E \"chr[12]?[0-9XYM]\.\")
                  for allelic in $(echo \"${{pppp}}\" |grep -v \"chr1\.\")
                  do
                     sed \"1,/^CONTIG/d\" ${{allelic}} >>${{inputp2}}/${{couple}}.chromall.allelic_counts
                  done
              done
            else
               ln -sf ${{inputp}}/{wildcards.tumor}.chromall.part.allelic_counts ${{tumor_alleclic}}
               ln -sf ${{inputp}}/${{normal}}.chromall.part.allelic_counts ${{normal_alleclic}}
            fi
            allelicparams=\" --allelic-counts ${{tumor_alleclic}} --normal-allelic-counts ${{normal_alleclic}} \"
         fi
         time {java} -Xmx20G -jar -Djava.io.tmpdir={resources.tmpdir} {GATK} ModelSegments \\
            --denoised-copy-ratios ${{prefix2}}.denoisedCR.tsv \\
            ${{allelicparams}} \\
            --output ${{outputp}} --output-prefix {params.tumor_vs_normal}_chromall
         time {java} -Xmx5000m -jar {GATK} CallCopyRatioSegments \\
            -I ${{prefix2}}_chromall.cr.seg \\
            -O ${{outputcnv}}
      else
         echo "interval file missing, which is required by CNV caller,then skipped this part... " > ${{outputcnv}}
      fi
      echo \"completed\" >{output.completed}
      pppp=\"parameter: {wildcards.tumor}\"
      echo \"[## sCNVcaller ] date: $(date), ${{pppp}}, status: finished ...\" >>{config[status]}
      """
