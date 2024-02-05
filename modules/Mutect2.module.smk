
configfile: "null.yaml"

GATK = config["Application"]["GATK"]
java = config["Application"]["java"]
genomeversion = config["genomeversion"]
bundlehg = config["bundlehg"]
indelsvcf_1 = config["indelsvcf_1"]
indelsvcf_2 = config["indelsvcf_2"]
dbsnpvcf = config["dbsnpvcf"]
reference = config["reference"]
afonly = config["afonly"]
normallist = config["normallist"]
small_exac = config["small_exac"]
datatype_gatk = config["datatype_gatk"]
gatkinterval = config["gatkinterval"] # 捕获区间bed文件, 全基因组测序不需要该参数
workp = config["workp2"]
chrom = config["chrom"]
mt2_fir2_need = config["mt2_fir2_need"]
mt2_gmr_need = config["mt2_gmr_need"]
mt2_pon_need = config["mt2_pon_need"]
mt2_parallerbychrom = config["mt2_parallerbychrom"]


rule Mutect2part1ignore:
   input:
     BqsrBam = "",
   output:
     ignore = "",
   shell:
     """
     echo \"ignore mutect2 part1 by pon_need set false by user...\" >{output.ignore}
     """


rule Mutect2part1Bychrom:
   input:
      BqsrBam = workp+"/07_bqsrresult/{normal}.sorted.markdup.BQSR.bam",
   output:
      vcfpart = workp+"/08_Mutect2/{normal}.chromosome.vcf/{normal}.{chromno}.Mt2.normal.part.vcf",
   resources:
      tmpdir = "null",
      mem_mb = 1000,
   threads: 8
   params:
      staus = "null",
   shell:
      """
      # time {java} -jar {GATK} BedToIntervalList \\
      #    -I bed/{genomeversion}.exon.bed \\
      #    -O {{interval_mt2}} \\
      #    -SD $(echo "{reference}" |grep -o "^.*\." |sed "s/\.$//").dict
      
      chromno=\"{wildcards.chromno}\"
      if [[ {datatype_gatk} == "WES" ]] && [[ {mt2_parallerbychrom} == "bychrom" ]]
      then
         intervalparams=\"-L ${{chromno}} -L {gatkinterval} --interval-set-rule INTERSECTION \"
      elif [[ {datatype_gatk} != "WES" ]] && [[ {mt2_parallerbychrom} == "bychrom" ]]
      then
         intervalparams=\"-L ${{chromno}} \"
      elif [[ {datatype_gatk} == "WES" ]] && [[ {mt2_parallerbychrom} != "bychrom" ]]
      then
         intervalparams=\"-L {gatkinterval}\"
      else
         intervalparams=\"\"
      fi
         time {java} -Djava.io.tmpdir={resources.tmpdir} -jar {GATK} Mutect2 \\
            -R {reference} -I {input.BqsrBam} -O {output.vcfpart} \\
            ${{intervalparams}} --germline-resource {afonly} \\
            --max-mnp-distance 0 \\
            --disable-read-filter MateOnSameContigOrNoMappedMateReadFilter 
      """


rule Mutect2part1:
    input:
       expand(workp+"/null.{chromno}.null", chromno = chrom)
    output:
       expand( workp +"/{normal}.null", normal = normallist)
       # pon = workp+"null"
    threads: 8
    params:
       status = "null",
    shell:
       """
       outputp=$(echo \"{output[0]}\" |grep -o \"/.*/\")
       sample_map=${{outputp}}/wholeNormal.sample_map.txt
       for normalone in {normallist}
       do
          Mtvcfone="${{outputp}}/${{normalone}}.Mt2.normal.vcf"
          if [[ "{mt2_parallerbychrom}" == "bychrom" ]]
          then
             mergeline=\"\"
             for vcfone in $(echo \"{input}\" |sed -E "s/\s+/\n/g" |grep \"/${{normalone}}\.\")
             do 
                mergeline=\"${{mergeline}} -I ${{vcfone}}\"
             done
             time {java} -jar {GATK} MergeVcfs ${{mergeline}} -O  ${{Mtvcfone}}
          else
             ln -sf $(echo \"{input}\" |sed -E "s/\s+/\n/g" |grep \"/${{normalone}}\.\") ${{Mtvcfone}}
          fi
          echo \"${{normalone}}	${{Mtvcfone}}	1001\" >${{sample_map}}
       done
       pon_db=${{outputp}}/normal_pondb
       Lchrom=\"\" ######漏洞记得修复~~
       for chromno in {chrom}; do Lchrom=\"${{Lchrom}} -L ${{chromno}}\" ;done
       if [ -d ${{pon_db}} ] ;then rm -r ${{pon_db}} ;fi
       time {java} -jar {GATK} GenomicsDBImport \\
          -R {reference} --sample-name-map ${{sample_map}} \\
          --genomicsdb-workspace-path ${{pon_db}} \\
          ${{Lchrom}} --batch-size 50 --reader-threads 5
       time {java} -jar {GATK} CreateSomaticPanelOfNormals \\
          -R {reference} --germline-resource {afonly} \\
          -V gendb://${{pon_db}} \\
          -O ${{outputp}}/wholeNoraml.Mt2.pon.vcf
       fi
       echo \"[## Mutect2.part1 ] date: $(date), ${{pppp}}, status: finished ...\" >>{params.status}
       """


rule Mutect2part2:
   input:
      normal = lambda wildcards: "null",
      tumor = workp+"/07_bqsrresult/{tumor}.sorted.markdup.BQSR.bam",
   output:
      completed = workp+"/08_Mutect2/{tumor}.part2.1.completed"
   threads: 8
   params:
      status="null",
   shell:
      """
      outputp=$(echo \"{output.completed}\" |grep -o \"/.*/\")
      normal=$(echo \"{input.normal}\" |sed \"s/\/.*\///\" |cut -d \".\" -f1)
      tumor=$(echo \"{input.tumor}\" |sed \"s/\/.*\///\" |cut -d \".\" -f1)
      tumor_vs_normal=${{tumor}}_vs_${{normal}}
      inputbam={input.tumor}
      prefix=${{outputp}}/${{tumor}}
      Lparams=\"\"
      if [[ "{datatype_gatk}" == "WES" ]]
      then 
         Lparams=\"-L {gatkinterval} --interval-set-rule INTERSECTION \"
      fi
      if [ -f {small_exac} ]
      then
         for couple in \"tumor\" \"normal\"
         do
            if [[ "${{couple}}" == \"normal\" ]]
            then
               inputbam={input.normal}
               prefix=${{outputp}}/${{normal}}
            fi

            time {java} -jar {GATK} GetPileupSummaries \\
               -R {reference} -I ${{inputbam}} \\
               -V {small_exac} -L {small_exac} \\
               ${{Lparams}} -O ${{prefix}}.getpileupsum.table
         done
         time {java} -jar {GATK} CalculateContamination \\
            -I ${{outputp}}/${{tumor}}.getpileupsum.table \\
            -matched ${{outputp}}/${{normal}}.getpileupsum.table \\
            -O ${{outputp}}/${{tumor_vs_normal}}.contamination.TN.table \\
            --tumor-segmentation ${{outputp}}/${{tumor_vs_normal}}.segmentation.TN.table
      fi
      echo \"completed!!\" >{output.completed}
      """


rule Mutect2part2bychrom:
   input:
      Mtnormal = "null",
      tumor = workp+"/07_bqsrresult/{tumor}.sorted.markdup.BQSR.bam",
      completed = workp+"/08_Mutect2/{tumor}.part2.1.completed",
   output:
      completed = workp+"/08_Mutect2/{tumor}.Mutect2.chromosome/{tumor}.{chromno}.part2.2.completed",
   resources:
      tmpdir = "null",
      mem_mb = 5000,
   threads: 16
   params:
      status = "null",
   shell:
      """
      normal=$( echo \"{input.Mtnormal}\" | sed \"s/\/.*\///\" |cut -d \".\" -f1 )
      inputp=$(echo \"{input.Mtnormal}\" |grep -o "^.*/" )
      outputp=$(echo \"{output.completed}\" |grep -o \"/.*/\")
      tumor={wildcards.tumor}
      chromno=\"{wildcards.chromno}\"
      tumor_vs_normal=${{tumor}}_vs_${{normal}}
      intervalparams=\"-L ${{chromno}}\"
      if [[ {datatype_gatk} == "WES" ]]
      then
         if [[ {mt2_parallerbychrom} == "bychrom" ]]
         then
            pppp=\"--interval-set-rule INTERSECTION\"
            intervalparams=\"${{intervalparams}} -L {gatkinterval} ${{pppp}} \"
         else
            intervalparams=\"-L {gatkinterval}\"
         fi
      else
         if [[ {mt2_parallerbychrom} != "bychrom" ]] ;then intervalparams=\"\" ;fi
      fi
      afonlyparams=\"\"
      if {mt2_gmr_need}; then afonlyparams=\"--germline-resource {afonly}\" ;fi
      ponparams=\"\"
      if {mt2_pon_need}; then ponparams=\"--panel-of-normals ${{inputp}}/wholeNoraml.Mt2.pon.vcf \" ;fi
      f1r2params=\"\"
      if {mt2_fir2_need} 
      then 
         f1r2params=\" --f1r2-tar-gz ${{outputp}}/${{tumor_vs_normal}}.${{chromno}}.f1r2.TN.part.tar.gz \"
      fi
      normalbam=$(echo \"{input.tumor}\" |sed \"s/\/{wildcards.tumor}\./\/${{normal}}\./\")
      time {java} -Djava.io.tmpdir={resources.tmpdir} -jar {GATK} Mutect2 \\
         -R {reference} -I ${{normalbam}} -I {input.tumor} \\
         --tumor-sample ${{tumor}} --normal-sample ${{normal}} \\
         ${{intervalparams}} ${{ponparams}} ${{afonlyparams}} ${{f1r2params}} \\
         `# --disable-read-filter MateOnSameContigOrNoMappedMateReadFilter` \\
         `# --af-of-alleles-not-in-resource 5.0E-8 ` \\
         `# --output-mode EMIT_VARIANTS_ONLY` \\
         --bam-output ${{outputp}}/${{tumor_vs_normal}}.${{chromno}}.bamout.TN.part.bam \\
         -O ${{outputp}}/${{tumor_vs_normal}}.${{chromno}}.Mt2.TN.part.vcf
      echo \"completed\" >{output.completed}
      """

rule Mutect2part3:
   input:
      expand(workp+"/08_Mutect2/{{tumor}}.Mutect2.chromosome/null.completed", chromno = chrom )
   output:
      completed = workp+"/08_Mutect2/{tumor}.part3.completed"
   threads: 8
   params:
      status = "null",
   shell:
      """
      inputp=$(echo \"{input[0]}\" | grep -o \"/.*/\")
      outputp=$(echo \"{output.completed}\"|grep -o \"/.*/\")
      tumor_vs_normal=$(ls ${{inputp}}|grep -o \"{wildcards.tumor}_vs_[^.]*\" |sort -u)
      prefix=${{outputp}}/${{tumor_vs_normal}}
      prefix2=${{inputp}}/${{tumor_vs_normal}}
      artifactparams=\"\"
      if {mt2_fir2_need}
      then
         if [[ {mt2_parallerbychrom} == "bychrom" ]]
         then
            pppp=$(ls ${{inputp}}*_vs_*.f1r2.TN.part.tar.gz )
            all_f1r2_input=`for f1r2file in ${{pppp}}; do printf -- \"-I ${{f1r2file}} \"; done`
            time {java} -jar {GATK} LearnReadOrientationModel \\
                ${{all_f1r2_input}} -O ${{prefix}}.artifact.TN.tar.gz
         else
            ln -sf ${{prefix2}}.chromall.f1r2.TN.part.tar.gz ${{prefix}}.artifact.TN.tar.gz
         fi
         artifactparams=\" --ob-priors ${{prefix}}.artifact.TN.tar.gz \"
      fi
      if [[ {mt2_parallerbychrom} == "bychrom" ]]
      then
         mergeline=\"\"
         mergeline2=\"\"
         stateExist=true
         for vcffile in $(ls ${{prefix2}}.chr*.Mt2.TN.part.vcf )
         do
            mergeline=\"${{mergeline}} -I ${{vcffile}}\"
            mergeline2="${{mergeline2}} --stats ${{vcffile}}.stats"
            if ${{stateExist}} && [ ! -f ${{vcffile}}.stats ] ;then stateExist=false ;fi
         done
         time {java} -jar {GATK} MergeVcfs ${{mergeline}} -O ${{prefix}}.Mt2.TN.vcf
         if ${{stateExist}}
         then
            time {java} -jar {GATK} MergeMutectStats ${{mergeline2}} -O ${{prefix}}.Mt2.TN.vcf.stats
         fi
      else
         ln -sf ${{prefix2}}.chromall.Mt2.TN.part.vcf ${{prefix}}.Mt2.TN.vcf
         ln -sf ${{prefix2}}.chromall.Mt2.TN.part.vcf.stats ${{prefix}}.Mt2.TN.vcf.stats
      fi
      statsparams=\"\"
      contamination_table_params=\"\"
      tumor_segmentation_params=\"\"
      if [ -f ${{prefix}}.contamination.TN.table ]
      then
         contamination_table_params=\"--contamination-table ${{prefix}}.contamination.TN.table\"
         tumor_segmentation_params=\"--tumor-segmentation ${{prefix}}.segmentation.TN.table \"
      fi
      if [ -f ${{prefix}}.Mt2.TN.vcf.stats ]; then statsparams=\"-stats ${{prefix}}.Mt2.TN.vcf.stats\" ;fi
      time {java} -jar {GATK} FilterMutectCalls \\
         -R {reference} -V ${{prefix}}.Mt2.TN.vcf \\
         ${{contamination_table_params}} ${{tumor_segmentation_params}} \\
         ${{artifactparams}} ${{statsparams}} \\
         -O ${{prefix}}.filtered0.TN.vcf
       
      # time {java} -jar {GATK} BwaMemIndexImageCreator \\
      #    -I {reference} \\
      #    -O $(echo "{reference}" |grep -o "^.*\." |sed "s/\.$//").bwaindex.img 

      echo \"completed!!\" >{output.completed}
      pppp=\"parameter: {wildcards.tumor}\"
      echo \"[## Mutect2.part3 ] date: $(date), ${{pppp}}, status: finished ...\" >>{params.status}
      """

