

def removetmp(wildcards):
   import os
   if os.path.exists("tmp"): os.system("rm -r tmp")
   return "/tmp"

rule Suplementary:
   input:
      expand( workp+"/" + resultnode2 + "/{sample}.final.completed", sample = samplelist )
   output:
      PPPP = workp+ "/12_finalresult/12_Supplementary.txt",
   params:
      Supple = config["Files"]["Supplementary"],
      MainModule = config["MainModule"],
      
   resources:
      tmpdir = removetmp,
   shell:
      """
      echo \"[## finalreault ] date: $(date), parameter: NULL, status: running ...\" >>{status}
      cp -rf {params.Supple} {output.PPPP}
      outputp=$(echo \"{output.PPPP}\" |grep -o \"/.*/\")
      if [ ! -d ${{outputp}}/QualityControl ]; then mkdir -p ${{outputp}}/QualityControl ;fi
      if [ ! -d ${{outputp}}/Copynumbervariant ]; then mkdir -p ${{outputp}}/Copynumbervariant ; fi
      if [ ! -d ${{outputp}}/SmallNucleotideVariation ]; then mkdir -p ${{outputp}}/SmallNucleotideVariation ;fi
      if [ -d {workp}/02_fastqc ] && [ $(ls {workp}/02_fastqc |grep "fastqc\.html$" |wc -l) -gt 0 ]
      then
         {config[Application][multiqc]} --outdir ${{outputp}}/QualityControl {workp}/02_fastqc
      fi
      if [[ {params.MainModule} == \"HaplotypeCaller\" ]]
      then
         echo \"Haplotypecaller模块的结果正在汇总中....\"
         ln -s {workp}/11_HaAnnotation/*.VQSR.funcotator.vcf* ${{outputp}}/SmallNucleotideVariation/

      elif [[ {params.MainModule} == \"Mutect2\" ]]
      then
         echo \"Mutect2模块的结果正在汇总中....\"
         ln -s {workp}/11_mutect2annotation/*_vs_*.Mt2.TN*.funcotator.vcf* \
            ${{outputp}}/SmallNucleotideVariation/
         if [ $(ls  {workp}/10_cnvcaller/ |grep \".*_vs_.*.called.seg.cnv\" |wc -l) -gt 0 ] 
         then
            ln -s {workp}/10_cnvcaller/*_vs_*.called.seg.cnv ${{outputp}}/Copynumbervariant/
         fi

      elif [[ {params.MainModule} == \"ECGEA\" ]]
      then
         echo \"ECGEA模块的结果正在汇总中....\"
         ln -s {workp}/07_TNhaplotyper/*_vs_*.filter2.TN.vcf* ${{outputp}}/SmallNucleotideVariation/
         ln -s {workp}/08_cnvkitresult/*_vs_* ${{outputp}}/Copynumbervariant/
      else
         echo \"这个模块脚本没办法识别呀，是不是名字错了~~\"
      fi
      echo \"[## finalresult ] date: $(date), parameter: NULL, status: finished ...\" >>{status}
      """
