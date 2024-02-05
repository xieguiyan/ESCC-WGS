## BQSR snakemake script

rule BQSRrun:
   input:
      SortedMarkdupBam = workp+"/05_markduplicate/{sample}.BWAMEM.sorted.markdup.bam",
      ignore = workp+"/06_realigntarget/{sample}.ignore",
   output:
      table = workp+"/07_bqsrresult/{sample}.recal_data.table",
      BqsrBam = workp+"/07_bqsrresult/{sample}.sorted.markdup.BQSR.bam",
   params:
      knownsitefile = dowithgatkparams.obtainBQSRknownsite(config),
   shell:
      """
      inputbam=$(echo \"{input.ignore}\"|grep -o \"/.*\.\")BWAMEM.sorted.markdup.realign.bam
      if [ ! -f ${{inputbam}} ];then inputbam={input.SortedMarkdupBam};fi
      Lparams=\"\"
      if [[ {datatype_gatk} == \"WES\" ]] ;then Lparams=\"-L {gatkinterval}\" ;fi
      knownsiteparams=\"\"
      for vcffile in {params.knownsitefile}
      do
         if [ -f ${{vcffile}} ];then knownsiteparams=\"${{knownsiteparams}} --known-sites ${{vcffile}} \" ;fi
      done
      time {java} -jar {GATK} BaseRecalibrator -R {reference} \\
         -I ${{inputbam}} ${{Lparams}} \\
         ${{knownsiteparams}} \\
         -O {output.table}
      {java} -jar {GATK} ApplyBQSR --bqsr-recal-file {output.table} -R {reference} \\
         -I ${{inputbam}} -O {output.BqsrBam}
      echo \"[## BQSR ] date: $(date), parameter: {wildcards.sample}, status: finished ...\" >>{status}
      """
