

rule snpeffrun:
   input:
      completed = workp+"/07_TNhaplotyper/{tumor}.completed",
      completed2 = workp + "/08_cnvkitresult/{tumor}.completed",
   output:
      completed = workp + "/09_ECannotation/{tumor}.completed",
   params:
      snpEff = config["Application"]["snpEff"],
      annote = lambda wildcards: workp + "/09_ECannotation/%s_vs%s.TN.filter2.annotate.vcf"%(wildcards.tumor, config["tumor_vs_normal"][wildcards.tumor])
   shell:
       """
       # java -jar snpEff.jar download Triticum_aestivum # 全名为Triticum_aestivum的数据
库文件
                                                         # 下载后会在snpEff目录下生成一>个名为data的文件夹，
                                                         # 里面会有一个bin文件
       outputp=$(echo \"{outputp.completed}\" |grep -o \"/.*/\" )a
       inputp=$(echo \"{input.completed}\" |grep -o \"/.*/\" )
       inputfilter2=$(ls ${{inputp}}/{wildcards.tumor}_vs_*.TN.PASS.vcf |sort -u)
       if [ ! -d ${{outputp}} ] ;then mkdir -p ${{outputp}} ;fi
       {java} -Xmx8g -jar {params.snpEff} {genomeversion} \\
           ${{inputfilter2}} > {params.annote}
       echo "completed!!" >{output.completed}
       echo \"[## annotate ] date: $(date), parameter: {wildcards.tumor}, status: finished ...\" >>{status}
       """
