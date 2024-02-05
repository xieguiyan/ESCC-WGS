#!/usr/bin/python3
import re
import os
import sys
import argparse
import time


def generatesbatchscript(workp, snakemakepath):
   print("###  生成sbatch投递脚本....")
   with open("%s/run.sh"%workp, "w+") as sbatchthis:
      sbatchthis.write("#!/bin/bash\n")
      sbatchthis.write("#SBATCH -p cn\n")
      sbatchthis.write("#SBATCH -J prowork\n")
      sbatchthis.write("#SBATCH -N 1\n")
      sbatchthis.write("#SBATCH -n 48\n")
      sbatchthis.write("#SBATCH -o %s/work3.o\n"%workp )
      sbatchthis.write("#SBATCH -e %s/work3.e\n"%workp )
      sbatchthis.write("# %s -s %s/work.smk --cores 48\n"%(snakemakepath, workp))
      sbatchthis.write("echo \"hello world!!\" \n")
      sbatchthis.close()


def check_sample_information(sample_info):
   print("正在确认样品信息的表格是否正确....")
   time.sleep(3)
   if os.path.exists(sample_info):
      try:
         with open(sample_info, "r+") as sampletable:
            header = sampletable.readline().strip().split("\t")
            if "tumor" not in header or "normal" not in header:
               print("当前的WGS流程，体细胞突变检测是tumor.vs.normal模式, 所以tunmor和normal两列信息都要有~~")
               if "tumor" not in header: print("找不到列名为tumor的那列")
               else: print("找不到列名为normal\n")
               exit(-1)
            sampletable.close()
      except:
         print("表格打不开呀~")
         exit(-1)
   else:
      print("没有找到这张表呢~~")
      exit(-1)
   print("表格格式没有问题")
   print("  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")

def checkpipemodules(abspath, modulelist ):
  print("核对下pipeline")
  for mdone in modulelist:
    if not os.path.exists( abspath + "/" + mdone ): 
       print("ERROR:: Autoconfig.py要和WGS的流程下'%s'在同个目录下"%mdone )
       exit(-1)
  print("  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")

def checksoftwareanddata(configfile, sd): # sd, software and database absolute path
   print("正在核对WGS流程需要的软件和工具::::::")   
   sddict2 = {
      "HumanGenomeFasta":"", "bundle":"", "GATK":"",
      "fastqc":"", "fastp":"", "bwa":"", "java":"",
      "samtools":"", "picard":"", "python":"", 
      "Rscript":"", "snakemake":"",
   }

   sddinfo = {
      "HumanGenomeFasta":"bwa比对的参考基因组库路径", 
      "bundle":"gatk依赖的数据资源路径(缺少请从官网下载)",
      "GATK":"gatk工具的路径(*/gatk*,jar)",
      "fastqc":"fastqc路径", 
      "fastp":"fastp路径", 
      "bwa":"bwa路径", 
      "java":"java路径",
      "samtools":"samtools路径", 
      "picard":"picard路径",
      "python":"python路径", 
      "Rscript":"Rscript路径，gatk作图会用到",
      "snakemake":"流程搭建工具",
   }
   if os.path.exists(sd): cfg = sd
   else: cfg = configfile
   with open(cfg, "r+") as configthis:
      for line in configthis.readlines():
         if "#" in line: line = line.split("#")[0]
         line = line.strip()
         if line == "": continue
         line = re.sub("^\s+", "", line)
         sdone = line.split(":")[0]
         if sdone in list( sddict2.keys() ):
            print(f"    核对{sdone}的路径~")
            time.sleep(1)
            sdpath = re.sub("\s+", "", line.split(":")[1]).strip()
            print(sdone, ",,, PATH: %s,"%sdpath)
            if sdpath != "" and os.path.exists(sdpath):
               sddict2[ sdone ] = sdpath
            else:
               sdpathlist = (re.sub("^\s+", "", ( os.popen('whereis %s'%sdone).read().strip()).split(":")[1])).split(" ")
               print( sdpathlist )
               if sdpathlist != [] and sdpathlist != ['']:
                  sddict2[ sdone ] = ffff[0]
               else: 
                  print( "在环境里找不到 %s，所以手动输入"%sdone )
                  sddict2[ sdone ] = os.path.abspath(input("在这里输入:") )
      configthis.close()
   with open(sd, "w+") as sdthis: #software and database table
      for key, value in sddict2.items():
         sdthis.write("%s: %s # %s \n"%(key, sddict2[key], sddinfo[key] ))
      sdthis.close()
   print("  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ")
   return sddict2


def main():
   print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>...........")
   print("   @@@@@@@@   Aotuconfig.py runnning...")
   parser = argparse.ArgumentParser(description='Demo of argparse')
   parser.add_argument('-w', "--workdir-path", dest = "workp", type=str, default="", help = "工作目录，运行时质控结果、比对结构和变异检测结果等中间和结果文件，您希望保存在服务器的哪个位置")
   parser.add_argument('-d', '--rawdata-path', dest = "rawdatap", type = str, default = "", help ="测序的原始数据，需要告诉程序它，您的数据放在哪里，然后它才摘到去哪找数据开始分析" )
   parser.add_argument('-s', '--sample-information', dest = "sample_info", type = str, default = "", help ="样本信息的表格(\\tab制表符分开)，就是告诉程序它，哪些样本要分析（列名\"sample\"），如果是体细胞变异检测的话，还需要说明癌症样本（列名\"tumor\"）和对应的癌旁(列名\"normal\")" )
   parser.add_argument('-m', '--manager', dest = "mana" , default = "slurm", help = "集群上的任务调度系统" )
   # parser.add_argument('-skp', '--snakemake-path', dest = "snake", default = "snakemake"，help = "如果在您的环境里没有检测到snakemake工具的话，用它添加·x" )
   parser.add_argument('-r', '--run', dest = "run", action='store_true', help="配置好运行的条件之后，跳过用户确认的环节直接开始运行")
   args = parser.parse_args()
   
   
   currentp = os.getcwd()
   pipeplinep = os.path.dirname(os.path.abspath( sys.argv[0] ))
   print(sys.argv[0], pipeplinep)
   checkpipemodules(pipeplinep, ["modules", "scripts", "selfpython", "rules", "work.WGS.yaml", "supplementary.txt", "work.smk" ])
   args.mana = (args.mana).lower()
   # if os.path.exists( args.snake ): snakemakepath = args.snake
  
   printHelp=False
   if args.workp == "" or args.rawdatap == "" or args.sample_info == "": printHelp = True
   if args.workp == "":
      print("\n### 分析的工作目录没有输入呢，这样的话，就不知道在哪个目录下进行分析")
      if pipeplinep != currentp:
         args.workp = os.path.abspath( input("  在这里输入 [直接回车的话会默认当前目录]: ") )
         args.workp = os.getcwd()
      else:
         args.workp = os.path.abspath( input("  在这里输入: ") )
      print(",", args.workp, ",")
   else:
      args.workp = os.path.abspath( args.workp )
   if args.rawdatap == "":
      print("### 分析的原始数据放在没有输入呢")
      args.rawdatap = os.path.abspath( input("  在这里输入：" ) ) 
      print(",", args.rawdatap, ",")
   else:
      args.rawdatap = os.path.abspath( args.rawdatap )
   if args.sample_info == "":
      print("### 样本信息的表还没输入呢")
      args.sample_info = os.path.abspath( input("  在这里输入：") )
      print(args.sample_info)
   else:
      args.sample_info = os.path.abspath(args.sample_info)
   check_sample_information(args.sample_info)
   if printHelp :
      print("\n\n  接下来，执行以下命令来创建工作脚本和配置文件:")
      print(f"Command line:  {pipeplinep}/Autoconfig.py -w {args.workp} -d {args.rawdatap} -s {args.sample_info} -r")
      for i in range(1, 101):
         print("\r", end="")
         print("进度: {}%: ".format(i), "▓" * (i // 2), end="")
         sys.stdout.flush()
         time.sleep(0.05)
      print("\n\n\033[31;40mWRAN: !!!!!上面这条命令行下次可以直接把必须的参数传递进去", end ="，")
      print("它就不会缺什么参数跟你要什么参数了!!!!!\033[0m\n")
      time.sleep(2)
   else:
      print("正在用您提供的参数生成工作脚本和配置文件....")
   config_file_name = "" 
   config_file_name = "work.WGS.yaml"
   work_script = "%s/work.smk"%pipeplinep
   work_config = "%s/%s"%(pipeplinep, config_file_name )
   sddict={ "Pipeline": pipeplinep, }
   sddict2 = checksoftwareanddata(work_config, pipeplinep+"/softwareanddatabase" )
   sddict = dict(**sddict, **sddict2)
   os.system('sed "s/^workp.*/workp = \'%s\'/" %s >%s/work.smk'%(args.workp.replace("/","\/"), 
      work_script.replace("/", "\/"), args.workp))
   os.system('sed "s/^workp: .*/workp: %s/" %s >%s/%s'%(args.workp.replace("/","\/"), 
      work_config.replace("/","\/"), args.workp, config_file_name))
   os.system('sed -i "s/^rawdata: .*/rawdata: %s/" %s/%s'%( args.rawdatap.replace("/","\/"), 
      args.workp, config_file_name))
   os.system('sed -i "s/^sampleinformation:.*$/sampleinformation: %s/"  %s/%s'%(
      (args.sample_info).replace("/","\/"), args.workp, config_file_name))
   
   generatesbatchscript( args.workp, sddict["snakemake"] )
 
   for key, value in sddict.items():
      os.system('sed -i -E "s/%s: .*$/%s: %s/" %s/%s'%(key, key, value.replace("/", "\/"), args.workp, config_file_name ))
   print("分析数据需要的脚本和配置文件生成好了，进入目录%s进行查看" %(args.workp) )
   print("可能,具体的参数设置，需要在 work.WGS.yaml 进行修改，")
   print("确认没有问题就可以开始分析，")
   print("``````````````````````")
   if args.run:
      print("\n\n由于你设置了’-r‘参数，流程自动投递~~")
      print( "        投递：sbatch %s/run.sh"%args.workp )
      os.system( "sbatch %s/run.sh"%args.workp )
      print( "        投递完成！！！！       " )
   else:
      print("由于你没有设置-r参数需要你手动投递：")
      print( "     Command: sbatch %s/run.sh"%args.workp )  
   print("~~~~~~~~~~~~~~~\nTHE END ")
   print("<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<.............")
   
   
   

if __name__ == "__main__":
   main()


#refer 
#1 https://zhuanlan.zhihu.com/p/388930050
#2 https://www.cnblogs.com/shaozelong/p/15911675.html
