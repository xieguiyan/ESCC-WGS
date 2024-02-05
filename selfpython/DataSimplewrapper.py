#!/usr/bin/python3

import os
import re
import random
import string


def updatestatus(status):
   with open(status, "a+") as FFFF:
      FFFF.write("###------------- update   --------- \n")
      FFFF.close()


def mkdir_tmpdirectory(wildcards):
  letters = string.ascii_letters + string.digits
  random_string = 'tmp/'
  while(True):
     for _ in range(20):
        random_string += random.choice(letters)
     if not os.path.exists(random_string):
        os.system(f"mkdir -p {random_string}")
        break
  return random_string


def ReadSampleinformation(sample_info):
   samplelist = []
   tumor_vs_normal = {}
   with open( sample_info, "r+") as FFFF: # 读取表格，虽然pandas好用，但是需要先安装库，
      headers = FFFF.readline().strip().split("\t")
      colnamedict = {"tumor":-1, "normal":-1, "sample":-1}
      for keys in colnamedict.keys():
         if keys in headers : colnamedict[ keys ] = headers.index( keys )
      for oneline in FFFF.readlines():
         allcells = oneline.strip().split("\t")
         if allcells == []: continue
         if colnamedict["normal"] >=0 and colnamedict["tumor"] >=0:
            tumor_vs_normal[ allcells[colnamedict["tumor"]] ] = allcells[colnamedict["normal"]]
         if colnamedict["sample"] >=0: samplelist.append( allcells[colnamedict["sample"]] )
      FFFF.close()
   return samplelist, tumor_vs_normal


def ObtainNodeandSamplelist(config):
   NodeandSample={
      "samplelist": [],
      "tumorlist": [],
      "normallist": [],
      "resultnode2": "NULL",
      "resultnode3": "NULL",  
   }
   tumor_vs_normal = {}
   Sampleinformation = config["Sampleinformation"]
   if isinstance(Sampleinformation, list):
      NodeandSample["samplelist"] = Sampleinformation
   elif isinstance(Sampleinformation, dict):
      tumor_vs_normal = Sampleinformation    
   elif isinstance(Sampleinformation, str) and \
      os.path.exists( Sampleinformation ):
      NodeandSample["samplelist"], tumor_vs_normal = ReadSampleinformation( Sampleinformation )
   else:
      print("ERROR: ")
      exit(-1)
   config["tumor_vs_normal"] = tumor_vs_normal
   NodeandSample["tumorlist"] = list(config["tumor_vs_normal"].keys())
   NodeandSample["normallist"] = list(config["tumor_vs_normal"].values())

    
   if ( "HaplotypeCaller" not in config["MainModule"] ):
      print("tumor vs normal model...")
      if config["tumor_vs_normal"] is None:
         print("Misssing matched condition(Tumor_vs_normal):{NULL}")
         exit(-1)
      NodeandSample["resultnode2"] = "11_mutect2annotation"
      if "ECGEA" in config["MainModule"]:
         resultnode2 = "09_ECannotation"
      NodeandSample["resultnode3"] = NodeandSample["resultnode2"]
      if tumor_vs_normal == {}:
         print("ERROR: ")
         exit(-1)
      NodeandSample["samplelist"] = NodeandSample["tumorlist"] + NodeandSample["normallist"]
   else:
      print("germline matution model...")
      if NodeandSample["samplelist"] == []:
         print("ERROR: ")
         exit(-1) 
      NodeandSample["resultnode2"] = "11_HaAnnotation" # HaplotypeCaller
     
   for key,value in NodeandSample.items(): config[key] = value
   return NodeandSample


def TrueorFalse(para):
  if isinstance(para, bool) and para: return "true"
  elif isinstance(para,str) and para[0].lower()== "t": return "true"
  else: return "false"


def obtainbwapara(config):
   marksecondary = TrueorFalse( config["Parameter"]["bwaparams"]["marksecondary"] )
   KKvalue = config["Parameter"]["bwaparams"]["Kvalue"]
   Kvalue = 0
   if isinstance(KKvalue, str):
      Knums = re.findall("[0-9]+", KKvalue)
      if Knums != [] and int(Knums[0]) >0: Kvalue = int(Knums[0])
   elif isinstance(KKvalue, int) and KKvalue >0: Kvalue =  KKvalue
   config["bwaKvalue"] = Kvalue
   print("bawpara:", marksecondary, Kvalue)
   return marksecondary
    
    

