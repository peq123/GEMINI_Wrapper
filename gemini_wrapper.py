"""

Author           : Paul Quinn
Email            : paul.quinn-2@postgrad.manchester.ac.uk
Github           : https://github.com/peq123/GEMINI_Wrapper
Supervisor       : Dr John Bowes  
Supervisor Email : j.bowes@manchester.ac.uk

Description:

Usage:



"""



### __future__ Imports are always first
from __future__ import absolute_import
from __future__ import print_function

import os, sys,subprocess,  argparse, glob, time, string #,resource
import logging
import logging.handlers
import gzip

from threading  import Thread

try:
    from queue import Queue, Empty
except ImportError:
    from Queue import Queue, Empty  # python 2.x







try:
  import psutil
  psutilavail = True
except:
  psutilavail = False


import shlex
if sys.version_info[0] < 3:
  import ConfigParser
else:
  import configparser
  


runtime = time.time()




"""
------------------------------------------------
            Classes Used by Script
------------------------------------------------
"""
# define Python user-defined exceptions
class ScriptError(Exception):
   """Base class for other exceptions"""
   pass

class VCFFile:
  """ stores information about the VCF file """
  def __init__(self,filepath,raw=False,ped={}):
    self.filepath = filepath
    self.filename = os.path.split(filepath)[1]
    self.fileext = ""
    self.ped = {}
    self.processed = False
    self.compressed = False
    self.raw = raw
    self.indexed = False
    self.indexproc = None
    self.varianttype = "" # for filename use only
    if not filepath == "":
      vtype = []
      for varianttype in ["snp","indel","pindel"]:
        if varianttype.lower() in filepath.lower():
          vtype.append(varianttype)
      self.varianttype = "_".join(vtype)
      self.retrieveinfo() 
      pass 
    if len(ped) > 0:
      self.ped = { k : ped.get(k,pedvalue) for k,pedvalue in self.ped.items() if k in self.ped.keys()}
      

    
  def get_filehandle(self):
    if self.filepath.endswith("vcf.gz"):
      self.fileext = ".vcf.gz"
      self.compressed = True
      return gzip.open(self.filepath,"rb")
    else:
      self.fileext = ".vcf"
      return open(self.filepath,"r")

  def retrieveinfo(self):
    with self.get_filehandle() as fh:
      # opened the file now read all comments until the column headers are found      
      for line in fh:
        if isinstance(line,str):
          line = line.strip()
        else:
          line = line.decode("utf-8").strip()
        if not line.startswith("#") or line.startswith("##"):
          # skip meta info
          if "snpeff" in line.lower():
            self.processed = True
          continue
        self.ped = { id:[id,id,"-1","-1","-1","-1"] for id in line.split() if id.lower() not in ["#chrom","pos","id","ref","alt","qual","filter","info","format"]}
        break # no more lines need to be read
  def index(self):

    file_start_time = time.time()
    tabix = ChildProcess("tabix -p vcf {0}".format(self.filepath))    
    self.indexed = True
    ChildProcess.poll_processes([tabix]) 
    tabix.reportusage(16)
    scriptlogger.log(19,"Time taken:{0}".format(time.time() - file_start_time))
    return tabix
  def remove(self):
    if not self.raw:
      try:
        os.remove(self.filepath)
      except:
        pass
      for p in glob.glob(self.filepath + "*i"):
        try:
          os.remove(p)
        except:
          pass
      
  def idcount(self):
    return len(self.ped.keys())
  
  def families(self):
    return set([v[0] for k,v in self.ped.items()] )

  def move(self,destination,overwrite=False,raw=False):
    destination = check_directory(destination)
    try: 
      relatedfilelist = glob.glob(self.filepath + "*i")
      if not self.filepath in relatedfilelist:
        relatedfilelist += [self.filepath]
        
      for f in relatedfilelist:
        fname = os.path.split(f)[1]
        if not overwrite:
          newpath = os.path.join(destination,"{0}_{1}".format(executiontime,fname ))
        else:
          newpath = os.path.join(destination,fname)
          if os.path.exists(newpath):
            os.remove(newpath)    
        os.rename(f,newpath)
        if f == self.filepath:
          self.filepath = newpath
    except:
      pass
    self.raw = self.raw or raw

  @classmethod
  def movelist(cls,action,destination,samplelist):
    destination  = check_directory(destination)
    if destination != "" :
      scriptlogger.info("Moving {0} file(s) to {1}".format(action,destination))
      for id in samplelist:
        for vcf in samplelist[id]:
          vcf.move(destination,overwrite=subargs.overwrite,raw=True)
  

class ChildProcess():
  currentopenfile = []
  cycle = 0 
  
  def __init__(self,cmd,communicate=False,stdin=None,stdout=subprocess.PIPE,stderr=subprocess.PIPE,pollout=False,pollerr=True):
    """ Opens subprocess. Depending on passed communicate, the process open object or final outputs will be returned """
    # set up attributes
    pollthread = None
    self.outq = None
    self.errq = None
    self.stdout = ""
    self.stderr = ""
    self.rss = 0 
    self.vms = 0
    self.cpu = 0 
    self.pollout = pollout
    self.pollerr = pollerr
    self.outcallback = None
    self.errcallback = None
    
    
    if type(cmd) is list:
      cmd = " ".join(cmd)
    cmd = " ".join(shlex.split(cmd,posix=False))
    self.cmd = cmd
    self.process = subprocess.Popen(cmd,stdin=stdin,stdout=stdout,stderr=stderr,shell=True)  
    if communicate:      
      (self.stdout,self.stderr) = self.process.communicate()
      if not isinstance(self.stdout,str):
        self.stdout = self.stdout.decode("utf-8")
      if not isinstance(self.stderr,str):
        self.stderr = self.stderr.decode("utf-8")
    
  def get_child_use(self,py):
    crss = 0 
    cvms = 0
    try:
      pychildren = py.children(recursive=True)
      if len(pychildren) > 0:
        for child in pychildren:
          crss += child.memory_info()[0]
          cvms += child.memory_info()[1]
    except:
      pass
    return crss,cvms

  def capture_stream(self,out, queue):
    for line in iter(out.readline, b''):
        queue.put(line)
    out.close()

  #def capture_memory(pid):


  ## threading used to non-block 
  def start_read_thread(self,p):
    q = Queue()
    try:      
      t = Thread(target=self.capture_stream, args=(p, q))
      t.daemon = True # thread dies with the program
      t.start()  
    except:
      pass
    return q
  
  def startpolling(self):
    self.cpu = time.time()
    if self.pollout :
      self.outq = self.start_read_thread(self.process.stdout)
    if self.pollerr:
      self.errq = self.start_read_thread(self.process.stderr)
    try:      
      t = Thread(target=self.poll)
      t.daemon = True # thread dies with the program
      t.start()  
    except:
      pass

  
  def poll(self):
    """ this function setups the polling for this process, used if need to capture memory """
    global psutilavail
    self.cpu = time.time()
    while self.process.returncode is None:
      
      self.process.poll()
      try:
        if psutilavail:          
          py = psutil.Process(self.process.pid)
          crss = py.memory_info()[0]
          cvms = py.memory_info()[1]
          # check if the the process spawns and children and calculate them in the mumuse
          childrss,childvms = self.get_child_use(py)
          crss += childrss
          cvms += childvms 
          if self.rss <= crss :            
            self.rss = crss
          if self.vms <= cvms :            
            self.vms = cvms

      except:
        pass
      if self.outq is not None:
         
        try:        
          while 1:
            line = self.outq.get_nowait()
            if not isinstance(line,str):
              line = line.decode("utf-8")
            if self.outcallback is None:
              self.stdout += line
            else:
              self.outcallback(line)
        except Empty:
          pass
    
      if self.errq is not None:
        try:        
          while 1:
            line = self.errq.get_nowait()
            if not isinstance(line,str):
              line = line.decode("utf-8")
            self.stderr += line
        except Empty:
            pass
    self.cpu = time.time() - self.cpu

  def reportusage(self,channel):
    scriptlogger.log(channel,"cmd:{0} cpu:{1} rss:{2} vms:{3}".format(self.cmd,self.cpu,self.rss,self.vms) )



  @classmethod    
  def openfile(cls,path,flags):
     if not isinstance(cls.currentopenfile,list):
       cls.currentopenfile = []
     cls.currentopenfile.append(path)
     return open(path,flags)

  @staticmethod
  def poll_processes(procs):
    """
    Poll passed processes
    """

    active =  len(procs)

  # start procs 
    
    for proc in procs:
      
      proc.startpolling()

    maxmsg = 0

    while active > 0 :  
      active = 0 
      message = []
      for i,proc in  enumerate(procs):      
        if proc.process.returncode == None:
          active += 1
      
      Spinner.spin()

      time.sleep(.25)
 


class Spinner():
  
  cycle = 0 
  @classmethod
  def spin(cls):    
    cls.cycle = ((cls.cycle + 1) % 4) if cls.cycle is not None else 0
    print("{0}\033[F".format(["-","\\","|","/"][cls.cycle]))
    
class ReportWorkSheet():
  def __init__(self,ws):
    self.worksheet = ws
    self.rowoffset = 0 
  
  def adddata(self,data):
    if isinstance(data,list):
      for i,d in enumerate(data):
        d = d.replace(";",";\n").replace(",",",\n").replace("\n\n","\n").strip()
        try:
          testfloat = float(d)
          self.worksheet.write(self.rowoffset,i,"={0}".format(d))
        except:
          self.worksheet.write(self.rowoffset,i,d)
      self.rowoffset += 1
"""
------------------------------------------------
              Globals
------------------------------------------------              
"""
config = {} # for script config
globalargs = None # for args parse
subargs = None # for args parse
wrappedcmds = []# ["","load","query","comp_hets","mendelian_error","de_novo","autosomal_recessive","autosomal_dominant"] 
uniquecmds = ["id","consolidate","analysis"]
scriptpath = os.path.dirname(os.path.realpath(__file__))
cfgfilepath = os.path.abspath(scriptpath + "/gemini_wrapper.cfg")
workingtmpfolder = "" # set after config parse
logpath = ""
geminiapi = False
bcftools = False
builtincmds = [] # defined later
notwrapped = [] # defined later
outpipe = None
builtinargdict ={}
## for clean up
samplelist = {}
currentopenfile = "" 
##
executiontime = "{0:04d}{1:02d}{2:02d}{3:02d}{4:02d}{5:02d}".format(time.localtime().tm_year,
                                                                          time.localtime().tm_mon,
                                                                          time.localtime().tm_mday,
                                                                          time.localtime().tm_hour,
                                                                          time.localtime().tm_min,
                                                                          time.localtime().tm_sec)

## initailise log format
formatter = logging.Formatter('%(levelname)s : %(message)s')

# initailise logger and clear any handles already set
scriptlogger = logging.getLogger()
if len(scriptlogger.handlers)> 0 :
  for handler in scriptlogger.handlers[:]:
      scriptlogger.removeHandler(handler)  
scriptlogger.setLevel(logging.DEBUG)
#scriptlogger.setFormatter(formatter)

# initialise logger

  
# reinitalise the terminal logs
termlog = logging.StreamHandler(sys.stderr)
termlog.setLevel(logging.DEBUG)
termlog.setFormatter(formatter)
scriptlogger.addHandler(termlog)
logging.addLevelName(25, "INFO")

# create custom logs for the programs that provide stats in output
logging.addLevelName(11, "BCFTOOLS")
logging.addLevelName(12, "SED")
logging.addLevelName(13, "VT")
logging.addLevelName(14, "SNPEFForVEP")
logging.addLevelName(15, "BGZIP")
logging.addLevelName(16, "TABIX")
logging.addLevelName(18, "GEMINI")
logging.addLevelName(19, "TIME")







"""
------------------------------------------------
            Config/Template Functions
------------------------------------------------
"""



def parse_config(thepath,isstring=False):
  """ Parses an ini style file using configparser """

  def removeemptyvalues(thedict):
    out = {}
    for k,v in thedict.items():
      if type(v) is dict:
        tmp = removeemptyvalues(v)
        if len(tmp) > 0:
          out[k.lower()] = tmp
      elif type(v) is not str or (type(v) is str and len(v) > 0):
        out[k.lower()] = v
    return out

  if sys.version_info[0] < 3:
    configobj = ConfigParser.RawConfigParser({},dict,allow_no_value=True)
  else:
    configobj = configparser.RawConfigParser({},dict,allow_no_value=True)
  if isstring:
    configobj.read_string(thepath)
  else:
    configobj.read(thepath)



  return  removeemptyvalues(configobj._sections)

def parse_envconfig():

  """ Parses environment configuration file """
  global config
  def mergeDict(dict1, dict2):
    ''' Merge dictionaries and keep values of common keys in list'''
    for key, value in dict1.items():
      if isinstance(value,dict):
        dict1[key] = mergeDict(value,dict2.get(key,value))
      else:        
        dict1[key] = dict2.get(key,value)
        # dont allow blank if default is not 
        if dict1[key] == "" and value != "" :
          dict1[key] = value
    return dict1
  defaults = {
              "env":{
                  "anaconda": "",
                  "gemini":"",
                  "bcftools":"",
                  "vt" : "",
                  "snpeff":"",
                  "vep":"",
                  "genome_reference":"",
                  "switch_to_anaconda":"True",
                  "permissions":664
                  },
              "load":{
                "consolidation":"family",
                "cores":"3"               
              },
              "paths": {
                "tmp_dir":os.path.join(scriptpath,"tmp"),
                "log_dir":scriptpath,
                "database_dir":os.path.join(scriptpath,"databases"),                
                "consolidated_dir":"",
                "prepped_dir":"",                
                "template_dir":"",
                "report_dir":""
              }
              
            }
  
  # ensure all values are present
  config = parse_config(cfgfilepath)
  config = mergeDict(defaults,config)

  # expand user paths 
  for id,p in config["paths"].items():    
    config["paths"][id] = check_directory(p)

  config["env"]["snpeff"] = parse_paths(config["env"]["snpeff"])
  config["env"]["genome_reference"] = parse_paths(config["env"]["genome_reference"])


  # check permissions 
  config["env"]["permissions"] = int(''.join( [ c if int(c) <=7 else "7" for c in config["env"].get("permissions","666")[-3::] ] ))



  # check consolidation 
  if not config["load"]["consolidation"] in ["none","sample","family","all"]:
    config["load"]["consolidation"] = "none"

def generate_envconfig(thepath):
  """ Generates configuration file without overriding existing parameters """
  
  
  # write the config file
  with open(thepath,"w") as configfile:
    for key, value in config.items():
      configfile.write("[{0}]\n".format(key))
      for key2, value2 in value.items():
        if key2.startswith("#"):
          configfile.write("{0} {1}\n".format(key2,value2))
        else:
          configfile.write("{0} = {1}\n".format(key2,value2))


def parse_args():
  """ Parses incoming command line interface """

  global globalargs, subargs, remainingargs,builtincmds,builtinargdict, termlog
  global scriptlogger

  def addsubcommand(indict):         
    kwargs = indict["kwargs"]
    currentparser = subparsers.add_parser(indict["cmd"], **kwargs)
    
    for i in indict["arguments"]:
      currentparser.add_argument(*i["args"],**i["kwargs"])

    return currentparser
  


  # reinitalise the parser
  parentparser = argparse.ArgumentParser(description="GEMINI Wrapper\nProvides a complete workflow solution with GEMINI at its heart", 
                                
                                add_help=False)

  

  #parser.add_argument("--manpage",                dest="manpage",    action="store_true",                             help='Show more detailed \"manpage\" ')


  ## first define global arguments
  parentparser.add_argument("-h", "--help" ,  dest="help",       action="store_true", help="Show this help message and exit.")
  parentparser.add_argument("--saveconfig",   dest="saveconfig", action="store_true", help="Save configuration. Useful if you have previously removed unnecessary parameters")
  parentparser.add_argument("--logfile",      dest="logfile",    action="store_true", help="Produce a logfile named time-of-execution.log in script directory")
  parentparser.add_argument("--usage",         dest="usage",     action="store_true", help="Produce a memory and cpu useage of child precesses")
  parentparser.add_argument("--overwrite",action="store_true",dest="overwrite",help="If provided any files with existing names will be overwritten. WARNING: Does not apply to file specified by --output.")
  parentparser.add_argument("--permissions",dest="permissions",default=config["env"]["permissions"],type=int,help="The permisions for newly created files" )

  ## build the mutually exclusive group for output
  


  group = parentparser.add_mutually_exclusive_group()
  group.add_argument(  "--debug",         dest="debug",     action="store_true",                            help="Displays debug messages to standard output.")
  group.add_argument(  "--verbose",         dest="verbose",     action="store_true",                            help="Displays info messages to standard output.")
  group.add_argument(  "--silent",          dest="silent",      action="store_true",                            help="Removes all messages from standard output")

  globalparser = argparse.ArgumentParser(description="GEMINI Wrapper\nProvides a complete workflow solution with GEMINI at its heart", 
                                
                                add_help=False,parents = [parentparser])

  
  globalparser.add_argument("cmd"  ,  nargs="?",default=""    ,          help="Command")
  globalparser.add_argument("--envreloaded" ,  dest="envreloaded",   help="Environment Reload Pass")
  globalparser.add_argument("--child" ,action="store_true",    dest="childprocess",   help="Marks script run as a child for parellelisation")
  
  globalargs,remainingargs = globalparser.parse_known_args()

  subcmdparser = argparse.ArgumentParser(description="GEMINI Wrapper\nProvides a complete workflow solution with GEMINI at its heart", 
                                
                                add_help=False,parents=[parentparser]) 
 
  oparser = argparse.ArgumentParser(add_help=False)
  oparser.add_argument("--out", dest="output"  , default="" ,        help="Directory path for the output (either database, xlsx report depending on wrap)")
  subparsers = subcmdparser.add_subparsers(dest="cmd")
  
  parentparser.add_argument("source" , action="store"  ,  nargs="+",     help="Shell pattern to load one or more files (VCF or database depending on command")
  
  
    #unknowncmd.add_argument("--stdout",dest="stdout",action="store_true",     help="Will write results to stdout rather than excel file ")
    



  # prepare the wrapper subcommands


  # "ID"                            
  parser_id  = subparsers.add_parser('id', help='Provide a file pattern and the script will extract the sample IDs from the vcf files and display them', 
                                description="Use this command to extact sample IDs from vcf files",
                                add_help=False,parents = [parentparser,oparser]) 
  
    #parser_id.add_argument("--rename",dest="rename",  nargs="?", help="Supply a CSV file containing current sample ids and the script will update the vcfs with the new ids. If the filename contains the ID, it will also be updated.")
  

  # TODO
  # CONSOLIDATE
  # PREP  
  parser_consolidate = subparsers.add_parser('consolidate', help='Provide a file pattern and the script will consolidate as requested', 
                                description="Use this command to extact sample IDs from vcf files, or rename samples by supplying a csv file.",
                                add_help=False,parents = [parser_id]) 
  parser_consolidate.add_argument("--consolidate",dest="consolidate", default=config["load"]["consolidation"].lower() ,choices=["none","sample","family","all"],help="Use BCF Tools to merge vcfs before processing. Merging is either done by sample, family or all.")

  #parser_prep  = subparsers.add_parser('prep', help='Provide a file pattern and the script will prep them to be loaded into GEMINI', 
  #                              
  #                              add_help=False,parents = [dbparser,parser_consolidate])#

  # "Load"
  parser_load  = subparsers.add_parser('load', help='Instructs GEMINI to load vcf data, first the wrapper will prep each the VCF(s) and if required generate PED files', 
                                
                                add_help=False,parents = [parser_consolidate])

  parser_load.add_argument("--cores",dest="cores", default=int(config.get("load",{}).get("cores",3)),type=int, help="Define how many cores GEMINI should ultilise in the load procedure")
  
  
  parser_load.add_argument("-t",dest="t", help="This should not be passed, it is handled by wrapper")
  parser_load.add_argument("-v",dest="v", help="This should not be passed, it is handled by wrapper")
  #parser_load.add_argument("--skipprep",dest="existingvalid",action="store_true", help="Use this flag to use existing processed files, useful if previous run halted without loading into GEMINI")

  #parser_load.add_argument("source" , action="store"  ,  nargs="+",                  help="Shell pattern to load one or more files")
  parser_load.add_argument("-p","--ped",dest="ped" , action="store" , help="PED file (or CSV file) containing sample inforamtion. Described on http://zzz.bwh.harvard.edu/plink/data.shtml")
  
  
  
  






  # gemini built-in common arguments 
  analysis = argparse.ArgumentParser(add_help=False,parents = [parentparser,oparser])
  analysis.add_argument("--format",dest="format", default="XLSX", choices=["XLSX","TSV","CSV"] , help="File format of results report")
  analysis.add_argument("--template",dest="queryfile",default="", help="Pass a text file containing analysis information",required=globalargs.cmd=="analysis")
  analysis.add_argument("--listheading",action="store_true", help="Will display all tables and columns in database")
  
  helpfunctions = {"id":parser_id,
                  "consolidate":parser_consolidate,
                  #"prep":parser_prep,
                  "load":parser_load       
                    }
  builtincmds = [
        {
      "cmd":"analysis",
      "kwargs": {"help":'Run an analysis using a file. Can run multiple commands from one file', "add_help":False,"parents":[analysis]},
      "arguments":[
         
      ]    
      },
      {
      "cmd":"query",
      "kwargs": {"help":'GEMINI analysis. See below ', "add_help":False,"parents":[analysis]},
      "arguments":[
         {
          "args":["-q"],
           "kwargs": {"help":"See below"}
        },
        {
          "args":["--gt-filter"],
           "kwargs": {"help":"See below",}
        },
        {
          "args":["--show-samples"],
           "kwargs": {"help":"See below","action":"store_true"}
        },
        {
          "args":["--show-families"],
           "kwargs": {"help":"See below","action":"store_true"}
        },
        {
          "args":["--family-wise"],
           "kwargs": {"help":"See below","action":"store_true"}
        },
        {
          "args":["--min-kindreds"],
           "kwargs": {"help":"See below"}
        },
        {
          "args":["--sample-delim"],
           "kwargs": {"help":"See below"}
        },
        {
          "args":["--header"],
           "kwargs": {"help":"Will always be passed","action":"store_true"}
        },
        {
          "args":["--sample-filter"],
           "kwargs": {"help":"See below"}
        },
        {
          "args":["--in"],
           "kwargs": {"help":"See below","nargs":"+"}
        },
        {
          "args":["--region"],
           "kwargs": {"help":"See below"}
        },
        {
          "args":["--carrier-summary-by-phenotype"],
           "kwargs": {"help":"See below"}
        },
        {
          "args":["--dgidb"],
           "kwargs": {"help":"See below","action":"store_true"}
        },
        {
          "args":["--use-bcolz"],
           "kwargs": {"help":"See below","action":"store_true"}
        }
      ]    
      },
      {   "cmd":"amend",
      "kwargs": {"help":'GEMINI Builtin analysis. See below ', "add_help":False,"parents":[parentparser]},
      "arguments":[
        {
          "args":["--sample"],
           "kwargs": {"help":"See below"}
        },
        {
          "args":["--clear"],
           "kwargs": {"help":"See below"}
        }
        ]
    },
    {
       "cmd":"actionable_mutations",
      "ignore":True,
      "kwargs": {"help":'GEMINI Builtin command. See below ', "add_help":False,"parents":[parentparser]},
      "arguments":[
        
      ]
    },
    
    {   "cmd":"annotate",
      "kwargs": {"help":'GEMINI Builtin analysis. See below ', "add_help":False,"parents":[parentparser]},
      "arguments":[
        {
          "args":["-f"],
           "kwargs": {"help":"See below"}
        },
        {
          "args":["-c"],
           "kwargs": {"help":"See below"}
        },
        {
          "args":["-a"],
           "kwargs": {"help":"See below","choices":["boolean","count", "extract"]}
        },
        {
          "args":["-e"],
           "kwargs": {"help":"See below"}
        },
       {
          "args":["-t"],
           "kwargs": {"help":"See below"}
        },
       {
          "args":["-o"],
           "kwargs": {"help":"See below"}
        },
       {
          "args":["--region-only"],
           "kwargs": {"help":"See below","action":"store_true"}
        }
      ]
    },{"cmd":    "autosomal_recessive",
      "kwargs": {"help":'GEMINI Builtin analysis. See below ', "add_help":False,"parents":[analysis]},
      "arguments":[{
          "args":["--columns"],
           "kwargs": {"help":"See below"}
        },
        {
          "args":["--filter"],
           "kwargs": {"help":"See below","dest":"gtfilterrequired","action":"append"}
        },
        {
          "args":["--min-kindreds"],
           "kwargs": {"help":"See below"}
        },
        {
          "args":["--familes"],
           "kwargs": {"help":"See below"}
        },
        {
          "args":["--allow-unaffected"],
           "kwargs": {"help":"See below","action":"store_true"}
        },
        {
          "args":["-d"],
           "kwargs": {"help":"See below","metavar":"MIN_SAMPLE_DEPTH"}
        },
        {
          "args":["--min-gq"],
           "kwargs": {"help":"See below"}
        },
        {
          "args":["--gt-pl-max"],
           "kwargs": {"help":"See below","metavar":"GT_PHRED_LL"}
        }]
    },
    {"cmd": "autosomal_dominant",
      "kwargs": {"help":'GEMINI Builtin analysis. See below ', "add_help":False,"parents":[analysis]},
      "arguments":[{
          "args":["--columns"],
           "kwargs": {"help":"See below"}
        },
        {
          "args":["--filter"],
           "kwargs": {"help":"See below","dest":"gtfilterrequired","action":"append"}
        },
        {
          "args":["--min-kindreds"],
           "kwargs": {"help":"See below"}
        },
        {
          "args":["--familes"],
           "kwargs": {"help":"See below"}
        },
        {
          "args":["--allow-unaffected"],
           "kwargs": {"help":"See below","action":"store_true"}
        },
        {
          "args":["-d"],
           "kwargs": {"help":"See below","metavar":"MIN_SAMPLE_DEPTH"}
        },
        {
          "args":["--min-gq"],
           "kwargs": {"help":"See below"}
        },
        {
          "args":["--gt-pl-max"],
           "kwargs": {"help":"See below","metavar":"GT_PHRED_LL"}
        }]
    },
    {"cmd":"comp_hets",
      "kwargs": {"help":'GEMINI Builtin analysis. See below ', "add_help":False,"parents":[analysis]},
      "arguments":[
         {
          "args":["--columns"],
           "kwargs": {"help":"See below"}
        },
        {
          "args":["--filter"],
           "kwargs": {"help":"See below"}
        },
        {
          "args":["--min-kindreds"],
           "kwargs": {"help":"See below"}
        },
        {
          "args":["--familes"],
           "kwargs": {"help":"See below"}
        },
        {
          "args":["--allow-unaffected"],
           "kwargs": {"help":"See below","action":"store_true"}
        },
        {
          "args":["-d"],
           "kwargs": {"help":"See below","metavar":"MIN_SAMPLE_DEPTH"}
        },
        {
          "args":["--min-gq"],
           "kwargs": {"help":"See below"}
        },
        {
          "args":["--gt-pl-max"],
           "kwargs": {"help":"See below","metavar":"GT_PHRED_LL"}
        },
        {
          "args":["--gene-where"],
           "kwargs": {"help":"See below","metavar":"WHERE"}
        },
        {
          "args":["--pattern-only"],
           "kwargs": {"help":"See below","action":"store_true"}
        },
        {
          "args":["--max-priority"],
           "kwargs": {"help":"See below"}
        }
      ]    
    },
    {"cmd":"mendel_errors",
      "kwargs": {"help":'GEMINI Builtin analysis. See below ', "add_help":False,"parents":[analysis]},
      "arguments":[{
          "args":["--columns"],
           "kwargs": {"help":"See below"}
        },
        {
          "args":["--filter"],
           "kwargs": {"help":"See below","dest":"gtfilterrequired","action":"append"}
        },
        {
          "args":["--min-kindreds"],
           "kwargs": {"help":"See below"}
        },
        {
          "args":["--familes"],
           "kwargs": {"help":"See below"}
        },
        {
          "args":["--lenient"],
           "kwargs": {"help":"See below","action":"store_true"}
        },
        {
          "args":["-d"],
           "kwargs": {"help":"See below","metavar":"MIN_SAMPLE_DEPTH"}
        },
        {
          "args":["--min-gq"],
           "kwargs": {"help":"See below"}
        },
        {
          "args":["--gt-pl-max"],
           "kwargs": {"help":"See below","metavar":"GT_PHRED_LL"}
        },
        {
          "args":["--only-unaffected"],
           "kwargs": {"help":"See below","action":"store_true"}
        }]
    },
    { "cmd": "de_novo",
      "kwargs": {"help":'GEMINI Builtin analysis. See below ', "add_help":False,"parents":[analysis]},
      "arguments":[{
          "args":["--columns"],
           "kwargs": {"help":"See below"}
        },
        {
          "args":["--filter"],
           "kwargs": {"help":"See below"}
        },
        {
          "args":["--min-kindreds"],
           "kwargs": {"help":"See below"}
        },
        {
          "args":["--familes"],
           "kwargs": {"help":"See below"}
        },
        {
          "args":["--lenient"],
           "kwargs": {"help":"See below","action":"store_true"}
        },
        {
          "args":["--allow-unaffected"],
           "kwargs": {"help":"See below","action":"store_true"}
        },
        {
          "args":["-d"],
           "kwargs": {"help":"See below","metavar":"MIN_SAMPLE_DEPTH"}
        },
        {
          "args":["--min-gq"],
           "kwargs": {"help":"See below"}
        },
        {
          "args":["--gt-pl-max"],
           "kwargs": {"help":"See below","metavar":"GT_PHRED_LL"}
        }]
    },
    
    {"cmd":  "x_linked_recessive",
      "kwargs": {"help":'GEMINI Builtin analysis. See below ', "add_help":False,"parents":[analysis]},
      "arguments":[{
          "args":["--columns"],
           "kwargs": {"help":"See below"}
        },
        {
          "args":["--filter"],
           "kwargs": {"help":"See below","dest":"gtfilterrequired","action":"append"}
        },
        {
          "args":["--min-kindreds"],
           "kwargs": {"help":"See below"}
        },
        {
          "args":["--familes"],
           "kwargs": {"help":"See below"}
        },
        {
          "args":["--allow-unaffected"],
           "kwargs": {"help":"See below","action":"store_true"}
        },
        {
          "args":["-d"],
           "kwargs": {"help":"See below","metavar":"MIN_SAMPLE_DEPTH"}
        },
        {
          "args":["--min-gq"],
           "kwargs": {"help":"See below"}
        },
        {
          "args":["--gt-pl-max"],
           "kwargs": {"help":"See below","metavar":"GT_PHRED_LL"}
        }]
    },
    {"cmd": "x_linked_dominant",
      "kwargs": {"help":'GEMINI Builtin analysis. See below ', "add_help":False,"parents":[analysis]},
      "arguments":[{
          "args":["--columns"],
           "kwargs": {"help":"See below"}
        },
        {
          "args":["--filter"],
           "kwargs": {"help":"See below","dest":"gtfilterrequired","action":"append"}
        },
        {
          "args":["--min-kindreds"],
           "kwargs": {"help":"See below"}
        },
        {
          "args":["--familes"],
           "kwargs": {"help":"See below"}
        },
        {
          "args":["--allow-unaffected"],
           "kwargs": {"help":"See below","action":"store_true"}
        },
        {
          "args":["-d"],
           "kwargs": {"help":"See below","metavar":"MIN_SAMPLE_DEPTH"}
        },
        {
          "args":["--min-gq"],
           "kwargs": {"help":"See below"}
        },
        {
          "args":["--gt-pl-max"],
           "kwargs": {"help":"See below","metavar":"GT_PHRED_LL"}
        }]
    },
    {"cmd":    "x_linked_de_novo",
      "kwargs": {"help":'GEMINI Builtin analysis. See below ', "add_help":False,"parents":[analysis]},
      "arguments":[{
          "args":["--columns"],
           "kwargs": {"help":"See below"}
        },
        {
          "args":["--filter"],
           "kwargs": {"help":"See below","dest":"gtfilterrequired","action":"append"}
        },
        {
          "args":["--min-kindreds"],
           "kwargs": {"help":"See below"}
        },
        {
          "args":["--familes"],
           "kwargs": {"help":"See below"}
        },
        {
          "args":["--allow-unaffected"],
           "kwargs": {"help":"See below","action":"store_true"}
        },
        {
          "args":["-d"],
           "kwargs": {"help":"See below","metavar":"MIN_SAMPLE_DEPTH"}
        },
        {
          "args":["--min-gq"],
           "kwargs": {"help":"See below"}
        },
        {
          "args":["--gt-pl-max"],
           "kwargs": {"help":"See below","metavar":"GT_PHRED_LL"}
        }]
    },
    {   "cmd":   "gene_wise",
      "kwargs": {"help":'GEMINI Builtin analysis. See below ', "add_help":False,"parents":[analysis]},
      "arguments":[
        {
          "args":["--min-filters"],
           "kwargs": {"help":"See below","dest":"minfilters","default":1,"type":int}
        },
        {
          "args":["--gt-filter"],
           "kwargs": {"help":"See below","dest":"gtfilterrequired","action":"append"}
        },
        {
          "args":["--gt-filter-required"],
           "kwargs": {"help":"See below","dest":"gtfilterrequired","action":"append"}
        },
        {
          "args":["--where"],
           "kwargs": {"help":"See below"}
        },
        {
          "args":["--filter"],
           "kwargs": {"help":"See below","dest":"gtfilterrequired","action":"append"}
        },
        {
          "args":["--columns"],
           "kwargs": {"help":"See below"}
        }
      ]
    },
    {   "cmd":     "pathways",
      "kwargs": {"help":'GEMINI Builtin analysis. See below ', "add_help":False,"parents":[analysis]},
      "arguments":[
        {
          "args":["-v"],
           "kwargs": {"help":"See below","dest":"version","metavar":"STRING"}
        },
        {
          "args":["--lof"],
           "kwargs": {"help":"See below"}
        }
      ]
    },    
    {   "cmd": "interactions",
      "kwargs": {"help":'GEMINI Builtin analysis. See below ', "add_help":False,"parents":[analysis]},
      "arguments":[
        {
          "args":["-g"],
           "kwargs": {"help":"See below"}
        },
        {
          "args":["-r"],
           "kwargs": {"help":"See below","default":1,"type":int}
        },        
        {
          "args":["--edges"],
           "kwargs": {"help":"See below"}
        },        
        {
          "args":["--var"],
           "kwargs": {"help":"See below"}
        }
      ]
    },
    {   "cmd":"lof_interactions",
      "kwargs": {"help":'GEMINI Builtin analysis. See below ', "add_help":False,"parents":[analysis]},
      "arguments":[
        {
          "args":["-r"],
           "kwargs": {"help":"See below","default":1,"type":int}
        },
        
        {
          "args":["--edges"],
           "kwargs": {"help":"See below"}
        },
        
        {
          "args":["--var"],
           "kwargs": {"help":"See below"}
        }
      ]
    },

    {   "cmd":"lof_sieve",
      "kwargs": {"help":'GEMINI Builtin analysis. See below ', "add_help":False,"parents":[analysis]},
      "arguments":{}
    },

    
    {   "cmd":"region",
      "kwargs": {"help":'GEMINI Builtin analysis. See below ', "add_help":False,"parents":[analysis]},
      "arguments":[
        {
          "args":["--reg"],
           "kwargs": {"help":"See below"}
        },
        {
          "args":["--gene"],
           "kwargs": {"help":"See below"}
        },
       {
          "args":["--header"],
           "kwargs": {"help":"See below","action":"store_true"}
        },
       {
          "args":["--columns"],
           "kwargs": {"help":"See below"}
        },
        {
          "args":["--filter"],
           "kwargs": {"help":"See below"}
        },
       {
          "args":["--show-samples"],
           "kwargs": {"help":"See below","action":"store_true"}
        }
          ]
    },
    {   "cmd":"windower",
      "kwargs": {"help":'GEMINI Builtin analysis. See below ', "add_help":False,"parents":[analysis]},
      "arguments":[
        {
          "args":["-w"],
           "kwargs": {"help":"See below"}
        },
        {
          "args":["-s"],
           "kwargs": {"help":"See below"}
        },
        {
          "args":["-t"],
           "kwargs": {"help":"See below","choices":["nucl_div","hwe"]}
        },
        {
          "args":["-o"],
           "kwargs": {"help":"See below","choices":["mean","median","min","max","collapse"]}
        }
      ]
    },
    {   "cmd":"stats",
      "kwargs": {"help":'GEMINI Builtin analysis. See below ', "add_help":False,"parents":[analysis]},
      "arguments":[
        {
          "args":["--tstv-coding"],
           "kwargs": {"help":"See below","action":"store_true"}
        },
        {
          "args":["--tstv-noncoding"],
           "kwargs": {"help":"See below","action":"store_true"}
        },
        {
          "args":["--snp-counts"],
           "kwargs": {"help":"See below","action":"store_true"}
        },
        {
          "args":["--sfs"],
           "kwargs": {"help":"See below","action":"store_true"}
        },
        {
          "args":["--mds"],
           "kwargs": {"help":"See below","action":"store_true"}
        },
        {
          "args":["--gts-by-sample"],
           "kwargs": {"help":"See below","action":"store_true"}
        },
        {
          "args":["--vars-by-sample"],
           "kwargs": {"help":"See below","action":"store_true"}
        },
        {
          "args":["--summarize"],
           "kwargs": {"help":"See below"}
        },
        {
          "args":["--gt-filter"],
           "kwargs": {"help":"See below"}
        }
     ]
    },
    {   "cmd":"burden",
      "kwargs": {"help":'GEMINI Builtin analysis. See below ', "add_help":False,"parents":[analysis]},
      "arguments":[
       {
          "args":["--nonsynonymous"],
           "kwargs": {"help":"See below","action":"store_true"}
        },
        {
          "args":["--cases"],
           "kwargs": {"help":"See below","nargs":"+"}
        },
        {
          "args":["--controls"],
           "kwargs": {"help":"See below","nargs":"+"}
        },
         
        {
          "args":["--calpha"],
           "kwargs": {"help":"See below","action":"store_true"}
        },
        {
          "args":["--min-aaf"],
           "kwargs": {"help":"See below"}
        },
        {
          "args":["--max-aaf"],
           "kwargs": {"help":"See below"}
        },
         
        {
          "args":["--save_tscores"],
           "kwargs": {"help":"See below","action":"store_true"}
        }
        
      ]
    },
    {   "cmd":"set_somatic",
      "kwargs": {"help":'GEMINI Builtin analysis. See below ', "add_help":False,"parents":[analysis]},
      "arguments":[
       {
          "args":["--min-depth"],
           "kwargs": {"help":"See below"}
        },
        {
          "args":["--min-qual"],
           "kwargs": {"help":"See below"}
        },
        {
          "args":["--min-somatic-score"],
           "kwargs": {"help":"See below"}
        },
        {
          "args":["--max-norm-alt-freq"],
           "kwargs": {"help":"See below"}
        },
        {
          "args":["--max-norm-alt-count"],
           "kwargs": {"help":"See below"}
        },
        {
          "args":["--min-norm-depth"],
           "kwargs": {"help":"See below"}
        },
        {
          "args":["--min-tumor-alt-freq"],
           "kwargs": {"help":"See below"}
        },
        {
          "args":["--min-tumor-alt-count"],
           "kwargs": {"help":"See below"}
        },
        {
          "args":["--min-tumor-depth"],
           "kwargs": {"help":"See below"}
        },
        {
          "args":["--chrom"],
           "kwargs": {"help":"See below"}
        },
        {
          "args":["--dry-run"],
           "kwargs": {"help":"See below","action":"store_true"}
        }
        
        
      ]
    },
    {   "cmd":"roh",
      "kwargs": {"help":'GEMINI Builtin analysis. See below ', "add_help":False,"parents":[analysis]},
      "arguments":[
       {
          "args":["--min-snps"],
           "kwargs": {"help":"See below"}
        },
        {
          "args":["--min-total-depth"],
           "kwargs": {"help":"See below"}
        },
        {
          "args":["--min-gt-depth"],
           "kwargs": {"help":"See below"}
        },
        {
          "args":["--min-size"],
           "kwargs": {"help":"See below"}
        },
        {
          "args":["--max-hets"],
           "kwargs": {"help":"See below"}
        },
        {
          "args":["--max-unknowns"],
           "kwargs": {"help":"See below"}
        },
        {
          "args":["-s"],
           "kwargs": {"help":"See below"}
        }
      ]
    },
    {   "cmd":"fusions",
      "kwargs": {"help":'GEMINI Builtin analysis. See below ', "add_help":False,"parents":[analysis]},
      "arguments":[
       {
          "args":["--in_cosmic_census"],
           "kwargs": {"help":"See below","action":"store_true"}
        },
        {
          "args":["--min-qual"],
           "kwargs": {"help":"See below"}
        },
        {
          "args":["--evidence_type"],
           "kwargs": {"help":"See below"}
        }
      ]
    },
    {   "cmd":"bcolz_index",
      "kwargs": {"help":'GEMINI Builtin analysis. See below ', "add_help":False,"parents":[parentparser]},
      "arguments":[
       {
          "args":["--cols"],
           "kwargs": {"help":"See below","action":"store_true"}
        }
      ]
    },
    
    {   "cmd":"update",
      "ignore":True,
      "kwargs": {"help":'GEMINI Builtin command. See below ', "add_help":False,"parents":[parentparser]},
      "arguments":[  
        {
          "args":["--devel"],
           "kwargs": {"help":"See below","action":"store_true"}
        },
        {
          "args":["--dataonly"],
           "kwargs": {"help":"See below","action":"store_true"}
        },
        {
          "args":["--nodata"],
           "kwargs": {"help":"See below","action":"store_true"}
        },
        {
          "args":["--extra"],
           "kwargs": {"help":"See below"}
        },
        {
          "args":["--tooldir"],
           "kwargs": {"help":"See below"}
        }

              ]
    },
    {   "cmd":"dump",
      "ignore":True,
      "kwargs": {"help":'GEMINI Builtin command. See below ', "add_help":False,"parents":[parentparser]},
      "arguments":[
        {
          "args":["--variants"],
           "kwargs": {"help":"See below","action":"store_true"}
        },
        {
          "args":["--genotypes"],
           "kwargs": {"help":"See below","action":"store_true"}
        },
        {
          "args":["--samples"],
           "kwargs": {"help":"See below","action":"store_true"}
        },
        {
          "args":["--header"],
           "kwargs": {"help":"See below","action":"store_true"}
        },
        {
          "args":["--sep"],
           "kwargs": {"help":"See below"}
        },
        {
          "args":["--tfam"],
           "kwargs": {"help":"See below","action":"store_true"}
        }

      ]
    },
    {   "cmd":"db_info",
      "ignore":True,
      "kwargs": {"help":'GEMINI Builtin command. See below ', "add_help":False,"parents":[parentparser]},
      "arguments":[]
    },
    {   "cmd":"qc",
      "ignore":True,
      "kwargs": {"help":'GEMINI Builtin command. See below ', "add_help":False,"parents":[parentparser]},
      "arguments":[{
          "args":["--mode"],
           "kwargs": {"help":"See below"}
        }, 
        {
          "args":["--chrom"],
           "kwargs": {"help":"See below"}
        }
        ]
    },
    {   "cmd":"examples",
      "ignore":True,
      "kwargs": {"help":'GEMINI Builtin command. See below ', "add_help":False,"parents":[parentparser]},
      "arguments":[
        
      ]

    }

  ]

  notwrapped = [
    {   "cmd":"browser",
      "kwargs": {"help":'GEMINI Builtin command. Is not wrapped ', "add_help":False,"parents":[parentparser]},
      "arguments":[]
    },
    {   "cmd":"load_chunk",
      "ignore":True,
      "kwargs": {"help":'GEMINI Builtin command. Is not wrapped ', "add_help":False,"parents":[parentparser]},
      "arguments":[]
    },
    {   "cmd":"merge_chunks",
      "ignore":True,
      "kwargs": {"help":'GEMINI Builtin command. Is not wrapped ', "add_help":False,"parents":[parentparser]},
      "arguments":[]
    }
  ]


  for v in builtincmds:
    v["parser"] =  addsubcommand(v)
    helpfunctions[v["cmd"]] = v["parser"]
  for v in notwrapped:
    v["parser"] =  addsubcommand(v)
    helpfunctions[v["cmd"]] = v["parser"]  


  if globalargs.help:
    print("""
--------------------------------------------------------------------------------------------

  The following help and commands are from the GEMINI wrapper.
  Parameters that can be recognised by either the wrapper or 
  by GEMINI, will also show specific GEMINI help below.

--------------------------------------------------------------------------------------------
     """)
  parentparser._positionals.title = "[sub-commands]"
  
  globalargs.permissions =  int(''.join( [ c if int(c) <=7 else "7" for c in str(globalargs.permissions)[-3::]]))

  wrappedcmds = [k for k in  helpfunctions.keys() if not k in uniquecmds]
  #print (parser_id.items())

  if not globalargs.help and globalargs.cmd != "" :
    subargs,remainingargs = subcmdparser.parse_known_args([globalargs.cmd] + remainingargs)
    if not (globalargs.cmd in uniquecmds or globalargs.cmd in wrappedcmds):
      remainingargs = [globalargs.cmd] + remainingargs
  # control help like this in order for the gemini help to also appear

  if globalargs.cmd == "" and not globalargs.saveconfig:
    globalargs.help = True

  

  if globalargs.help:
    if globalargs.cmd in helpfunctions.keys() :
      helpfunctions[globalargs.cmd].print_help()
    elif not globalargs.cmd == "":
      #unknowncmd.print_help()]
      pass
    else:

      subcmdparser.print_help() 
        
    pass
    print("\n\nAdditional Help")
    print("The following can be subsituted in paths for saved directories\n")
    print("{database_dir}     for the database directory set in the configuration file" )
    print("{log_dir}          for the log directory set in the configuration file" )
    print("{tmp_dir}          for the tmp directory set in the configuration file" )
    print("{consolidated_dir} for the consolidated directory set in the configuration file" )
    print("{prepped_dir}      for the prepped directory set in the configuration file" )
    print("{template_dir}     for the template directory set in the configuration file" )
    print("{user}             for the current username" )
    print("{user_dir}         for the current user's directory" )
    print("{id}               for the sample/family id (only in load command and depends on consolidation mode)" )
    print("{time}             for the time at which the script is run (only in load command)" )
  # logic to print message to standard output depending on verbose or silent flags
  debuglevel = 25
  if globalargs.debug:
      debuglevel = logging.DEBUG
      globalargs.verbose = True # ensure all verbose messages are given 
  elif globalargs.usage:
    debuglevel = 19
    globalargs.verbose = True # ensure all verbose messages are given 
  elif globalargs.verbose:
      # show everything above debug level if verbose
      debuglevel = logging.INFO
  elif globalargs.silent:
      # show everything above warning if user requested silent
      debuglevel = logging.ERROR

  
  
  
  
  # check if we are to create log file
  if globalargs.logfile and not globalargs.childprocess:
    global logpath
    #user wants to save log to file
    if globalargs.envreloaded is None or globalargs.envreloaded == "":
      logpath = os.path.join(config.get("paths",{}).get("log_dir",scriptpath),executiontime)
      # if there is already a log with this timestamp
      # append a number to the end 
      if len(glob.glob(logpath + ".log") ) > 0:
        logpath = logpath + "-" + str( len(glob.glob(logpath + ".log")) + 1)
      logpath += ".log"
    else:
      logpath = globalargs.envreloaded
    lfh = logging.handlers.WatchedFileHandler(logpath)
    if globalargs.debug:
      lfh.setLevel(logging.DEBUG)
    elif globalargs.usage:
      lfh.setLevel(11)
    else:
      lfh.setLevel(logging.INFO) # always minimum verbose level
    lfh.setFormatter(formatter)
    scriptlogger.addHandler(lfh)
  
  
  termlog.setLevel(debuglevel)


def byteliteraldecode(invar):
  return invar.decode('UTF-8')

def get_human_readable_size(num):
    exp_str = [ (0, 'B'), (10, 'KB'),(20, 'MB'),(30, 'GB'),(40, 'TB'), (50, 'PB'),]               
    i = 0
    rounded_val = round(float(num) / 2 ** exp_str[2][0], 2)
    return "{0:0.2f} {1}".format(rounded_val, exp_str[2][1])

# check to see if the run os has environment modules available
# if it does try to load the "conda" environment

def environment_load(themodule):
  """ load in a the environment, either as a linux environment module or by setting the PATH variable  """
  
  themodule = themodule.strip() # remove any extra whitespace
  if themodule == "": 
    return ""
  elif os.path.exists(parse_paths(themodule)):
    if os.path.isdir(parse_paths(themodule)):
      # if the passed module is exists in directory tree assume this is to be added to path variable
      os.environ["PATH"] += ":" + themodule
    return ""
    
  elif "MODULESHOME" in os.environ:     
    proc = ChildProcess('/usr/bin/modulecmd python load {0}'.format(themodule),communicate=True)
    exec (proc.stdout)
      # check to see if the module has been loaded
    if not themodule in os.environ["LOADEDMODULES"]:
      pass
      raise ScriptError("Configuration lists {0} to be loaded but this system does not have this path or module".format(themodule))
  else:
    raise ScriptError("Configuration lists {0} to be loaded but this system does not have this path or module".format(themodule))
  

def environment_test(cmd,notnone=False):
  """ Test whether the configuration values have been loaded """
  if " " in cmd:    
    proc = ChildProcess(cmd,communicate=True)
    if (proc.stderr== "") :
      return True
  else:
    if cmd and notnone:
      return True    
    if os.path.exists(parse_paths(cmd)):
      return True
    seen = set()
    for dir in os.environ.get("PATH", os.defpath).split(os.pathsep):
      normdir = os.path.normcase(dir)
      if normdir not in seen:
        seen.add(normdir)
        name = os.path.join(dir, cmd)
        if  os.path.exists(name) and os.access(name, os.F_OK | os.X_OK) and not os.path.isdir(name):
          return True
    
  return False

def environment_force(cmd,errmsg,notnone=False):
  """ Throws error if the cmd fails """
  if not environment_test(cmd,notnone):
    raise ScriptError("{0} not found. Please check configuration or read the documentation for further help.".format(errmsg))

def find_files(ext):
  infilelist = [] # dict for uncompressed and compressed files
  def expandglob(pathlist):    
    """ expands glob path """
    pathout = []
    base_dir = []
    for testpath in pathlist:
      if "{" in testpath:
        testpath = parse_paths(testpath)
      if os.path.exists(testpath):

        if os.path.isdir(testpath):
          # if a folder is passed we need to check for files in that folder
          newpathout,newbase_dir = expandglob(glob.glob(os.path.join(testpath,"*")))
          pathout += newpathout
          base_dir += newbase_dir
        elif len([p for p in ext if testpath.endswith(p)]) > 0 :
          
          pathout += [testpath]
          base_dir += [os.path.dirname(testpath)]
        
    return pathout,base_dir

  # determine passed folder 
  # only way to do it as bash expands shell wildcards 
  if type(subargs.source) is not list:
    subargs.source = [subargs.source]
  if len(subargs.source + remainingargs)== 0:
    return [],""
  
  #base_dir = os.path.dirname(subargs.source[0])
  
  return expandglob(subargs.source + remainingargs)

def parse_vcfs(inlist,base_dir="",filenames_only=False):
  """ opens the files and reads the column names """


  import string
  out = []
  # work through list and extract ids from vcf files
  # extract out all sample ids 
      
      
  for filepath in inlist:
    cvcf = VCFFile(filepath,raw=True)
    out.append(cvcf)
    message = "\n{0:<80}\n   Sample(s): ".format(cvcf.filepath)#,"Ready for GEMINI" if vcf.processed else "Pre processing required")
        # run through the ped of each file
    for id in cvcf.ped.keys():
      message += " {0} ".format(id)
    if globalargs.cmd == "id":
      scriptlogger.log(25,message)  
    else:
      scriptlogger.debug(message)  
    
  return out #set(out)

def get_seperator(format="TSV"):
  formats = {"TSV":"\t","SSV":" ","CSV":","}
  return formats.get(format," ")

def generate_blank_ped(ids,path,format="TSV"):
  """ generates a blank ped file at the given path (will overwrite) """
  sep = get_seperator(format)
    
  with open(path,"w") as thefile:
    thefile.write("""
#The standard PED file format must strictly contain the columns in the following order:
#FamilyID SampleID PaternalID MaternalID Sex Phenotype 
#Any columns to the right of this should contain additional genotype information
#The identified sample id are pre filled in column 1 and column 3
#This is important if a sample spans multiple vcf files
#FamilyID{0}SampleID{0}PaternalID{0}MaternalID{0}Sex{0}Phenotype{0}

""".format(sep))
    for id in ids:
      thefile.write("-9{0}{1}{0}-9{0}-9{0}-9{0}-9\n".format(sep,id))
  
def parse_ped(filepath):
  """ Opens the given filepath and parses the format. returns dict of samples """
  filepath = os.path.expanduser(filepath)
  filename = os.path.split(filepath)[1]
  sep = get_seperator(filename.split(os.path.extsep)[-1].upper())
  samples = {}
  sampleidlist = []
  warnmultiple = False
  invalidrows = 0 
  #try:
  if 1:
    with open( filepath ,"r") as pedfile:
      lines = 0 
      for line in pedfile:
        lines += 1 
        if line.replace(",","").replace(" ","").replace("\t","").strip() =="":
          continue
        sampleinfo = line.strip().split(sep)
        if line.startswith("#"):
          if len(samples.keys()) <2: 
            # comment is before first sample line
            samples["#HEADER#"] = [x.lower() for x in sampleinfo]
            # this will set mask to sample id if file prefix not supplied 
        else:      
          # this removes fileprefix from info as GEMINI doesn't need it!                    
          # 0 = family id
          # 1 = sample id
          
          if sampleinfo[1] in samples.keys():
            warnmultiple = True
          elif len(sampleinfo) < 6:
            invalidrows +=1
          else:
            if sampleinfo[0] in ["-1","-9"]: # if family missing, set to sampleid
              sampleinfo[0] = sampleinfo[1]
            samples[sampleinfo[1]] = sampleinfo

    if warnmultiple: 
      scriptlogger.warn("PED file contains samples with multiple definitions. The first instance in the file has been used")
    if invalidrows > 0:
      scriptlogger.warn("PED file contains invalid lines")
    samples.pop("#HEADER#","")          
    return samples
  #except:
  #  raise ScriptError("Supplied PED file does not match the correct file format. Please visit http://zzz.bwh.harvard.edu/plink/data.shtml for more information")
def check_directory(path,failto=""):
  if path == "":
    return failto
  
  path = os.path.abspath(os.path.expanduser(path))
  if not os.path.exists(path):
    os.makedirs(path)
    chmod(newfilepath)
    #os.mkdir(path)
  elif not os.path.isdir(path):
    return os.path.split(path)[0]
  return path

def chmod(path):
  """ sets the persmissions on supplied path as provided by cli or config file """
  try:
    ChildProcess("chmod {0} {1}".format(globalargs.permissions,path),communicate=True)
  except:
    pass

def initialprep(sampleid,vcfs):
  """ takes a sample id and vcf/vcf.gz that contain sample. Consolidates the data and returns vcf.gz using bcf tool """
  # from a list of filepaths
  tmpfilelist = []
  newfilepath_template = os.path.join(workingtmpfolder,"{id}_{type}.vcf.gz")
  scriptlogger.info("Initial preparation of VCF(s) containing sample: {0}".format(sampleid))
  
  for vcf in vcfs:
   
    # running bcftools view -O v will extract and decompress vcf 
    file_start_time = time.time()

    workingfilepath = vcf.filepath

    newfilepath = newfilepath_template.format(id=sampleid,type=vcf.varianttype)
    #zless = ChildProcess("zless {0}".format(workingfilepath))
    with ChildProcess.openfile(newfilepath,"w") as outfile: # # check if needs to keep outfile 

      #if vcf.filename.endswith(".vcf.gz") or len(vcf.ped) > 1 :        
      piper = ChildProcess("bcftools view -s {0} {1}".format(sampleid,workingfilepath)) # decompress and extract 
      sed = ChildProcess("sed 's/ID=AD,Number=1/ID=AD,Number=./'",stdin=piper.process.stdout)  # ensures the AD tag matches in all files    
      bgzip = ChildProcess("bgzip -c",stdin=sed.process.stdout,stdout=outfile)
      
      
      ChildProcess.poll_processes([piper,sed , bgzip,])

    chmod(newfilepath)
    piper.reportusage(11) # numbers assoicate program with info in logfile
    sed.reportusage(12)
    bgzip.reportusage(15)

    scriptlogger.log(19,"Time taken:{0}".format(time.time() - file_start_time))

    newvcf = VCFFile(newfilepath,ped=vcf.ped)
    newvcf.index()
    tmpfilelist.append(newvcf)
    
    
  
  return tmpfilelist
  
def consolidate_sample(sampleid,vcfs):
  scriptlogger.info("Consolidating VCFs containing sample: {0}".format(sampleid))
  tmpfilelist = []
  file_start_time = time.time()
  
  # generate list of vcf files with list comprehension
  newfilepath = os.path.join(workingtmpfolder,"{id}_{type}_consolidated.vcf.gz".format(id=sampleid,type='_'.join(set([vcf.varianttype for vcf in vcfs]))))
  bcftools = ChildProcess("bcftools concat -a -O z -o {0} {1}".format(newfilepath,' '.join([vcf.filepath for vcf in vcfs])))
  
  # close the stdout sequentially to run the file through the sequence
  # start the threads
  ChildProcess.poll_processes([bcftools])
  newvcf = VCFFile(newfilepath)
  chmod(newfilepath)
  bcftools.reportusage(11)
  scriptlogger.log(19,"Time taken:{0}".format(time.time() - file_start_time))

  for vcf in vcfs:
    newvcf.ped.update(vcf.ped)
    vcf.remove()  
  newvcf.index()
  tmpfilelist.append(newvcf)
  return tmpfilelist

def consolidate_group(vcfs,familyid=""):
  file_start_time = time.time()
  if familyid != "":
    scriptlogger.info("Consolidating VCFs containing family: {0}".format(familyid))
    newfilepath = os.path.join(workingtmpfolder,"{0}_consolidated.vcf.gz".format(familyid))
  else:
    scriptlogger.info("Consolidating VCFs containing all samples")
    newfilepath = os.path.join(workingtmpfolder,"all_consolidated.vcf.gz")
  tmpfilelist = []

  #bcftools merge -o 014_merged.vcf -O v 05_merged.vcf.gz 02_merged.vcf.gz
  if len(vcfs)> 1:
    bcftools = ChildProcess("bcftools merge -O z -o {0} {1}".format(newfilepath,' '.join([vcf.filepath for vcf in vcfs])))
  
  # close the stdout sequentially to run the file through the sequence
  # start the threads

    
  else:
    bcftools = ChildProcess("cp {0} {1}".format(' '.join([vcf.filepath for vcf in vcfs]),newfilepath))

  ChildProcess.poll_processes([bcftools])
  scriptlogger.log(11,bcftools.stderr)

  bcftools.reportusage(11)
  scriptlogger.log(19,"Time taken:{0}".format(time.time() - file_start_time))

  newvcf = VCFFile(newfilepath)
  chmod(newfilepath)
  for vcf in vcfs:
    newvcf.ped.update(vcf.ped)
    vcf.remove()

  
  newvcf.index()
  tmpfilelist.append(newvcf)
    
    
  return tmpfilelist 

def finalprep(id,vcfs):
  #filepath = os.path.expanduser(filepath)
  #outpath = os.path.expanduser(outpath)
  tmpfilelist = []
  for vcf in vcfs:
    #
    #  Pre process selected VCF files as part of the GEMINI workflow as described https://gemini.readthedocs.io/en/latest/index.html
    # VCF=/path/to/my.vcf
    #   NORMVCF=/path/to/my.norm.vcf.gz
    #  REF=/path/to/human.b37.fasta
    #   SNPEFFJAR=/path/to/snpEff.jar
    #   # decompose, normalize and annotate VCF with snpEff.
    #   # NOTE: can also swap snpEff with VEP
    #   zless $VCF \
    #    | sed 's/ID=AD,Number=./ID=AD,Number=R/' \
    #     | vt decompose -s - \
    #     | vt normalize -r $REF - \
    #     | java -Xmx4G -jar $SNPEFFJAR GRCh37.75 \
    #     | bgzip -c > $NORMVCF
    # 
    # Additionally this needs to create a ped file so that gemini knows what each sample is! 
    # https://gemini.readthedocs.io/en/latest/content/preping.html#describing-samples-with-a-ped-file

    
    newpath = os.path.join(config["paths"]["tmp_dir"],vcf.filename.replace(".vcf.gz","_prepped.vcf.gz"))
    

    previousfiles = glob.glob(newpath + "*")
    
    # clear any previous runs
    if len(previousfiles) > 0:
      for f in previousfiles:
        os.remove(f)
    
    scriptlogger.log(25,"Final preparation of {0}".format(vcf.filepath))
    file_start_time = time.time()
    with ChildProcess.openfile(newpath,"w") as outfile:# # check if needs to keep
    # prep preping programs
      bcftools = ChildProcess("bcftools view {0}".format(vcf.filepath)) # decompress
      sed = ChildProcess("sed 's/ID=AD,Number=./ID=AD,Number=R/'",stdin=bcftools.process.stdout) # should already be done
      vtd = ChildProcess("vt decompose -s -",stdin=sed.process.stdout)
      vtn = ChildProcess("vt normalize -r {0} - -".format(config["env"]["genome_reference"]),stdin=vtd.process.stdout)

      # TODO need to add VEP
      
      if  config["env"]["snpeff"] != "":
        annotator = ChildProcess("java -Xmx8G -jar {0} GRCh37.75".format(config["env"]["snpeff"]),stdin=vtn.process.stdout)
      else:
        annotator = ChildProcess(config["env"]["vep"],stdin=vtn.process.stdout)
      bgzip = ChildProcess("bgzip -c",stdin=annotator.process.stdout,stdout=outfile)
  

      ChildProcess.poll_processes([bcftools,sed,vtd,vtn,annotator,bgzip])


    scriptlogger.log(11,bcftools.stderr)
    scriptlogger.log(12,vtd.stderr)
    scriptlogger.log(12,vtn.stderr)
    scriptlogger.log(13,annotator.stderr)
    scriptlogger.log(14,bgzip.stderr)

    bcftools.reportusage(11)
    sed.reportusage(12)
    vtd.reportusage(13)
    vtn.reportusage(13)
    annotator.reportusage(14)
    bgzip.reportusage(15)

    scriptlogger.log(19,"Time taken:{0}".format(time.time() - file_start_time))
    vcf.remove()
    newfile = VCFFile(newpath,ped=vcf.ped)
    newfile.index()
    tmpfilelist.append(newfile)
  return tmpfilelist


def extractwrappedpcmds(filecmds,clicmds,args):
  # the default is the filecmds!  
  # extract the wrap commands and update from cli list
  filecmds = {k:v for k,v in filecmds.items() if k in args.keys()}
  filecmds.update({k:v for k,v in subargsdict.items() if k in args.keys() and ((not isinstance(v,bool) and v is not None) or args[k]["action"] == "store_" + str(v).lower() )  })
  if "header" in  filecmds.keys():
    filecmds["header"] = True # force header to true to ensure that 
  if "format" in  filecmds.keys():
    filecmds["format"] = None # force header to true to ensure that 
  return filecmds
def escape_chars(i):
  return "\"{0}\"".format(i) if any(c in i for c in [" ",">","<","|"]) else i 
def rebuildcli(args,info):
  """ Rebuilds cli parameters from list """
  cli = []
  for k,v in args.items():
    if isinstance(v,bool): 
      if (info["dest"][k]["action"] == "store_" +str(v).lower()):
        cli += [info["info"]["arguments"][          info["dest"][k]["index"]        ]["args"][0]]
    elif isinstance(v,list):
      cli += [info["info"]["arguments"][          info["dest"][k]["index"]        ]["args"][0] ] + [escape_chars(i) for i in v]
    elif v is not None:
      cli += [info["info"]["arguments"][          info["dest"][k]["index"]        ]["args"][0]]  +  [escape_chars(v)]
  return cli
def perform_wrap(cmds,communicate,output=None):
  """ send command to GEMINI """
  def process_gemini_output(line):
    line = line.strip()
    if isinstance(output,ReportWorkSheet):
      output.adddata(line.split("\t"))
    else:
      print(line,file=output)


  file_start_time = time.time()
  
  gemini = ChildProcess(["gemini"] + cmds,communicate=communicate,pollout=True)
  # poll process returns info on all process provided as dict...therefore extract out gemini info
  if not communicate: 
    if not globalargs.cmd in ["load"]:
      gemini.outcallback = process_gemini_output
    ChildProcess.poll_processes([gemini])
  gemini.reportusage(18)
  scriptlogger.log(19,"Time taken:{0}".format(time.time() - file_start_time))
  return gemini


def get_input(message):
  if sys.version_info[0] < 3:
    return raw_input(message)
  else:
    return input(message)

def parse_paths(inpath):
 
  try:
    return os.path.expanduser(inpath.format(database_dir=config["paths"]["database_dir"],
                                    tmp_dir=config["paths"]["tmp_dir"],
                                    log_dir=config["paths"]["log_dir"],
                                    consolidated_dir=config["paths"]["consolidated_dir"],
                                    prepped_dir   =config["paths"]["prepped_dir"]  ,         
                                    template_dir=config["paths"]["template_dir"],
                                    report_dir=config["paths"]["report_dir"],
                                    user_dir="~",
                                    user=os.environ["USER"]))
  except:
    raise ScriptError("Error parsing the following path: {0}".format(inpath))
#
#------------------------------------------------
#              Main Entry 
#------------------------------------------------
#

if __name__ == "__main__":   
  
  try:
    parse_envconfig()
    parse_args()

      
    if globalargs.envreloaded is None:
      scriptlogger.debug("Loaded config: {0}".format(config))
      scriptlogger.debug("Loaded args: {0}".format(sys.argv))

    # force value
    if not os.path.exists(cfgfilepath):
      print("No configuration file found at script location, assuming frist run.")
      print("Before you can use this script please read the readme.md")
      print("and then setup the configuration at {0}".format(cfgfilepath))
      generate_envconfig(cfgfilepath)
    elif globalargs.saveconfig:
      generate_envconfig(cfgfilepath)

    # before loading enviroments, validate some of the config
    #if subargs.consolidate and not subargs.consolidate in ["none","sample","family","all"]:
    #
    #   raise ScriptError("Config error with merge option")

    oldenvpath = os.environ["PATH"].split(":")

    environment_load(config["env"]["gemini"])
    environment_force("gemini","GEMINI")
    #environment_load(config["env"]["tabix"])
    environment_force("tabix","GEMINI")

    environment_load(config["env"]["anaconda"])  
    
    ## if anaconda environment has been supplied and is loaded
    # and anaconda is not the current environment
    # and psutilavail is false or config says to switch
    
    if globalargs.help:      
      if not globalargs.cmd in uniquecmds:
        
        # if the help parameter is passed we need to send this through to gemini, with the sub command (if any)
        print("""


----------------------------------------------------------------------------------------------------

    The following help and commands are from the GEMINI program.
    Using these parameters will override parameters set by within the wrapper script
    
    Please note...
    
    "usage: gemini" : Accessing gemini directly will be prevent wrapper script functions 
                      from working 
    
----------------------------------------------------------------------------------------------------

      """)
        proc = perform_wrap([globalargs.cmd,"-h"],True)
        print(proc.stdout)
        print(proc.stderr)
      exit()   
    
    
    if "ANACONDAHOME" in os.environ.keys() and not "anaconda" in os.environ["_"] and (not psutilavail or config["env"]["switch_to_anaconda"]):
      
      for i in range(len(sys.argv[1::])):
        sys.argv[1+ i] = escape_chars(sys.argv[1+ i])

      os.system(" ".join(["python"] + [sys.argv[0]] + ["--envreloaded","\"{0}\"".format(logpath) ] + sys.argv[1::]))
      exit()
    
    if not psutilavail:
      scriptlogger.info("Memory and CPU not available. Please either use an anaconda environment or install psutil manually or via pip")
    # load the rest of the environment for loading data
    if globalargs.cmd  in ["load","prep"] :      
      # load and check vt
      environment_load(config["env"]["vt"])
      environment_force("vt","vt")

      # check snpEff and throw if error
      # TODO add vep 
      snpeffpath = config["env"]["snpeff"]
      veppath = parse_paths(config["env"]["vep"])
      
      if environment_test(snpeffpath):
        environment_force("java -Xmx8G -jar {0} -version".format(parse_paths(config["env"]["snpeff"])),"snpEff")
      elif environment_test(veppath):
        config["env"]["snpeff"] = ""
        raise ScriptError("VEP is currently not supported, please use snpEff")
        environment_force(veppath,"VEP")
        
      else:
        raise ScriptError("Program needs either snpEFF or VEP to provide functional annotation")
      # check reference and throw if error
      environment_force(config["env"]["genome_reference"],"Genome reference")#"~/exome_seq_platform/bin/HG19/hg19.inuse.fasta")
      # check bcftools and throw if error
      environment_load(config["env"]["bcftools"])

      environment_force("bcftools --version-only","bcftools")
    
    ### Do some directory checks before we start
    if globalargs.cmd in ["load","id","cconsolidate","prep"]:
      # check if the user has passed -db parameter
      subargs.output = os.path.join(parse_paths(os.path.split(subargs.output)[0]),os.path.split(subargs.output)[1])
      # check if passed element is a file or directory
      check_directory(os.path.split(subargs.output)[0])

    if globalargs.cmd in ["load"]:
      #check user has supplied a database
      if os.path.isdir(subargs.output):
        subargs.output = os.path.join(subargs.output,"")

      if os.path.split(subargs.output)[1] == "":
        raise ScriptError("No database specified")
      


    # ensure the tmp folder exists 
    workingtmpfolder = check_directory(config["paths"]["tmp_dir"],os.path.join(scriptpath,"tmp"))
    
    if globalargs.cmd in ["id","consolidate","prep","load"]:
      
      # first get all vcf files

      infilelist,base_dir = find_files((".vcf",".vcf.gz"))
      if len(infilelist) == 0:
        raise ScriptError("No vcf files found")
      infilelist = parse_vcfs(infilelist)

    # ID cmd for generating sample ID list  
        
    if globalargs.cmd == "id" and subargs.output != "":
      generate_blank_ped(set([id for vcf in infilelist for id in vcf.ped.keys()]),subargs.output ,"CSV")

    if globalargs.cmd in ["load","consolidate","prep"] and not subargs.ped:       
      subargs.ped = os.path.join(workingtmpfolder,"tmpped.csv")
      generate_blank_ped(set([id for vcf in infilelist for id in vcf.ped.keys()]),subargs.ped,"CSV")
      scriptlogger.log(25,"No PED file provided, temporary PED file (as CSV) created here: {0}".format(subargs.ped))
      scriptlogger.log(25,"It is strongly advised you edit this file to ensure the accurate loading into GEMINI")
      
      if not globalargs.silent:
        input = get_input("Would you like to edit this PED file now? (y/n) ").lower()
      else:
        input = "n"
      if input =="y":
        scriptlogger.debug("User selected to edit file")
        os.system("nano " + subargs.ped)
        
        # cli creation
        pass
    realstart = time.time()
    currentsteptime = time.time()

    #scriptlogger.info("Current memory by script: {0}MB".format(round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / rusage_denom,2)))
    if globalargs.cmd in ["load","consolidate","prep"]:
      # parse the ped and compare against selected files
      peds = parse_ped(subargs.ped)

      scriptlogger.log(25,"Please be aware that the time taken to prep files depends upon how many were selected and the options provided")
      

      # combine ped with selected files      
      for vcf in infilelist:
        vcf.ped = { k : peds.get(k,pedvalue) for k,pedvalue in vcf.ped.items() if k in vcf.ped.keys()}
        samplecount = len(vcf.ped.keys())

        for id,info in vcf.ped.items():
          if not id in samplelist:
            samplelist[id] = [vcf]        
          else:
            samplelist[id] += [vcf]
      
      # pre processing steps
      total = sum([len(samplelist[id]) for id in samplelist])
      scriptlogger.info("Sample data loaded in {0}".format(time.time() - currentsteptime ) )
      currentsteptime = time.time()
      scriptlogger.log(25,"Initial preparation of {0} sample(s) across {1} file(s) ".format(len(samplelist),total))

      for id in samplelist.keys():
        # 1) foramt AD tag (required for gemini) + also ensure proper consolidation
        samplelist[id] = initialprep(id,samplelist[id])

        if not subargs.consolidate == "none" and len(samplelist[id]) > 1:          
          # only consolidate if sample appears across multiple files...no point wasting time and effort
          samplelist[id] = consolidate_sample(id,samplelist[id])
      
      # 2) consolidate groups if requested
      if len(samplelist) > 1:
        if subargs.consolidate == "family":
          # validids should now contains family unit
          for famid in set([ vcfs[0].ped[id][0] for id,vcfs in samplelist.items() if vcfs[0].ped[id][0] != vcfs[0].ped[id][1]]):
            currentvcfs = []
            for sid in [ sid for sid,vcfs in samplelist.items() if vcfs[0].ped.get(sid,["0"])[0] == famid ]:
              currentvcfs.extend(samplelist.pop(sid))
            
            if len(currentvcfs) > 0:
              samplelist[famid] = consolidate_group(currentvcfs,familyid=famid)
        elif subargs.consolidate == "all":
          currentvcfs = []
          for id in [id for id in samplelist.keys()]:
            currentvcfs.extend(samplelist.pop(id))
          samplelist["all"] = consolidate_group(currentvcfs)
            #consolidate_group()

        
      if not subargs.consolidate == "none":
        VCFFile.movelist("consolidated",config["paths"]["consolidated_dir"],samplelist)
        scriptlogger.log(25,"Consolidated file(s) to {0} file(s)".format(sum([len(samplelist[id]) for id in samplelist])))

            
    if globalargs.cmd in ["load","prep"]:  
      # 3) decompose > normalize > annotate > index
      total = sum([len(samplelist[id]) for id in samplelist ])
      complete = 0
      if total > 1:
        scriptlogger.log(25,"{0} files to be prepped. Please be patient".format(total))
      for id in samplelist:
        samplelist[id] = finalprep(id,samplelist[id])
        complete += len(samplelist[id])
        if total > 1:
          scriptlogger.log(25,"{0} / {1} files prepped. Please be patient".format(complete,total))
    
      
      # add parents to list ped list
      for fid,vcfs in samplelist.items():
        for vcf in vcfs:
          for pid,mid in [ [values[2],values[3]] for id,values in vcf.ped.items() if values[2] in peds.keys() or values[3] in peds.keys() ]:
            vcf.ped[pid] = peds[pid]
            vcf.ped[mid] = peds[mid]

      #VCFFile.movelist("prepped",config["paths"]["prepped_dir"],samplelist)
      
    if globalargs.cmd in ["load"]:
      # check if the user has passed -db parameter
      #TODO expand replaceable texts to other directories and ids
      base_dir = parse_paths(os.path.split(subargs.output)[0])
      filename = os.path.split(subargs.output)[1]

      dbname_template = []

      if os.path.exists(subargs.output):
        passedfile =  os.path.isfile(subargs.output)    
      else:
        passedfile = ("." in filename)

      # build database filepath
      dbname_template = list(set(dbname_template))
      if len(samplelist)>1 and not "{id}" in dbname_template :        
        dbname_template += ["{id}"]

      if subargs.consolidate == "none" and not "{type}" in dbname_template :
        dbname_template += ["{type}"]
      
      if passedfile:
        dbname_template += [filename]
      else:
        dbname_template += [".db"]
      dbname_template = "_".join(dbname_template)
      dbname_template = dbname_template.replace("_.",".")
      
      
      pedfilepath = os.path.join(workingtmpfolder,"load.ped")
      scriptlogger.log(25,"Beginning load into GEMINI. The time taken depends upon how many files/variants are processed by the \"--cores\" parameter.")
      total = total = sum([len(samplelist[id]) for id in samplelist ])
      complete = 0
      if total > 1:
        scriptlogger.log(25,"{0} files to be loaded. Please be patient".format(total))
      for id in samplelist:
        for vcf in samplelist[id]:
          complete += 1
          # load vcf into a gemini database...
          #
          # dbname
          dbpath = os.path.join(base_dir,dbname_template.format(id=id,time=executiontime,type=vcf.varianttype))

          # generate ped for this vcf
          with open(pedfilepath,"w") as pedfile:
            
            for l in vcf.ped.values():
              pedfile.write("{0}\n".format(" ".join(l) ))
          ChildProcess.currentopenfile.append(dbpath)
          # error contains pids and how many variants processed....useful stats?
          result = perform_wrap(["load","-v",vcf.filepath,"--cores",str(subargs.cores),"-p",pedfilepath,"-t","snpEff"]+remainingargs + [dbpath],False)
          scriptlogger.log(18,result.stderr.replace("""/opt/gridware3/el7/apps/gcc/python-packages/gemini/share/anaconda/lib/python2.7/site-packages/gemini/config.py:61: YAMLLoadWarning: calling yaml.load() without Loader=... is deprecated, as the default Loader is unsafe. Please read https://msg.pyyaml.org/load for full details.
  config = yaml.load(in_handle)""",""))
          scriptlogger.log(25,"Finished loading GEMINI database : " + dbpath)
          chmod(dbpath)
          while(ChildProcess.currentopenfile.count(dbpath)):
            ChildProcess.currentopenfile.remove(dbpath)
          
          if total > 1:
            scriptlogger.log(25,"{0} / {1} files loaded. Please be patient".format(complete,total))
    
      



    #
    #   ANALYSIS COMMANDS
    #

    if globalargs.cmd in [cmd["cmd"] for cmd in builtincmds] or globalargs.cmd in ["analysis"]:
      
     
      
      # this should handle all other commands
      infilelist,base_dir = find_files((".db"))
      if len(infilelist) == 0:
        raise ScriptError("No database files found")
      

      if subargs.listheading :
        try:
          import json
        except:
          raise ScriptError("Python Package Required: json")
        dbpath = infilelist[0]

        #result = perform_wrap(["query","-q","\"select * from {0} limit 1\"".format(subargs.listheading),"--header"]+remainingargs + [dbpath],True)     
        result = perform_wrap(["query","-q","\"SELECT tbl_name,sql FROM sqlite_master where type='table' \"","--header","--format","json"]+remainingargs + [dbpath],True)     
        print("\nShowing the all tables in " + dbpath + "\n")
        
        columns = 5
        for r in result.stdout.split("\n"):
          if r == "":
             continue
          r = json.loads(r)
          print("\nShowing the headings for table : " + r["tbl_name"] + "")
          headings = [h.strip().split()[0] for h in r["sql"].split("\n")[1:-1]]
          for i in range(0,len(headings),columns):
            
            l = "\t".join([ "{0:<25}".format(h.replace("\"","")) if h not in ["PRIMARY","CHECK"] else "" for h in headings[i:i+columns] ])
            if l.strip() == "":
              break
            print(l)
# NOT USING THIS FORMAT NOW....keep here in case
 #       print("\nTo use these with the query command file, please prefix the table name eg. table.column")
 #       print("\nSections in query file are:\n")
 #       print("[sql]")
 #       print("sql statement\n")
 #       print("[gemini]")
 #       print("GEMINI sub command, to allow use of built in analysis\n")
 #       print("[columns]")
 #       print("one column per line\n")
 #       print("[filter]")
 #       print("Lines are automatically combined with \"AND\" unless \"OR\" exists at beginning or end of line\n ")
 #       print("[extra]")
 #       print("SQL structure definitions such as \"limit\", \"order by\" etc")
#        print("for complex queries it is advised to use the sql section\n ")
#        print("[CLI]")
#        print("Additional command line arguments GEMINI can accept\n")
#        print("Please look at https://gemini.readthedocs.io/en/latest/content/querying.html for advanced options")
      else:
         # first built the arg list for the builtin commands
        for v in builtincmds:
          builtinargdict[v["cmd"]] = {"dest":{
            arg["args"][0].replace("-","_").strip("_") if "dest" not in arg["kwargs"] else arg["kwargs"]["dest"]: 
              {"action": arg["kwargs"]["action"] if "action" in arg["kwargs"].keys() else "",
              "index" : i }
                for i,arg in enumerate(v["arguments"])
          },
          "info":v}  
        

        if subargs.output == "":
          scriptlogger.log(25,"No output file selected, will output to stdout. (All log messages print on stderr)")
          subargs.format = "TSV"
          outpipe = sys.stdout 
        else:
          check_directory(os.path.split(parse_paths(subargs.output))[0])

          opath = os.path.split(subargs.output)
          fileext = opath[1].split(".")[-1]
          if fileext.lower() != "xlsx":
            subargs.format= "tsv"
            outpipe = open(subargs.output,"w")
          else:
            subargs.format= "xlsx"
            
        subargs.format = subargs.format.lower()
        if not subargs.output == "" and subargs.format == "xlsx":
          try:
            import xlsxwriter
          except:
            raise ScriptError("Python Package Required: xlsxwriter")
          usingxlsx = True
          output = xlsxwriter.Workbook(subargs.output)
          infows = output.add_worksheet("Information")
        else:
          usingxlsx = False
          
          infows = {} # ensure following commands don't cause errors
          
        
        
        subargs.format = subargs.format.lower()
        orgcli = " ".join([ p if not " " in p else "\"{0}\"".format(p) for p in sys.argv if not p == "--envreloaded"])
        if usingxlsx:
          infows.write(0,0,"Information about query")
          infows.write(1,0,"CLI Parameters" )
          infows.write(1,1, orgcli)
        else:
          print("""
Information about query 
 - CLI Parameters       {0}""".format(orgcli  ),file=outpipe)

        subargsdict = vars(subargs) # convert object into dictionary 
        wrapped = [{"cmd": globalargs.cmd,"args":extractwrappedpcmds({},subargsdict,builtinargdict[globalargs.cmd]["dest"])}]
        if not subargs.queryfile == "":


          subargs.queryfile = parse_paths(subargs.queryfile)
          if usingxlsx:
            infows.write(2,0, "Passed QueryFile")
            infows.write(2,1,subargs.queryfile)
          else:
            print(""" - Passed QueryFile     {0}""".format(subargs.queryfile),file=outpipe)
          
        
          querybricks = ""
          sqlbrick = False
            # ensure all lines do not have a space at beginning to ensure accurate parse
          queryconfig = {}
          def parse_individual_section(input):
            r = parse_config(querybricks,True)
            if not len(r) == 1:
              raise ScriptError("Issue with template file")
            cmd = list(r.keys())[0]
            if cmd in queryconfig.keys():
              queryconfig[cmd].append(r[cmd])
            else:
              queryconfig[cmd] = [r[cmd]]  
          try:
          
            with open(subargs.queryfile,"r") as queryfile:
              for l in queryfile:
                if l.strip().startswith("[") and not querybricks.strip() == "":
                  parse_individual_section(querybricks)
                  querybricks = "" 
                querybricks += l.strip() + "\n"  
                if l.strip().startswith("["):
                  querybricks += "falsecmd" + "\n" #  this ensure all sections are parsed correctly

            parse_individual_section(querybricks)
          except ScriptError as se:
            raise ScriptError(se)  #pass it up the chain
          

          querylinks = {"variantsgene_summary":"variants.gene=gene_summary.gene",
                        "variantsvariant_impacts":"variants.variant_id=variant_impacts=variant_id",
                        "variant_impactsgene_summary":"variant_impacts.gene=gene_summary.gene",
                        "variant_impactsgene_detailed": "variant_impacts.gene=gene_detailed.gene",
                        "variantsgene_detailed":"variants.gene=gene_detailed.gene"}
          query = {"tables":["variants"]}
          wrapped = []
          destlist = {}
    
          def addtowrap(cmd,filecmds):
            # add "to be wrapped" list
            for cmdinfo in filecmds:
              # protect the command line from extra quotes / out of place special characters
              cmdinfo = shlex.split(" ".join([ k + "=" + v.strip() if v is not None else k for k,v in cmdinfo.items() if not k=="falsecmd" ] )) 

              args,unknown = builtinargdict[cmd]["info"]["parser"].parse_known_args(["--template","faketemplate"] + cmdinfo + [" fakesource "])
              pcmds = extractwrappedpcmds(vars(args),subargsdict,builtinargdict[cmd]["dest"])
              wrapped.append({"cmd": cmd,"args":pcmds})


          if globalargs.cmd == "analysis": # do some thing different for the analysis command
            # here we check all commands in the builtincommands and parse apporpriately 
            for k in queryconfig:
              if k in builtinargdict.keys():
                addtowrap(k,queryconfig[k])
                # command is valid so 
          elif globalargs.cmd in queryconfig.keys():
            # built command list from file
            addtowrap(globalargs.cmd,queryconfig[globalargs.cmd])
           # filecm
          elif len(builtinargdict[globalargs.cmd]["dest"]) > 0 :
            scriptlogger.log(25,"File contains no parameters that match this command.")
          else:
            wrapped.append({"cmd": globalargs.cmd,"args":extractwrappedpcmds({},subargsdict,builtinargdict[globalargs.cmd]["dest"])})


          # this code was suppose to make things simpler for the end user..if time implement
          """
          # check for gemini analysis command
          if len(queryconfig.get("gemini",{}).keys()) > 0 and  globalargs.cmd == "query":
            # take the first line as the command
            globalargs.cmd = list(queryconfig["gemini"].keys())[0].lower().strip()          
            allowwrap = (globalargs.cmd != "query") # default to allowwrap is based upon whether the query is called or not
          if len(queryconfig.get("sql",{})) > 0 and  globalargs.cmd == "query":
              subargs.sql = ' '.join( [ k + "=" + v.strip() if v is not None else k for k,v in queryconfig["sql"].items()  ]  )
              allowwrap = True  
          else:

            subargs.sql = ["select"]
            if not globalargs.cmd == "query" and not subargs.columns == "":
              print("loaded from cli")
              query["columns"] = "--columns \"" + subargs.columns + "\""
            if len(queryconfig.get("columns",{}).keys()) > 0:             
              query["columns"] = [",".join(queryconfig.get("columns",{}).keys())]
              subargs.sql += query["columns"]
              if globalargs.cmd == "query":
                query["tables"] = set([ t.split(".")[0].lower() for t in queryconfig.get("columns",{}).keys()])
                if len(query["tables"]) == 0:
                  query["tables"] = ["variants"]
                subargs.sql += ["from"] + [",".join(query["tables"])]
            else:
              subargs.sql += ["* from variants"]
              query["columns"] = ["chrom","start","end"] # default for built in tools
            
            if not globalargs.cmd == "query" and not subargs.filter == "":
              query["columns"] = "--filter \"" + subargs.columns + "\""
            if len(queryconfig.get("filter",{}).keys()) > 0 or len(query["tables"]) > 1:
              subargs.sql += ["where"]
              query["filter"]  =  [" and ".join(
                              # parse user supplied parameters
                              [ k + "=" + v.strip() if v is not None else k for k,v in queryconfig.get("filter",{}).items()  ] +
                              # parse table joins 
                              [ querylinks[l] for l in [i+j for i in query["tables"] for j in query["tables"]] if l in querylinks]
                              )]
              subargs.sql += query["filter"]                 
            if len(queryconfig.get("extra",{}).keys()) > 0:
              subargs.sql +=[ k.strip() for k in queryconfig["extra"].keys()  ]
            subargs.sql = " ".join(subargs.sql)
            for d in [["and or","or"],["or and","or"]]:
              subargs.sql = subargs.sql.replace(d[0],d[1])
          """
         
        if usingxlsx:
          
          infows.write(7,0,"Database Index")
          infows.write(7,1, "Command")
          infows.write(7,2, "Command Index")
          infows.write(7,3,  "Database Path")
          infows.write(7,4, "Command Parameters")
        else:
          print("\n ======================================================================\n",file=outpipe)
        
        
        rowptr = 1
        scriptlogger.log(25,"There are {0} command(s) to be run on {1} database(s)".format(len(wrapped),len(infilelist)))   
        for dbindex,dbpath in enumerate(infilelist):
          cmdcounts = {}
          scriptlogger.log(25,"Running commands(s) on " + dbpath)
          for  wrap in  wrapped:
            if wrap["cmd"] in cmdcounts.keys():
              cmdcounts[wrap["cmd"]] += 1
            else:
              cmdcounts[wrap["cmd"]] = 1 
            if isinstance(wrap["args"],dict):
              wrap["args"]  = rebuildcli(wrap["args"],builtinargdict[wrap["cmd"]])
              wrap["argstring"] = " ".join(wrap["args"] )
              scriptlogger.info("Final parameters to be passed to GEMINI for command {0}: \n{1}".format(wrap["cmd"],wrap["argstring"]))    

            
            if usingxlsx:
            
              ws  = ReportWorkSheet(output.add_worksheet("DB" + str(dbindex + 1) + " - " + wrap["cmd"] + str(cmdcounts[wrap["cmd"]])))
              infows.write(7+rowptr,0,(dbindex+1))
              infows.write(7+rowptr,1,wrap["cmd"])
              infows.write(7+rowptr,2,cmdcounts[wrap["cmd"]])
              infows.write(7+rowptr,3,dbpath)
              infows.write(7+rowptr,4,wrap["cmd"] + " " + wrap["argstring"])
              rowptr +=  1
            else:
              ws = outpipe
              print("""
Database  : {0}
Command   : {1}             
              """.format(dbpath,wrap["cmd"] + " " + wrap["argstring"]),file=outpipe)

                      
            result = perform_wrap([ wrap["cmd"]]+wrap["args"] + [dbpath],subargs.output == "",output=ws)     
            
            if subargs.output == "":
              print(result.stdout.replace("\t",get_seperator(subargs.format)))

            if result.stderr != "":
              scriptlogger.log(18,result.stderr.replace("""/opt/gridware3/el7/apps/gcc/python-packages/gemini/share/anaconda/lib/python2.7/site-packages/gemini/config.py:61: YAMLLoadWarning: calling yaml.load() without Loader=... is deprecated, as the default Loader is unsafe. Please read https://msg.pyyaml.org/load for full details.
  config = yaml.load(in_handle)""",""))

            if not usingxlsx:
              print("\n ======================================================================\n",file=outpipe)


        if usingxlsx:
          
          output.close()
        
        if not subargs.output == "":
          scriptlogger.log(25,"Report saved : " + subargs.output)
        #else:
        #  scriptlogger.log(25,"No valid query passed to script")

  except KeyboardInterrupt:
    pass
  except ScriptError as ex:
    scriptlogger.fatal(ex)
  except IOError as ex:
    scriptlogger.fatal(ex)
  except OSError as ex:
    scriptlogger.fatal(ex)
  except Exception as ex:
    str(ex)
    import traceback
    traceback.print_exc()
  finally:
    # clean up temp files
    for id,vcfs in samplelist.items():
      for vcf in vcfs:
        vcf.remove()
        pass
    

    # in case the program was exited mid file write
    for p in ChildProcess.currentopenfile:
      if os.path.exists(p):
        pass
        os.remove(p)
      

else:
  print("GEMINI Wrapper script requries to be run as main script to ensure correct parsing of passed parameters")



