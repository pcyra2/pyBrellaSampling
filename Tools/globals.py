# global WorkDir    
# global verbosity
# global DryRun
# global parmfile
verbosity=0
WorkDir="./"
DryRun=True
parmfile="complex.parm7"

def init(v=0, wd="./", dr=True, parm="complex.parm7"):
    global WorkDir    
    global verbosity
    global DryRun
    global parmfile
    verbosity=v
    WorkDir=wd
    DryRun=dr
    parmfile=parm
    