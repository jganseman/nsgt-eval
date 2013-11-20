#!/usr/bin/env python

import sys,os
from mat2doc import *

plotengine=sys.argv[1]

cwd=os.getcwd()+'/'

#newlocals={}	
#execfile(conffile,globals(),newlocals)

# Get global configuration
#globalconf=newlocals['conf']

# Assemble configuration object
conf=ConfType()
	
conf.g=ConfType()
conf.g.root=cwd

# Global
#conf.g.plotengine=plotengine

# Target, none in this case
conf.t={}	

fname=sys.argv[2]

P=FunPrinter(conf,fname)

if plotengine=='octave':
   execplot('octave','octave',P.parsed['code'],'dummy',cwd,fname,sys.argv[3])
else:
   execplot('matlab','matlab',P.parsed['code'],'dummy',cwd,fname,sys.argv[3])


#plotexec(conf.g,conf.g.octexec,P.parsed['code'],conf.g.root,conf.g.workdir,fname,sys.argv[2])

