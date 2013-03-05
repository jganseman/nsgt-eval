#!/usr/bin/python

import sys,os
cwd=os.getcwd()+'/'

# ------- Configuration parameters -------------

projectname='nsg'

# Configure HTML placement at remote server
host='nholighaus,nsg@web.sourceforge.net'
www='/home/project-web/nsg/htdocs'


# ------- do not edit below this line ----------

# Import the localconf file, it should be place in the same directory
# as this function is called from.

sys.path.append(cwd)

import localconf

sys.path.append(localconf.mat2docdir)

# Get the data from localconf
project=eval('localconf.'+projectname)
conffile=project['dir']+'/mat2doc/mat2docconf.py'
filesdir=localconf.filesdir

f=file(project['dir']+projectname+'_version')
versionstring=f.read()[:-1]
f.close()

import printdoc, notes
#from datesuffix import *


todo=sys.argv[1]


# Release for users to download
if 'releasemat' in todo:
    #printdoc.git_repoexport(project['dir'],'master',projectname,filesdir)
    os.system('rm -rf '+project['mat']+'*')
    os.system('cp -r '+project['dir']+' '+project['mat'])
    os.system('rm -r -f '+project['mat']+'.git')

    # Remove unwanted files
    os.system('rm -rf '+project['mat']+'mat2doc')
    os.system('rm -rf '+project['mat']+'test_bench')
    #os.system('rm -rf '+project['mat']+'timing')

    printdoc.printdoc(projectname,'mat')
    
    fname=filesdir+projectname+'-'+versionstring
    print fname

    # Create the Unix src package
    os.system('tar zcvf '+fname+'.tgz '+projectname+'/')

    # Create the Windows src package
    os.system('rm '+fname+'.zip')
    printdoc.unix2dos(filesdir+projectname)
    os.system('zip -r '+fname+'.zip '+projectname+'/')
    
    
if 'pdf' in todo:
    printdoc.printdoc(projectname,'tex')


if 'php'==todo:
    printdoc.printdoc(projectname,'php')

    s='rsync -av '+project['php']+'/ '+host+':'+www+''
    print s
    os.system(s)    

if 'html' in todo:
    if 'rebuild' in todo:
        printdoc.printdoc(projectname,'html','rebuild')
    else:
        printdoc.printdoc(projectname,'html')

if 'sendserver' in todo:
    s='rsync -av '+project['php']+'/ '+host+':'+www+''
    os.system(s)  

if 'package' in todo:
    # 1) release mat

    #printdoc.git_repoexport(project['dir'],'master',projectname,filesdir)
    os.system('rm -rf '+project['mat']+'*')
    os.system('svn export --force '+project['dir']+' '+project['mat'])
	
    # Remove unwanted files
    os.system('rm -rf '+project['mat']+'mat2doc')
    os.system('rm -rf '+project['mat']+'test_bench')
    #os.system('rm -rf '+project['mat']+'timing')

    printdoc.printdoc(projectname,'mat')
    
    fname=filesdir+projectname+'-'+versionstring

    # Create the Unix src package
    os.system('tar zcvf '+fname+'.tgz '+projectname+'/')

    # Create the Windows src package
    os.system('rm '+fname+'.zip')
    printdoc.unix2dos(filesdir+projectname)
    os.system('zip -r '+fname+'.zip '+projectname+'/')

    # 2) print html
    printdoc.printdoc(projectname,'html')
    
    # 3) print tex
    printdoc.printdoc(projectname,'tex')

    # 4) create a package



        
    #os.system('rsync -av '+notehtml+' '+host+':'+noteswww);

