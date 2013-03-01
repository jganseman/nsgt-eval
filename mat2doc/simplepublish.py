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

if 'verify' in todo:
    printdoc.printdoc(projectname,'verify')

# Release for other developers to download
if 'develmat' in todo:
    #printdoc.git_stageexport(project['dir'],project['mat'])
    os.system('rm -rf '+project['mat']+'*')
    os.system('svn export --force '+project['dir']+' '+project['mat'])
    printdoc.printdoc(projectname,'mat')

    fname=filesdir+projectname+'-devel-'+versionstring

    # Create the Unix src package
    os.system('tar zcvf '+fname+'.tgz '+projectname+'/')

    # Create the Windows src package
    os.system('rm '+fname+'.zip')
    printdoc.unix2dos(filesdir+projectname)
    os.system('zip -r '+fname+'.zip '+projectname+'/')

# Release for users to download
if 'releasemat' in todo:
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
    
if 'pdf' in todo:
    printdoc.printdoc(projectname,'tex','rebuild')


if 'php'==todo:
    printdoc.printdoc(projectname,'php')

    s='rsync -av '+project['php']+'/ '+host+':'+www+''
    print s
    os.system(s)    

if 'phplocal'==todo:
    printdoc.printdoc(projectname,'phplocal')


if 'phprebuild'==todo:
    printdoc.printdoc(projectname,'php','rebuild')
    s='rsync -av '+project['php']+'/ '+host+':'+www+''
    os.system(s)    

if 'sendserver'==todo:
    s='rsync -av '+project['php']+'/.. '+host+':'+www+''
    os.system(s)  


        
    #os.system('rsync -av '+notehtml+' '+host+':'+noteswww);

