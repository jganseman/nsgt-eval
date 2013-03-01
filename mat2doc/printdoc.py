#!/bin/env python

import sys, os, os.path, shutil, re, commands
from mat2doc import *
from notes import *

# Append the current path so that we can load localconf
cwd=os.getcwd()+os.sep
sys.path.append(cwd)

import localconf

# rm -Rf
# Does not remove the directory itself
def rmrf(s):
    for root, dirs, files in os.walk(s, topdown=False):
        for name in files:
            os.remove(os.path.join(root, name))
        for name in dirs:
            os.rmdir(os.path.join(root, name))


# Make sure that a Git repo is on a specific branch
def assert_git_on_branch(repo,branchname):
    buf=commands.getoutput('git --git-dir='+repo+'.git branch | grep "*"')
    if not (buf.strip()[2:]==branchname):
        print 'Git repo',repo,'is not on branch "'+branchname+'". Stopping.'
        sys.exit(0)

def autostage(repo):
    #s = 'git --git-dir='+repo+'.git/ add -u'
    #os.system(s)
    pass

def git_stageexport(repo,outdir):
    rmrf(outdir)
    os.system('git --git-dir='+repo+'/.git/ checkout-index --prefix='+outdir+' -a')

def git_repoexport(repo,branchname,packagename,baseoutdir):
    rmrf(baseoutdir+packagename)
    s='git --git-dir='+repo+'/.git/ archive --prefix='+packagename+'/ '+branchname+' | (cd '+baseoutdir+' && tar xf -)'
    print s
    os.system(s)

def dos2unix(path):
    for root, dirs, files in os.walk(path, topdown=False):
        for name in files:
            name=os.path.join(root, name)
            if 'ASCII' in commands.getoutput('file '+name):
                os.system('dos2unix '+name)

def unix2dos(path):
    for root, dirs, files in os.walk(path, topdown=False):
        for name in files:
            name=os.path.join(root, name)
            if 'ASCII' in commands.getoutput('file '+name):
                os.system('unix2dos '+name)                


def printdoc(projectname,targetname,rebuildmode='auto'):

    project=eval('localconf.'+projectname)
    conffile=project['dir']+'/mat2doc/mat2docconf.py'

    newlocals={}
    execfile(conffile,globals(),newlocals)

    # Get global configuration
    globalconf=newlocals['conf']

    targetconf=newlocals[targetname]
    
    conf=ConfType()

    # Global
    conf.g=globalconf
    conf.g.root=project['dir']+os.sep
    conf.g.plotengine=localconf.plotengine
    conf.g.tmpdir=localconf.tmpdir+os.sep
    conf.g.bibfile=conf.g.root+'mat2doc/project'

    # Target
    conf.t=targetconf
    conf.t.dir=project[targetname]+os.sep

    if conf.t.basetype=='php' or conf.t.basetype=='tex' or conf.t.basetype=='html':

        fileext='.'+conf.t.basetype

        # These should not be neccesary to set, as they depend on
        # reStructuredTex, so they are impossible to change
        # Still needed for printing the code
        conf.t.hb='<H2>'
        conf.t.he='</H2>'

        print "Creating list of files"
        # Search the Contents files for all files to process
        allfiles=[]
        lookupsubdir={}
        for fname in conf.t.indexfiles:
            P=ContentsPrinter(conf,fname)

            # Create list of files with subdir appended	
            subdir,fname=os.path.split(fname)
            for name in P.files:
                allfiles.append(os.path.join(subdir,name))
                lookupsubdir[name]=subdir

        conf.lookupsubdir=lookupsubdir

        print "Writing internal refs"
        f=open(conf.g.root+'mat2doc/ref-internal.txt','w')

        for funname in lookupsubdir.keys():
            f.write('.. |'+funname+'| replace:: `'+funname+'`\n')
            f.write('.. _'+funname+': '+conf.t.urlbase+lookupsubdir[funname]+
                    os.sep+funname+'.php\n')

        # flush the file, because we need it again very quickly
        f.flush()
        f.close()

        # Print Contents files
        lookupsubdir={}
        for fname in conf.t.indexfiles:
            P=ContentsPrinter(conf,fname)
            if do_rebuild_file(conf.g.root+fname+'.m',
                               conf.t.dir+fname+fileext,
                               rebuildmode):
                P.write_the_file()


        for fname in allfiles:
            if do_rebuild_file(conf.g.root+fname+'.m',
                               conf.t.dir+fname+fileext,
                               rebuildmode):
                print 'Rebuilding '+conf.t.basetype+' '+fname

                P=matfile_factory(conf,fname)
                P.write_the_file()
                

    if conf.t.basetype=='mat':

        for root, dirs, files in os.walk(conf.t.dir, topdown=False):
            # Walk through the .m files
            for mfile in filter(lambda x: x[-2:]=='.m',files):
                print 'MAT '+os.path.join(root,mfile)
                print_matlab(conf,os.path.join(root,mfile),
                             os.path.join(root,mfile))
        

    if conf.t.basetype=='verify':
        # The code below is so old, that it does not work.
        #

        # Do this with functional programming and lambda functions
        # just because I can.

        # Get all distributed files.
        # The "distriblist" function has
        # been removed, replace the following line with something
        # modern.
        #namelist=distriblist(conf.g.root)

        # Keep only regular files (remove directories)
        namelist=filter(lambda x:x[1]=='f',namelist)

        # Skip the filetype from the data
        namelist=map(lambda x:x[0],namelist)

        # Keep only .m files.
        namelist=filter(lambda x:x[-2:]=='.m',namelist)
        #namelist=filter(lambda x:x[-6:]!='init.m',namelist)
        #namelist=filter(lambda x:x[-10:]!='Contents.m',namelist)

        for name in namelist:
            ignored=0
            for s in conf.t.ignore:
                if re.compile(s).search(name, 1):
                    ignored=1

            if ignored==1:
                print 'IGNORED',name
            else:
                print name
                f=open(conf.g.root+name);
                buf=f.read()
                f.close()

                for target in conf.t.targets:
                    if buf.find(target)==-1:
                        print '    ',target,'is missing'

                for notappear in conf.t.notappears:
                    pos=buf.find(notappear)
                    if pos>1:
                        endpos=buf.find('\n',pos)
                        print '    ',buf[pos:endpos]




def printnoteshtml(noteprefix,notesdir,notehtml):

    # Parse the authors file in the scripts directory
    print notesdir+'scripts/authors'
    authdict = parseauthors(notesdir+'scripts/authors')

    # Get information from all the 'config' files
    allnotesdict=parseconfigfiles(noteprefix,notesdir,authdict)
    notes=allnotesdict.keys()

    # Clear the target directory
    rmrf(notehtml)

    #keys=getcurrentnotes(allnotesdict)
    keys=allnotesdict.keys()
    keys.sort()

    createindexpage(noteprefix,notesdir,allnotesdict,keys,notehtml+'by_number.php')
    
    keys.sort(key=lambda x: allnotesdict[x]['year'])

    createindexpage(noteprefix,notesdir,allnotesdict,keys,notehtml+'by_year.php')

    keys.sort()
    keys.sort(key=lambda x: allnotesdict[x]['type'])

    createindexpage(noteprefix,notesdir,allnotesdict,keys,notehtml+'by_type.php')

    keys.sort()
    keys.sort(key=lambda x: allnotesdict[x]['author'][0]['name'].split(' ')[:-1])

    createindexpage(noteprefix,notesdir,allnotesdict,keys,notehtml+'by_author.php')

    for note in notes:                
        notename=noteprefix+note
        shutil.copy2(notesdir+note+os.sep+notename+'.pdf',
                     notehtml+notename+'.pdf')

        if allnotesdict[note]['bibentry']:
            shutil.copy2(notesdir+note+'/bibentry',
                         notehtml+notename+'.bib')

        if allnotesdict[note]['poster']:
            shutil.copy2(notesdir+note+'/poster.pdf',
                         notehtml+notename+'_poster.pdf')

        if allnotesdict[note]['slides']:
            shutil.copy2(notesdir+note+'/slides.pdf',
                         notehtml+notename+'_slides.pdf')

        if allnotesdict[note]['rr']:
            shutil.copy2(notesdir+note+'/rr.zip',
                         notehtml+notename+'_rr.zip')

if 'printdoc' in sys.argv[0]:
    if len(sys.argv)==1:
        print 'Usage: printdoc.py [project] [target] <rebuildmode>'
        sys.exit()
    
    if len(sys.argv)==2:
        printdoc(sys.argv[1],sys.argv[2])
    else:
        printdoc(sys.argv[1],sys.argv[2],sys.argv[3])
