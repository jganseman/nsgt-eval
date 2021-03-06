# -------------------------------------------
# Global configuration of the mat2doc system
# -------------------------------------------

# When writing this file, certain variables are already defined:
#
#   self.root points to the project directory


import localconf

conf=ConfType()

def mycopyrightfun(self):
    vf=file(self.root+'nsg_version');
    v=vf.readline()
    vf.close
    
    f=file(self.root+'mat2doc/copyrightplate')
    buf=f.readlines()
    f.close

    copyright=['Copyright (C) 2013 Nicki Holighaus.\n','This file is part of NSG toolbox version '+v]
    copyright.extend(buf)
    
    return copyright

conf.copyright=mycopyrightfun

conf.urlext='php'

contentsfiles=['/Contents','demos/Contents','core_routines/Contents','wrappers/Contents','windows/Contents','helpers/Contents','plotting/Contents','iteratives/Contents','generators/Contents','operators/Contents']

#'demos/Contents'

# ------------------------------------------
# Configuration of PHP for Sourceforge
# ------------------------------------------

php=PhpConf()

php.indexfiles=contentsfiles
php.includedir='include/'
php.urlbase='/doc/'
php.codedir=localconf.nsg['mat']

# ------------------------------------------
# Configuration of PHP for local server
# ------------------------------------------
phplocal=PhpConf()

phplocal.indexfiles=contentsfiles
phplocal.includedir='include/'
phplocal.urlbase='/nsglocal/'
phplocal.codedir=localconf.nsg['mat']

# ------------------------------------------
# Configuration of HTML
# ------------------------------------------

html=HtmlConf()
html.urlbase='http://nsg.sourceforge.net/'
html.indexfiles=contentsfiles
html.codedir=localconf.nsg['mat']

# ------------------------------------------
# Configuration of LaTeX
# ------------------------------------------

tex=TexConf()

tex.indexfiles=contentsfiles
tex.urlbase='http://nsg.sourceforge.net/'
tex.codedir=localconf.nsg['mat']
tex.imagetype='png'

# ------------------------------------------
# Configuration of Matlab
# ------------------------------------------

mat=MatConf()
mat.urlbase='http://nsg.sourceforge.net/'
