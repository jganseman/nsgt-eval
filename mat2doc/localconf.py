# Setting in this file are specific for this particular user on this
# machine

# The path definitions below are not necessary, they are just an easy
# way of writing the rest of this file
publishdir ='~/Documents/publish/nsg'

# ---- Global configuration options for all projects  ------------------

# Must be defined if you use Octave, absolute path of the Octave interpreter
octexec='/usr/bin/octave'

# Must be defined if you use Matlab, absolute path of the Matlab interpreter
matexec='/usr/bin/matlab'

# Which program should be used for plotting: 'matlab' or 'octave'
plotengine='matlab'

# Where to put temporary files
tmpdir='/tmp/'

# Directory where the mat2doc files are stored
mat2docdir='~/Documents/git/nsgtoolbox/mat2doc'

# Directory where files for upload are produced
filesdir='~/Documents/publish/nsg'

# ---- Project specific setups -----------------------------------------
#
# This file must define a structure for each project, with the
# following fields. Assume we consider the project 'myproject'
#
#   Source for the project. Must always be defined
#     myproject['dir'] = '/some/dir/myproject'
#   
#   For each target, you must define a field. This is where the output
#   is stored:
#
#     myproject['mat'] = '/some/dir/myprojectmat/'
#     myproject['php'] = '/some/dir/myprojectphp/'
#     myproject['pdf'] = '/some/dir/myprojectpdf/'

nsg={}
nsg['dir']      = '~/Documents/git/nsgtoolbox/'
nsg['mat']      = publishdir+'nsgfiles/'
nsg['php']      = publishdir+'nsgphp/doc'
nsg['phplocal'] = publishdir+'nsglocal/doc'
nsg['tex']      = publishdir+'nsgtex/'
nsg['html']     = publishdir+'nsghtml/'
