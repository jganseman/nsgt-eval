function str = getDataDirectory()
%GETDATADIRECTORY Returns the place of the data directory
%   The data for the Score-informed Source Separation is not present in the
%   Mercurial repository. Use this function to define a path to the data
%   and the examples.

str = '';

% When running on a mac, load the c4dm-scratch volume first.
% do this by Finder -> Go -> Connect to server landin.eecs.qmul.ac.uk
if ismac
    str = '/Volumes/c4dm-scratch/jga/data/';        % USE THIS if at C4DM
    str = '/Users/jg/Dropbox/data/';            % USE THIS if at home

% When running on Unix, currently I assume you're on a C4DM server    
elseif isunix
    str = '/c4dm-scratch/jga/data/';
    
% On a PC, heh? fill in here wherever you put the example data.    
elseif ispc
    str = 'C:\Users\JG\Dropbox\data\';
    
end

end

