% NSGT_STARTUP.M - Nicki Holighaus 01.03.12
%
% This script file adds the NSGT Toolbox to the MATLAB path.

% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/3.0/ or send a letter to 
% Creative Commons, 444 Castro Street, Suite 900, Mountain View, California, 94041, USA.

basepath = which('nsgt_startup');
basepath=basepath(1:end-15);

addpath(basepath);

bp=[basepath,filesep];

d = dir(basepath);

for ii=1:length(d)
  
  % We only look for directories
  if ~d(ii).isdir
    continue;
  end;
  
  % Skip the default directories . and ..
  if (d(ii).name(1)=='.')
    continue;
  end;
  
  addpath([bp,d(ii).name]);
end