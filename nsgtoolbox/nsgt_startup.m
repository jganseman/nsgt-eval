function nsgt_startup()
%NSGT_STARTUP  Set paths for using NSGToolbox
%   Usage: nsgt_startup()
%
%   This script file adds NSGToolbox folders to the MATLAB path.

% Author: Nicki Holighaus
% Date: 20.03.13

basepath = which('nsgt_startup');
basepath=basepath(1:end-15);

addpath(basepath);

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
  
  addpath([basepath,filesep,d(ii).name]);    
end

basepath = [basepath,filesep,'wrappers'];
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
  
  addpath([basepath,filesep,d(ii).name]);    
end