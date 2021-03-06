function [f_sliced,sl_len,tr_area] = slicing(f,sl_len,tr_area,Ls) 
%SLICING  Cut a signal into uniform slices with half-overlap
%   Usage:  f_sliced = slicing(f,sl_len,tr_area,Ls)
%           f_sliced = slicing(f,sl_len,tr_area)
%
%   Input parameters:
%         f         : Signal to be sliced
%         sl_len    : Slice length (in samples, must be even)
%         tr_area   : Length of each transition area (in samples,
%                     optional, default is $ceil( sl\_len/16 )$)
%         Ls        : Length of *f* (optional)
%   Output parameters:
%         f_sliced  : Matrix containing the slices of f as columns
%         sl_len    : Possibly corrected slice length (in samples, even)
%         tr_area   : Possibly corrected transition area length (in 
%                     samples, even)
%    
%   This function cuts a signal into compactly supported pieces of length
%   *sl_len* using a uniform partition of unity composed of Tukey windows
%   with plateau area $sl\_len/2-tr\_area$ and transition areas of length
%   *tr_area*. The resulting signal slices are stored in the columns of the
%   output.
%    
%   See also:  slicq, islicq, unslicing
%
%   References:  dogrhove12

% Author: Nicki Holighaus
% Date: 26.04.13

if nargin < 4
    Ls = length(f);
    if nargin < 2
        error('Not enough input arguments');
    end
end

if nargin < 3 || tr_area > ceil(sl_len/2)
    tr_area = 2*ceil(sl_len/16);
end

if mod(sl_len,2) == 1
    sl_len = 2*ceil(sl_len/2);
    warning('SLICING:oddlength','sl_len was increased by 1 as to be even');
end

if sl_len > Ls
    error('Reduce slice length or increase signal length');
end

if mod(tr_area,2) == 1
    tr_area = 2*ceil(tr_area/2);
    warning('SLICING:oddlength',['tr_area was increased by 1 ',...
        'as to be even']);
end

rows = ceil(2*Ls/sl_len);
hopsize = sl_len/2;

f = [f;zeros(hopsize*rows-Ls,1)];

f_sliced = zeros(sl_len,rows);

% Construct Tukey window for slicing
tw = winfuns('hann',2*tr_area);
tw = [tw(tr_area+1:end);ones(sl_len/2-tr_area,1);tw(1:tr_area)];

% The columns of the following matrix contain the indices of (possibly)
% non-zero entries in [odd slices; even slices; h_0 (0-centered Tukey
% windows used for slicing)]

loc = [sl_len-tr_area/2+1:sl_len,1:hopsize+tr_area/2;...
    hopsize-tr_area/2+1:sl_len,1:tr_area/2;...
    -tr_area/2+1:hopsize+tr_area/2];

f_sliced(loc(1,:),1) = ...
    f([end-tr_area/2+1:end,1:hopsize+tr_area/2]).*tw;

for kk=2:2:rows-1
    f_sliced(loc(2,:),kk) = f((kk-1)*hopsize+loc(3,:)).*tw;
end
for kk=3:2:rows-1
    f_sliced(loc(1,:),kk) = f((kk-1)*hopsize+loc(3,:)).*tw;
end

f_sliced(loc(mod(rows-1,2)+1,:),rows) = ...
    f([end-hopsize-tr_area/2+1:end,1:tr_area/2]).*tw;