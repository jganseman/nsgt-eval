function [c,g,shift,M,Ls] = OnsetNsgt(f,thre,short,max_win,win_length)

if nargin < 5
    win_length = 2048;
    if nargin < 4 
        max_win = 10;
        if nargin < 3
            short = 192;
            if nargin < 2
                error('%s: Not enough input arguments',upper(mfilename));
            end
        end
    end
end
        
if min(size(f)) > 1
    error('%s: Multichannel signals are not supported',upper(mfilename));
end

Ls = length(f);

positions = onsetdet(f,win_length,thre);

blocks = diff(positions);
idx = find(blocks < 4/3*short);
while numel(idx) > 0
    positions(1+idx) = [];
    blocks = diff(positions);
    idx = find(blocks < 4/3*short);
end

[g,shift,M] = nsgsclwin(positions,short,max_win,Ls);

c = nsgt_real(f,g,shift,M);