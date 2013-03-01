function S = nsgfrmmat(g,shift,M,L)

% NSGFRMMAT.M - Nicki Holighaus 13.02.13
%---------------------------------------------------------------
% S = nsgfrmmat(g,shift,M,L) creates the frame operator matrix of size LxL
% associated to the nonstationary Gabor system [g,shift,M] in sparse matrix
% format.
%---------------------------------------------------------------
%
% INPUT :  g ......... Cell array of window functions.
%          shift ..... Vector of shifts between the center frequencies.
%          M ......... Vector of lengths of the window functions.
%
% OUTPUT : S ......... Sparse frame operator matrix.
%

% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/3.0/ or send a letter to 
% Creative Commons, 444 Castro Street, Suite 900, Mountain View, California, 94041, USA.


if nargin < 2
   error('Too few input arguments');
end
    
N = length(g);
timepos = cumsum(shift)-shift(1);

Lg = cellfun(@length,g);
g = cellfun(@fftshift,g,'UniformOutput',0);

if nargin < 4
    L = timepos(end)+shift(1);
    if nargin < 3
        M = Lg;
    end    
end
    
B = floor((Lg-1)./M);
S0 = sparse(L,L);

for ll = 1:N
    
    win_range = mod(timepos(ll)+(-floor(Lg(ll)/2):ceil(Lg(ll)/2)-1),L)+1;
    S0(win_range,1) = S0(win_range,1) + M(ll)*abs(g{ll}).^2;
    
    for kk = 1:B(ll)
        temp0 = M(ll)*g{ll}(1+kk*M(ll):end).*conj(g{ll}(1:end-kk*M(ll)));
        win_begin = win_range(1:end-kk*M(ll));
        win_end = win_range(1+kk*M(ll):end);
        
        S0(win_begin,1+kk*M(ll)) = S0(win_begin,1+kk*M(ll)) + conj(temp0);
        S0(win_end,L+1-kk*M(ll)) = S0(win_end,L+1-kk*M(ll)) + temp0;
    end
end

[i,j] = find(S0); 
ind0 = sub2ind([L,L],i,j);
S = sparse(i,mod(i+j-2,L)+1,S0(ind0),L,L);