function S = nsgfrmmat(g,shift,M,Ls)
%NSGFRMMAT  Sparse nonstationary Gabor frame operator matrix
%   Usage:  S = nsgfrmmat(g,shift,M,Ls)
%           S = nsgfrmmat(g,shift,M)
%           S = nsgfrmmat(g,shift)
%
%   Input parameters:
%         g         : Cell array of window functions
%         shift     : Vector of shifts between the center frequencies
%         M         : Vector of lengths of the window functions
%   Output parameters:
%         S         : Sparse frame operator matrix
%
%   Creates the frame operator matrix of size $LsxLs$ associated to the 
%   nonstationary Gabor system *g*, *shift*, *M* in sparse matrix format.
%

% Author: Nicki Holighaus
% Date: 04.03.13

if nargin < 2
   error('Too few input arguments');
end
    
N = length(g);
timepos = cumsum(shift)-shift(1);

Lg = cellfun(@length,g);
g = cellfun(@fftshift,g,'UniformOutput',0);

if nargin < 4
    Ls = timepos(end)+shift(1);
    if nargin < 3
        M = Lg;
    end    
end
    
B = floor((Lg-1)./M);
S0 = sparse(Ls,Ls);

for ll = 1:N
    
    win_range = mod(timepos(ll)+(-floor(Lg(ll)/2):ceil(Lg(ll)/2)-1),Ls)+1;
    S0(win_range,1) = S0(win_range,1) + M(ll)*abs(g{ll}).^2;
    
    for kk = 1:B(ll)
        temp0 = M(ll)*g{ll}(1+kk*M(ll):end).*conj(g{ll}(1:end-kk*M(ll)));
        win_begin = win_range(1:end-kk*M(ll));
        win_end = win_range(1+kk*M(ll):end);
        
        S0(win_begin,1+kk*M(ll)) = S0(win_begin,1+kk*M(ll)) + conj(temp0);
        S0(win_end,Ls+1-kk*M(ll)) = S0(win_end,Ls+1-kk*M(ll)) + temp0;
    end
end

[i,j] = find(S0); 
ind0 = sub2ind([Ls,Ls],i,j);
S = sparse(i,mod(i+j-2,Ls)+1,S0(ind0),Ls,Ls);