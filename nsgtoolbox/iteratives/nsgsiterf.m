function [fr,res,Nit]=nsgsiterf(c,g,shift,M,varargin)
%NSGSITERF  Iterative nonstationary Gabor filterbank synthesis
%   Usage:  [fr,res,Nit]=nsgsiterf(c,g,shift,M,Ls,varargin)
%           [fr,res,Nit]=nsgsiterf(c,g,shift,M,varargin)
%           [fr,res,Nit]=nsgsiterf(c,g,shift,M,Ls)
%           [fr,res,Nit]=nsgsiterf(c,g,shift,M)
%           [fr,res]=nsgsiterf(...)
%           fr=nsgsiterf(...)
%
%   Input parameters:
%         c         : Nonstationary Gabor coefficients
%         g         : Cell array of filters
%         shift     : Vector of shifts between the center frequencies
%         M         : Number of time steps
%         Ls        : Original signal length
%         varargin  : Optional input pairs (see table below)
%   Output parameters: 
%         fr        : Synthesized output signal
%         res       : Vector of relative residuals
%         Nit       : Number of iterations
%
%   Given coefficients *c* and nonstationary Gabor frame specified by *g*, 
%   *shift* and *M*, this routine approximates the synthesis formula 
%   associated to the canonical dual frame. 
%
%   The synthesized signal *fr* is obtained by first synthesizing with 
%   respect to *g*, *shift* and *M* using |nsigtf| followed by iteratively
%   applying the inverse frame operator to the result using the conjugate 
%   gradients method. The following equivalence is used:
%
%   ..     N-1 M(n)-1                            N-1 M(n)-1
%       fr=sum  sum  c{n}(m)S^{-1}g_{n,m}=S^{-1}(sum  sum  c{n}(m)g_{n,m}),
%          n=0  m=0                              n=0  m=0
%
%   ..  math:: fr = \sum_{n=0}^{N-1}\sum_{m=0}^{M(n)-1}c\{n\}(m)S^{-1}g_{n,m} = S^{-1}\left(\sum_{n=0}^{N-1}\sum_{m=0}^{M(n)-1}c\{n\}(m)g_{n,m}\right),
%
%   where `N=numel(shift)`. The conjugate gradients algorithm uses the 
%   frame operator, or rather its efficient realization by applying |nsgtf| 
%   and |nsigtf| consecutively.
%
%   Convergence speed of the conjugate gradients algorithm depends on the 
%   condition number of the frame operator, which can be improved by
%   preconditioning. Currently, only a diagonal preconditioner using the 
%   inverse of the frame operator diagonal is implemented.
%
%   Note: The algorithm only converges if *g*, *shift* and *M* form a
%   frame.
%
%   Optional input arguments arguments can be supplied like this::
%
%       nsgsiter(c,g,shift,M,'tol',tol)
%
%   The arguments must be character strings followed by an
%   argument:
%
%     'tol',tol      Error tolerance
%
%     'Mit',Mit      Maximum number of iterations
%
%     'prec',prec    Preconditioning switch
%
%   See also:  nsigtf, nsgtf, nsgaiterf
%
%   References:  nebahoso13 gr93

% Author: Nicki Holighaus
% Date: 04.03.13

if nargin < 2
    error('Not enough input arguments');
end

% Set default parameters
tol = 10^-10;   % Error tolerance
Mit = 200;      % Maximum number of iterations
prec = 0;

if nargin >= 3
    if isnumeric(varargin{1}) && numel(varargin{1}) == 1
        Ls = varargin{1};
        varargin = varargin(2:end);
    end
    Lvar = length(varargin);
    if mod(Lvar,2)
        error('Invalid input argument');
    end
    for kk = 1:2:Lvar
        if ~ischar(varargin{kk})
            error('Invalid input argument');
        end
        switch varargin{kk}
            case {'tol'}
                tol = varargin{kk+1};
            case {'Mit'}
                Mit = varargin{kk+1};
            case {'prec'}
                prec = varargin{kk+1};
            otherwise
                error(['Invalid input argument: ', varargin{kk}]);
        end
    end
end

N = length(shift);
posit = cumsum(shift);
L = posit(end);
posit = posit-shift(1);

frmop = @(x) nsigt(nsgt(x,g,shift,M),g,shift,L);

fr = nsigtf(c,g,shift,L);

fr = fft(fr);

if prec == 0
    [fr,tmp1,tmp2,Nit,res] = pcg(frmop,fr,tol,Mit);
else 
    % Construct the diagonal of the frame operator matrix explicitly
    diagonal=zeros(L,1);
    for ii = 1:N
        Lg = length(g{ii});

        win_range = mod(posit(ii)+(-floor(Lg/2):ceil(Lg/2)-1),L)+1;
        diagonal(win_range) = diagonal(win_range) + ...
            (fftshift(g{ii}).^2)*M(ii);   
    end
    D = spdiags(diagonal,0,L,L);
    [fr,tmp1,tmp2,Nit,res] = pcg(frmop,fr,tol,Mit,D);    
end

fr = ifft(fr);

if  exist('Ls','var')
    fr = fr(1:Ls,:);
end

if nargout>1
   res=res/norm(fr(:));
end