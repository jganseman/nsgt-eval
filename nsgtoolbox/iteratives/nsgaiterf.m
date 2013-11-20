function [c,Ls,res,Nit]=nsgaiterf(f,g,shift,M,varargin)
%NSGAITERF  Iterative nonstationary Gabor filterbank analysis
%   Usage:  [c,Ls,res,Nit]=nsgaiterf(f,g,shift,M,varargin)
%           [c,Ls,res,Nit]=nsgaiterf(f,g,shift,M)
%           [c,Ls,res]=nsgaiterf(...)
%           [c,Ls]=nsgaiterf(...)
%           c=nsgaiterf(...)
%
%   Input parameters:
%         f         : Input signal
%         g         : Cell array of filters
%         shift     : Vector of shifts between the center frequencies
%         M         : Number of time steps
%         varargin  : Optional input pairs (see table below)
%   Output parameters: 
%         c         : Filterbank coefficients
%         Ls        : Input signal length
%         res       : Vector of relative residuals
%         Nit       : Number of iterations
%
%   Given a function *f* and nonstationary Gabor filterbank frame specified 
%   by *g*, *shift* and *M*, this routine approximates the frame 
%   coefficients associated to the canonical dual frame. 
%
%   The approximated coefficients are obtained by first applying the 
%   inverse frame operator to *f* iteratively using the conjugate gradients 
%   method, followed by computing the analysis coefficients of $S^{-1}f$ 
%   with respect to *g*, *shift* and *M* with |nsgtf|, followed. The 
%   following equivalence is used:
%
%   ..  c{n}(m) = < f , S^{-1} g_{n,m} > = < S^{-1} f , g_{n,m} >
%
%   ..  math:: c\{n\}(m) = \langle f, \mathbf{S}^{-1} g_{n,m} \rangle = \langle \mathbf{S}^{-1} f, g_{n,m} \rangle
%
%
%   The conjugate gradients algorithm uses the frame operator, or rather 
%   its efficient realization by applying |nsgtf| and |nsigtf| 
%   consecutively.
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
%       nsgaiterf(f,g,shift,M,'tol',tol)
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
%   See also:  nsgtf, nsigtf, nsgsiterf
%
%   References:  nebahoso13 gr93

% Author: Nicki Holighaus
% Date: 24.04.13

if nargin < 2
    error('Not enough input arguments');
end

% Set default parameters
tol = 10^-10;   % Error tolerance
Mit = 200;      % Maximum number of iterations
prec = 0;

if nargin >= 3
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

f = fft(f);
Ls = length(f);
N = length(shift);
posit = cumsum(shift)-shift(1);


frmop = @(x) nsigt(nsgt(x,g,shift,M),g,shift,Ls);

if prec == 0
    [f,tmp1,tmp2,Nit,res] = pcg(frmop,f,tol,Mit);
else 
    % Construct the diagonal of the frame operator matrix explicitly
    diagonal=zeros(Ls,1);
    for ii = 1:N
        Lg = length(g{ii});

        win_range = mod(posit(ii)+(-floor(Lg/2):ceil(Lg/2)-1),Ls)+1;
        diagonal(win_range) = diagonal(win_range) + ...
            (fftshift(g{ii}).^2)*M(ii);   
    end
    D = spdiags(diagonal,0,Ls,Ls);
    [f,tmp1,tmp2,Nit,res] = pcg(frmop,f,tol,Mit,D);    
end

f = ifft(f);
c = nsgtf(f,g,shift,M);

if nargout>1
   res=res/norm(f(:));
end