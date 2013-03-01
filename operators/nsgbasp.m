function G = nsgbasp(g,shift,M,Ls,phaselock);

% NSGBASP.M - Nicki Holighaus 01.06.2012
%
% G = nsgbasp(g,shift,M,Ls,phaselock);
% 
% Given the cell array 'g' of windows and the time shift vector 'shift',
% this function computes the corresponding non-stationary gabor analysis
% matrix.
%
% !Attention!: 	While this routine can be used to gain some insight into 
%		the structure of frame-related operators, it is not suited
%		for use with transform lentghs over a few thousand samples.
% 
% Input: 
%
%           g           : Cell array of frequency side analysis windows
%           shift       : Vector of time shifts
%           M           : Number of frequency channels (optional)
%	    Ls		: Transform length
%	    phaselock   : This optional 0/1 switch specifies the phaselock 
%			  convention: 0: non-phaselocked (default) / 
%                         1: phaselocked 
%
% Output:
%	    G		: Frame analysis operator corresponding to the
%			  input arguments
%
% Usage 	: nsgbasp(g,shift,M(optional),Ls(optional),phaselock(optional))
%                 nsgbasp(g,shift,M(optional),phaselock(optional))
%                 nsgbasp(g,shift,phaselock(optional))
%
% Note that, to allow the second and third mode of usage, M and Ls (if specified) 
% should always be greater than 1.
%
% More information can be found at:
% http://univie.ac.at/nonstatgab/

% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/3.0/ or send a letter to 
% Creative Commons, 444 Castro Street, Suite 900, Mountain View, California, 94041, USA.

% Check input

if nargin == 4
    if Ls == 1;
	phaselock = 1;
	Ls = sum(shift);
    elseif Ls == 0;
        phaselock = 0;
        Ls = sum(shift);
    else 
        phaselock = 0;
    end
end
if nargin == 3
    Ls = sum(shift);
    if sum(M) == 1
	phaselock = 1;
	for kk = 1:length(shift)
            M(kk) = length(g{kk}); 
        end
    elseif sum(M) == 0
        phaselock = 0;
        for kk = 1:length(shift)
            M(kk) = length(g{kk}); 
        end
    else 
        phaselock = 0;
    end
end
if nargin < 3
    if nargin < 2 
	error('Not enough input arguments');
    end
    phaselock = 0;
    Ls = sum(shift);
    for kk = 1:length(shift)
        M(kk) = length(g{kk}); 
    end
end

if max(size(M)) == 1
    M = M(1)*ones(length(shift),1);
end

G = zeros(sum(M),Ls);

% Setup the necessary parameters 
N = length(shift);

timepos = cumsum(shift);
NN = timepos(end);
timepos = timepos-shift(1);


% Construct the analysis operator matrix explicitly
rows = [1;1+cumsum(M)];
for ii = 1:N
    Lg = length(g{ii});

    win_range = mod(timepos(ii)+(-floor(Lg/2):ceil(Lg/2)-1),NN)+1;
    G(rows(ii),win_range) = g{ii}([Lg-floor(Lg/2)+1:Lg,1:ceil(Lg/2)]); 

    if phaselock == 1
        modulates = exp(-2*pi*i*(1:M(ii)-1)'*(-floor(Lg/2):ceil(Lg/2)-1)/M(ii));
        G(rows(ii)+(1:M(ii)-1),win_range) =  bsxfun(@times,modulates,G(rows(ii),win_range));
   elseif phaselock == 0
	modulates = exp(-2*pi*i*(1:M(ii)-1)'*(0:Ls-1)/M(ii));
        G(rows(ii)+(1:M(ii)-1),:) =  bsxfun(@times,modulates,G(rows(ii),:));
    else
	error('Invalid phaselock parameter');
    end
end