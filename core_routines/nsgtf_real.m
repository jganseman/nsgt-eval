function [c,Ls] = nsgtf_real(f,g,shift,M)
%NSGTF_REAL  Nonstationary Gabor filterbank for real signals
%   Usage: [c,Ls] = nsgtf_real(f,g,shift,M)
%          [c,Ls] = nsgtf_real(f,g,shift)
%          c = nsgtf_real(...)
%
%   Input parameters: 
%         f         : A real-valued signal to be analyzed (For multichannel
%                     signals, input should be a matrix which each
%                     column storing a channel of the signal).
%         g         : Cell array of Fourier transforms of the analysis 
%                     windows
%         shift     : Vector of frequency shifts
%         M         : Number of time channels (optional)
%                     If M is constant, the output is converted to a
%                     matrix
%   Output parameters:
%         c         : Transform coefficients (matrix or cell array)
%         Ls        : Original signal length (in samples)
%
%   Given the cell array *g* of windows, the time shift vector *shift*, and
%   channel numbers *M*, `nsgtf_real` computes the corresponding 
%   nonstationary Gabor filterbank of *f*, using only the filters with at 
%   least partially supported on the positive frequencies. Let 
%   $P(n)=\sum_{l=1}^{n} shift(l)$, then the output 
%   `c = nsgtf_real(f,g,shift,M)` is a cell array with 
%
%   ..         Ls-1                                      
%      c{n}(m)= sum fft(f)(l)*conj(g\{n\}(l-P(n)))*exp(2*pi*i*(l-P(n))*m/M(n))
%               l=0                                      
%
%   .. math:: c\{n\}(m) = \sum_{l=0}^{Ls-1} \hat{f}[l]\overline{g\{n\}[l-P(n)]}e^{-2\pi i(l-P(n))m/M(n)},
%
%   where `m` runs from `0` to `M(n)-1` and `n` from 1 to `N`, where
%   $g\{N\}$ is the final filter at least partially supported on the
%   positive frequencies. All filters in *g*, *shift* that are completely
%   supported on the negative frequencies are ignored.
%
%   For more details, see |nsgtf|.
%
%   See also:  nsigtf_real, nsdual, nstight

% Author: Nicki Holighaus, Gino Velasco
% Date: 23.04.13

% Check input arguments
if nargin < 5
    SL = 0;
end

if nargin <= 2
    error('Not enough input arguments');
end

[Ls,CH]=size(f);

if Ls == 1
    f = f.';
    Ls = CH;
    CH = 1;
end

if CH > Ls
    disp(['The number of signal channels (',num2str(CH),') ',...
        'is larger than']);
    disp(['the number of samples per channel (',num2str(Ls),').']);
    reply = input('Is this correct? ([Y]es,[N]o)', 's');
    switch reply
        case {'N','n','No','no',''}
            reply2 = input('Transpose signal matrix? ([Y]es,[N]o)', 's');
            switch reply2
                case {'N','n','No','no',''}
                    error('Invalid signal input, terminating program');
                case {'Y','y','Yes','yes'}
                    disp('Transposing signal matrix and continuing ',...
                        'program execution');
                    f = f.';
                    X = CH; CH = Ls; Ls = CH; clear X;
                otherwise
                    error('Invalid reply, terminating program');
            end
        case {'Y','y','Yes','yes'}
            disp('Continuing program execution');
        otherwise
            error('Invalid reply, terminating program');
    end
end

N=length(shift);    % The number of frequency slices

if nargin == 3
    M = zeros(N,1);
    for kk = 1:N
        M(kk) = length(g{kk});
    end
end

if max(size(M)) ==  1
    M = M(1)*ones(N,1);
end

% some preparation

f = fft(f);

posit = cumsum(shift)-shift(1); % Calculate positions from shift vector

% A small amount of zero-padding might be needed (e.g. for scale frames)

fill = sum(shift)-Ls;
f = [f;zeros(fill,CH)];

Lg = cellfun(@length,g);
N = find(posit-floor(Lg/2) <= (Ls+fill)/2,1,'last');
c=cell(N,1); % Initialisation of the result

% The actual transform
for ii = 1:N
    idx = [ceil(Lg(ii)/2)+1:Lg(ii),1:ceil(Lg(ii)/2)];
    win_range = mod(posit(ii)+(-floor(Lg(ii)/2):ceil(Lg(ii)/2)-1),...
        Ls+fill)+1;
    
    if M(ii) < Lg(ii) % if the number of frequency channels is too small,
        % aliasing is introduced
        col = ceil(Lg(ii)/M(ii));
        temp = zeros(col*M(ii),CH);
        
        temp([end-floor(Lg(ii)/2)+1:end,1:ceil(Lg(ii)/2)],:) = ...
            bsxfun(@times,f(win_range,:),g{ii}(idx));
        temp = reshape(temp,M(ii),col,CH);
        
        c{ii} = squeeze(ifft(sum(temp,2)));
        % Using c = cellfun(@(x) squeeze(ifft(x)),c,'UniformOutput',0);
        % outside the loop instead does not provide speedup; instead it is
        % slower in most cases.
    else
        temp = zeros(M(ii),CH);
        temp([end-floor(Lg(ii)/2)+1:end,1:ceil(Lg(ii)/2)],:) = ...
            bsxfun(@times,f(win_range,:),g{ii}(idx));
        
        c{ii} = ifft(temp);
    end
end

if max(M) == min(M)
    c = cell2mat(c);
    c = reshape(c,M(1),N,CH);
end