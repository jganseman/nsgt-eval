function [g,shift,M] = nsgerbwin(bins,sr,Ls)
%NSGERBWIN  ERBlet dictionary generator
%   Usage:  [g,shift,M]=nsgerbwin(bins,sr,Ls)
%
%   Input parameters: 
%         bins      : Desired bins per ERB
%         sr        : Sampling rate of f (in Hz)
%         Ls        : signal length
%   Output parameters: 
%         g         : Cell array of ERBlets
%         shift     : Vector of frequency shifts
%         M         : Number of time channels
%
%   Creates a set of windows for the ERBlet nonstationary Gabor transform. 
%
%   EXTERNALS: blackharr
%

% Author: Thibaud Necciari, Nicki Holighaus
% Date: 04.03.13

df = sr/Ls; % frequency resolution in the FFT

fmin = 0;
fmax = sr/2;

% Convert fmin and fmax into ERB
erblims = 9.2645*sign([fmin,fmax]).*log(1+abs([fmin,fmax])*0.00437); 

% Determine number of freq. channels
Nf = bins*ceil(erblims(2)-erblims(1));

% Determine center frequencies
erbs = linspace(erblims(1),erblims(2),Nf)';
fc = (1/0.00437)*sign(erbs).*(exp(abs(erbs)/9.2645)-1);
fc(1)=fmin;
fc(end)=fmax;

% Concatenate "virtual" frequency positions of negative-frequency windows
fc = [fc ; flipud(fc(1:end-1))];    

gamma = 24.7*(4.37*fc*1E-3 +1); % ERB scale

% Convert center frequencies in Hz into samples 

posit = round(fc/df);% Positions of center frequencies in samples
posit(Nf+1:end) = Ls-posit(Nf+1:end);% Extension to negative freq.
shift = [Ls-posit(end); diff(posit)];% Hop sizes in samples 

% Compute desired essential (Gaussian) support for each filter
Lwin = 4*round(gamma/df);

% Blackman-Harris windows are slightly broader than Gaussians, this is 
% offset by the factor 1.1

M = round(Lwin/1.1);

% Compute cell array of analysis filters
g = arrayfun(@(x) blackharr(x)/sqrt(x),M,'UniformOutput',0);

g{1}=1/sqrt(2)*g{1};
g{end}=1/sqrt(2)*g{end};