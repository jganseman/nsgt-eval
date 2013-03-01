function [g,shift,M] = nsgerbwin(V,sr,L)

% NSGERBWIN.M
%---------------------------------------------------------------
% [g,shift,M]=nsgerbwin(V,sr,L) creates a set of windows for the 
% ERBlet nonstationary Gabor transform. 
%---------------------------------------------------------------
%
% INPUT : V ......... Voices per ERB
%	      sr ........ Sampling rate (in Hz)
%         L ......... Length of signal (in samples)
%
% OUTPUT : g ......... Cell array of window functions.
%          shift ..... Vector of shifts between the center frequencies.
%          M ......... Vector of lengths of the window functions.
%
%
% AUTHOR(s) : Thibaud Necciari, Nicki Holighaus, 2012
%
% EXTERNALS : blackharr

% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/3.0/ or send a letter to 
% Creative Commons, 444 Castro Street, Suite 900, Mountain View, California, 94041, USA.

df = sr/L; % frequency resolution in the FFT

fmin = 0;
fmax = sr/2;

% Convert fmin and fmax into ERB
erblims = 9.2645*sign([fmin,fmax]).*log(1+abs([fmin,fmax])*0.00437); 

% Determine number of freq. channels
Nf = V*ceil(erblims(2)-erblims(1));

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
posit(Nf+1:end) = L-posit(Nf+1:end);% Extension to negative freq.
shift = [L-posit(end); diff(posit)];% Hop sizes in samples 

% Compute desired essential (Gaussian) support for each filter
Lwin = 4*round(gamma/df);

% Blackman-Harris windows are slightly broader than Gaussians, this is 
% offset by the factor 1.1

M = round(Lwin/1.1);

% Compute cell array of analysis filters
g = arrayfun(@(x) blackharr(x)/sqrt(x),M,'UniformOutput',0);

g{1}=1/sqrt(2)*g{1};
g{end}=1/sqrt(2)*g{end};