
%EXAMP_NSGTF CQ/VQ-type Nonstationary Gabor filterbank usage demo 
%
%   This script sets up a nonstationary Gabor filterbank frame with the 
%   specified parameters, computes windows and corresponding canonical dual 
%   windows and a test signal, and plots the windows and the energy of the 
%   coefficients.
%
%   FIGURE 1 windows + dual windows
%
%    This figure shows the window functions used and the corresponding
%    canonical dual windows. 
%
%   FIGURE 2 spectrogram (absolute value of coefficients in dB)
%
%    This figure shows a (color coded) image of the nsgtf coefficient
%    modulus. 
%
%   EXTERNALS:  NSGTF, NSIGTF, NSDUAL

% Author: Gino Velasco, Nicki Holighaus
% Date: 04.03.13


%% Setup parameters and load the signal.

fmin = 130; % Minimum desired frequency (in Hz)

fmax = 22050; % Maximum desired frequency (in Hz)
              % fmax is taken to be the Nyquist frequency if not indicated
              
%fmax = floor(fmin*2^(floor(log2(22050/fmin))));               

bins = 12; % Number of bins per octave

%bins = [12; 24; 36; 48; 12]; % Number of bins per octave (in Hz)

%% Test signals

[s,fs] = wavread('glockenspiel.wav'); name = 'Glockenspiel';

%[s,fs] = wavread('your_own_signal.wav'); name = 'Your own signal';

Ls = length(s); % Length of signal (in samples)

%% Window design
%  Define a set of windows for the nonstationary Gabor transform with
%  resolution evolving over frequency. In particular, the centers of the
%  windows correspond to geometrically spaced center frequencies.

[g,shift,M] = nsgcqwin(fmin,fmax,bins,fs,Ls);

% Compute corresponding dual windows.

gd = nsdual(g,shift,M);

% Plot the windows and the corresponding dual windows

figure;

subplot(211); plot_wins(g,shift);

subplot(212); plot_wins(gd,shift);

%% Calculate the coefficients

c = nsgtf(s,g,shift,M);

%% Plot the resulting spectrogram

figure;

plotnsgtf(c,shift,fs,2,60);

%% Test reconstruction
s_r = nsigtf(c,gd,shift,Ls);

% Print relative error of reconstruction.
rec_err = norm(s-s_r)/norm(s);

fprintf(['Relative error of reconstruction (should be close to zero.):'...
    '   %e \n'],rec_err);
