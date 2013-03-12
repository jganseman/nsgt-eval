%DEMO_NSGT Onset Detection type Nonstationary Gabor transform usage demo 
%
%   This script sets up a nonstationary Gabor frame with the specified
%   parameters, computes windows and corresponding canonical dual windows
%   and a test signal, and plots the windows and the energy of the 
%   coefficients.
%
%   .. figure::
%    
%      onset detection results
%
%      This figure shows a regular spectrogram with marked onsets and the
%      spectral flux function with marked onsets.
%
%   .. figure::
%
%      windows + dual windows
%
%      This figure shows the window functions used and the corresponding
%      canonical dual windows. 
%
%   .. figure::
%
%      spectrogram (absolute value of coefficients in dB)
%
%      This figure shows a (color coded) image of the nsgtf coefficient
%      modulus. 
%
%   EXTERNALS:  NSGT, NSIGT, NSDUAL, ONSETDET
%

% Author: Gino Velasco, Nicki Holighaus
% Date: 04.03.13


%% Setup onset detection parameters and load the signal.

win_length = 4096; % Window length for the onset analysis

thr = 0.7; % Onset detection threshold

area = 8; % Determines the size of the area over which local maxima are taken

multi = 3; % Area multiplier for the peak picking algorithm 

%% Test signals

[s,fs] = wavread('glockenspiel.wav'); name = 'Glockenspiel';

%[s,fs] = wavread('your_own_signal.wav'); name = 'Your own signal';

Ls = length(s); % Length of signal (in samples)

%% Window design
%  Define a set of windows for the nonstationary Gabor transform with
%  resolution evolving over time. In particular, short windows of 192 
%  samples length should be chosen at onset positions
pos = onsetdet(s,win_length,thr,area,multi,1.5,1);

[g,shift] = nsgsclwin(pos,192,8,Ls);

% Compute corresponding dual windows.

gd = nsdual(g,shift);

% Plot the windows and the corresponding dual windows

figure;

subplot(211); plot_wins(g,shift);

subplot(212); plot_wins(gd,shift);

%% Calculate the coefficients

c = nsgt(s,g,shift);

%% Plot the resulting spectrogram

figure;

plotnsgt(c,shift,fs,2,60);

%% Test reconstruction
s_r = nsigt(c,gd,shift,Ls);

% Print relative error of reconstruction.
rec_err = norm(s-s_r)/norm(s);

fprintf(['Relative error of reconstruction (should be close to zero.):'...
    '   %e \n'],rec_err);
