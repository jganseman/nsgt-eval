%DEMO_NSGTF Sliced constant-Q usage and comparison demo
%
%   This script sets up different nonstationary Gabor filterbank frames 
%   with the specified parameters, computes windows and corresponding 
%   canonical dual windows as well as the transform and reconstruction of a 
%   test signal, and plots the windows and the energy of the coefficients.
%
%   .. figure::
%
%      windows + dual windows (constant-Q)
%
%      This figure shows the window functions used in the constant-Q 
%      filterbank and the corresponding canonical dual windows. 
%
%   .. figure::
%
%      windows + dual windows (sliCQ)
%
%      This figure shows the window functions used in the sliced constant-Q
%      filterbank and the corresponding canonical dual windows. 
%
%   .. figure::
%
%      constantQ/sliCQ spectrogram (absolute value of coefficients in dB)
%
%      This figure shows (color coded) images of the constant-Q and sliced
%      constant-Q coefficient modulus. 
%
%   See also:  slicq, islicq, nsgtf, nsigtf, nsdual, nsgcqwin, 

% Author: Gino Velasco, Nicki Holighaus
% Date: 04.03.13


%% Setup parameters for the constant-Q.

fmin = 130; % Minimum desired frequency (in Hz)

fmax = 22050; % Maximum desired frequency (in Hz)
% fmax is taken to be the Nyquist frequency if not indicated

%fmax = floor(fmin*2^(floor(log2(22050/fmin))));

bins = 12; % Number of bins per octave

% Use this to test the variabe-Q transform
% bins = [12; 24; 36; 48; 12]; % Number of bins per octave (in Hz)

%% Setup parameters for the sliCQ

sl_len = 10000;
tr_area = 1024;

%% Load the test signal.

[s,sr] = wavread('glockenspiel.wav'); name = 'Glockenspiel';

%[s,sr] = wavread('your_own_signal.wav'); name = 'Your own signal';

Ls = length(s); % Length of signal (in samples)

%% Window design
%  Define a set of windows for the nonstationary Gabor transform with
%  resolution evolving over frequency. In particular, the centers of the
%  windows correspond to geometrically spaced center frequencies.

% Conpute constant-Q filters

[gCQ,shiftCQ,MCQ] = nsgcqwin(fminCQ,fmaxCQ,binsCQ,sr,Ls,'winfun','modblackharr');

% Compute corresponding dual windows.

gdCQ = nsdual(gCQ,shiftCQ,MCQ);

%% Calculate the coefficients

cCQ = nsgtf(s,gCQ,shiftCQ,MCQ);

[cSCQ,gSCQ,shiftSCQ,MSCQ] = ...
    slicq(s,fmin,fmax,bins,sl_len,tr_area,sr);

[s_rSCQ,gdSCQ] = islicq(cSCQ,gSCQ,shiftSCQ,MSCQ,Ls,sl_len,tr_area);

%% Plot the windows and spectrograms

% constant-Q windows
figure;

subplot(211); plot_wins(gCQ,shiftCQ);

subplot(212); plot_wins(gdCQ,shiftCQ);

% sliCQ windows

figure;

subplot(211); plot_wins(gSCQ,shiftSCQ);

subplot(212); plot_wins(gdSCQ,shiftSCQ);

% spectrograms

figure;

subplot(211); plotnsgtf(cCQ,shiftCQ,sr,fmin,fmax,bins,2,60);

subplot(212); plotslicq(cSCQ,shiftSCQ,sr,fmin,fmax,bins,2,60);

%% Test reconstruction
s_r = nsigtf(cCQ,gdCQ,shiftCQ,Ls);

% Print relative error of constant-Q reconstruction.
rec_err = norm(s-s_r)/norm(s);

fprintf(['Relative error of constant-Q reconstruction (should be close '...
    'to zero.): %e \n'],rec_err);

s_r = nsigtf(cERB,gdERB,shiftERB,Ls);

% Print relative error of sliCQ reconstruction.
rec_err = norm(s-s_rSCQ)/norm(s);

fprintf(['Relative error of sliCQ reconstruction (should be close '...
    'to zero.): %e \n'],rec_err);