%DEMO_SLICQ Sliced constant-Q usage and comparison demo
%
%   This script sets up constant-Q and sliced CQ nonstationary Gabor 
%   filterbank frames with the specified parameters, computes filters and 
%   corresponding canonical dual filters as well as the transform and 
%   reconstruction of a test signal, and compares the analyis and synthesis
%   filters and spectrograms.
%
%   .. figure::
%
%      filters + dual filters (constant-Q)
%
%      This figure shows the filter functions used in the constant-Q 
%      filterbank and the corresponding canonical dual filters. 
%
%   .. figure::
%
%      filters + dual filters (sliCQ)
%
%      This figure shows the filter functions used in the sliced constant-Q
%      filterbank and the corresponding canonical dual filters. 
%
%   .. figure::
%
%      constant-Q/sliCQ spectrogram (absolute value of coefficients in dB)
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
%  Define a set of filters for the nonstationary Gabor transform with
%  resolution evolving over frequency. In particular, the centers of the
%  filters correspond to geometrically spaced center frequencies.

% Conpute constant-Q filters

[gCQ,shiftCQ,MCQ] = nsgcqwin(fmin,fmax,bins,sr,Ls,'winfun','modblackharr');

% Compute corresponding dual filters.

gdCQ = nsdual(gCQ,shiftCQ,MCQ);

%% Calculate the coefficients

cCQ = nsgtf(s,gCQ,shiftCQ,MCQ);

[cSCQ,gSCQ,shiftSCQ,MSCQ] = ...
    slicq(s,fmin,fmax,bins,sl_len,tr_area,sr);

[s_rSCQ,gdSCQ] = islicq(cSCQ,gSCQ,shiftSCQ,MSCQ,Ls,sl_len,tr_area);

%% Plot the filters and spectrograms

% constant-Q filters
figure;

subplot(211); plot_wins(gCQ,shiftCQ);

subplot(212); plot_wins(gdCQ,shiftCQ);

% sliCQ filters

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

% Print relative error of sliCQ reconstruction.
rec_err = norm(s-s_rSCQ)/norm(s);

fprintf(['Relative error of sliCQ reconstruction (should be close '...
    'to zero.): %e \n'],rec_err);