%DEMO_NSGTF Nonstationary Gabor filterbank usage demo 
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
%      constant-Q spectrogram (absolute value of coefficients in dB)
%
%      This figure shows a (color coded) image of the constant-Q 
%      coefficient modulus. 
%
%   .. figure::
%   
%      windows + dual windows (ERBlet)
%
%      This figure shows the window functions used in the ERBlet filterbank
%      and the corresponding canonical dual windows. 
%   
%   .. figure::
%   
%      ERBlet spectrogram (absolute value of coefficients in dB)
%   
%      This figure shows a (color coded) image of the ERBlet coefficient
%      modulus.
%   
%   .. figure::
%   
%      windows + dual windows (Wavelet)
%   
%      This figure shows the window functions used in the Wavelet 
%      filterbank and the corresponding canonical dual windows. 
%   
%   .. figure::
%   
%      Wavelet spectrogram (absolute value of coefficients in dB)
%   
%      This figure shows a (color coded) image of the Wavelet coefficient
%      modulus. 
%
%   See also:  nsgtf, nsigtf, nsdual, nsgcqwin, nsgerbwin, nsgwvltwin

% Author: Gino Velasco, Nicki Holighaus
% Date: 04.03.13


%% Setup parameters for the constant-Q.

fminCQ = 130; % Minimum desired frequency (in Hz)

fmaxCQ = 22050; % Maximum desired frequency (in Hz)
% fmaxCQ is taken to be the Nyquist frequency if not indicated

%fmaxCQ = floor(fmin*2^(floor(log2(22050/fmin))));

binsCQ = 12; % Number of bins per octave

% Use this to test the variabe-Q transform
% binsCQ = [12; 24; 36; 48; 12]; % Number of bins per octave (in Hz)

%% Setup parameters for the ERBlets.

binsERB = 2;
Qvar = .5;

%% Setup parameters for the Wavelets.

fminWV = 130;
binsWV = 6;
fac = 2^(2/binsWV)-2^(-2/binsWV);
bwWV = fminWV*fac;

%% Load the test signal.

[s,sr] = wavread('glockenspiel.wav'); name = 'Glockenspiel';

%[s,sr] = wavread('your_own_signal.wav'); name = 'Your own signal';

s = s(1:88200);

Ls = length(s); % Length of signal (in samples)

%% Window design
%  Define a set of windows for the nonstationary Gabor transform with
%  resolution evolving over frequency. In particular, the centers of the
%  windows correspond to geometrically spaced center frequencies.

[gCQ,shiftCQ,MCQ] = nsgcqwin(fminCQ,fmaxCQ,binsCQ,sr,Ls);

[gERB,shiftERB,MERB] = nsgerbwin(binsERB,sr,Ls,'Qvar',Qvar);

[gWV,shiftWV,MWV] = nsgwvltwin(fminWV,bwWV,binsWV,sr,Ls);

% Compute corresponding dual windows.

gdCQ = nsdual(gCQ,shiftCQ,MCQ);

gdERB = nsdual(gERB,shiftERB,MERB);

gdWV = nsdual(gWV,shiftWV,MWV);

%% Calculate the coefficients

cCQ = nsgtf(s,gCQ,shiftCQ,MCQ);

cERB = nsgtf(s,gERB,shiftERB,MERB);

cWV = nsgtf(s,gWV,shiftWV,MWV);

%% Plot the windows and spectrograms

% constant-Q

figure;

subplot(211); plot_wins(gCQ,shiftCQ);

subplot(212); plot_wins(gdCQ,shiftCQ);

figure;

plotnsgtf(cCQ,shiftCQ,sr,2,60);

% ERBlet

figure;

subplot(211); plot_wins(gERB,shiftERB);

subplot(212); plot_wins(gdERB,shiftERB);

figure;

plotnsgtf(cERB,shiftERB,sr,2,60);
 
% Wavelet

figure;

subplot(211); plot_wins(gWV,shiftWV);

subplot(212); plot_wins(gdWV,shiftWV);

figure;

plotnsgtf(cWV,shiftWV,sr,2,60);

%% Test reconstruction
s_r = nsigtf(cCQ,gdCQ,shiftCQ,Ls);

% Print relative error of constant-Q reconstruction.
rec_err = norm(s-s_r)/norm(s);

fprintf(['Relative error of constant-Q reconstruction (should be close '...
    'to zero.): %e \n'],rec_err);

s_r = nsigtf(cERB,gdERB,shiftERB,Ls);

% Print relative error of ERBlet reconstruction.
rec_err = norm(s-s_r)/norm(s);

fprintf(['Relative error of ERBlet reconstruction (should be close '...
    'to zero.): %e \n'],rec_err);

s_r = nsigtf(cWV,gdWV,shiftWV,Ls);

% Print relative error of Wavelet reconstruction.
rec_err = norm(s-s_r)/norm(s);

fprintf(['Relative error of Wavelet reconstruction (should be close '...
    'to zero.): %e \n'],rec_err);