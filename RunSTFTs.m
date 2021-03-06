% Author: Joachim Ganseman, University of Antwerp, all rights reserved
% december 2013

% This demo tests different implementations of the STFT and shows how they
% need to be run. The accuracy of the STFT-ISTFT loop is tested (it should
% be near perfect). 

% Many STFT implementations are flawed, 
% often because they do not use windows in Matlab with the 'periodic' flag
% (see http://www.dsprelated.com/dspbooks/sasp/Overlap_Add_Decomposition.html#21934 )

% On top of that, many do not implement a unitary transform, therefore the
% STFT-ISTFT loop does not result in an identical signal. Given that
% metrics exist that depend on gain, this can seriously influence results.

% We adapted the STFT implementations we found, such that the full transform loop 
% is unitary. This is done by scaling the output with a window correction factor 
% (see http://www.mathworks.com/matlabcentral/newsreader/view_thread/21751 )
% The window correction factor is sqrt(mean(window^2)), the RMS of the window.

% For some implementations, that don't use windows both ways, we just
% added a correction factor after the ISTFT. Those implementations are for
% the moment only guaranteed to work with Hann windows.

% for more on windowed STFTs and the properties that windows must adhere to:
% http://www.katjaas.nl/FFTwindow/FFTwindow.html
% http://www.katjaas.nl/FFTwindow/FFTwindow&filtering.html


%%
clear;
addpath(genpath('./'));     % add subdirectories to path


%% Parameters

fftsize = 1024;
hopsize = fftsize / 4;      % 75% is a must if STFT coefs are going to be 
                            % changed -> needs also window at synthesis

bins_per_octave = 48;       % make this a multiple of 12 for music


% Read file
datadir = getDataDirectory();       % directory with example files
[origMix,samplerate] = wavread([datadir 'pianoclip4notes.wav']);

  % if stereo, make mono
  if size(origMix, 2) > 1
     origMix = sum(origMix, 2) ./2 ; 
  end


% Because the signal does not start with 0 but our windows do, we zero-pad
% the beginning as to not lose any information.
% TODO: can be done with less padding: See Hodgkinson STFT implementation.
input = [ zeros(fftsize,1) ; origMix ; zeros(fftsize,1) ];
% afterwards, chop as follows:   
% InvOrigMix = InvOrigMix( fftsize +1 : fftsize +length(origMix) );

  
disp('--- Starting ---')
disp('--- Smaragdis and labROSA, use Hann at FFT and IFFT ---')
% use an overlap of 75% , since we double window (at analysis and synthesis)

disp('*** STFT, Smaragdis implementation (periodic Hann):')  
  
tic
% Do STFT
SmarSTFT = stft(input, fftsize, hopsize, 0, 'hann');

% invert
InvSmarSTFT = stft(SmarSTFT, fftsize, hopsize, 0, 'hann');
InvSmarSTFT = InvSmarSTFT( fftsize +1 : fftsize +length(origMix) )';
toc
% calculate reconstruction error
rec_err = norm(origMix-InvSmarSTFT)/norm(origMix);
fprintf(['  Normalized Reconstruction Error:'...
    '   %e \n'],rec_err);

% property of STFT distributed with Smaragdis' PLCA code: does not rescale.
% This could be attributed to two errors in the implementation:
% - not using periodic windows, which violated the STFT COLA condition
% - not using a window correction factor, so the transform was not unitary.
% This has been fixed in this repository.



% Now test implementation from LabRosa's (Dan Ellis et al) pvsample demo.

% this implementation computes less frames, ends up with a signal
% smaller than the original. So, zero pad before forward fft here.

% NOTE: this implementation rescales 1-on-1 with these windows at Analysis/Synthesis: 
% - ana: hann, syn: rect, overl: 50%
% - ana: hann, syn: hann, overl: 25% 

disp('*** STFT, LabRosa implementation (periodic hann):')
% parameters: mixture, fftsize, windowvector, hopsize, samplerate
tic
EllisSTFT = stft_ellis(input, fftsize, hann(fftsize, 'periodic')', hopsize, samplerate);
InvEllisSTFT = istft_ellis(EllisSTFT, fftsize, hann(fftsize, 'periodic')', hopsize)';  %ones(fftsize, 1)', hopsize)';      
InvEllisSTFT = InvEllisSTFT(fftsize+1:fftsize+length(origMix));
toc

rec_err = norm(origMix-InvEllisSTFT)/norm(origMix);
fprintf(['  Normalized Reconstruction Error:'...
    '   %e \n'],rec_err);

% Problems with this implementation:
% - Renormalization to account for overlap is not in the ISTFT:
%     different hopsizes will not lead to unity. f = f / (sz/hp);
% - Window correction factor needed to be applied for general windows
% This has been fixed in this repository.


disp('--- CATbox and Hodgkinson, use Hann at FFT and rectangular at IFFT ---')
% The following implementations also work with a 50% overlap (unwindowed synthesis)

% Now try an implementation from Dubnov's CATbox
disp('*** STFT, CATbox implementation (periodic hann):')  
% parameters: signal, window, overlap, fftsize
tic
DubnovSTFT = stft_catbox(input, hann(fftsize, 'periodic'), fftsize-hopsize, fftsize);
% parameters: stft, HOPS PER WINDOW , fftsize, 'perfect'/'smooth'
InvDubnovSTFT = istft_catbox(DubnovSTFT, fftsize / hopsize, fftsize, 'perfect')';
InvDubnovSTFT = InvDubnovSTFT(fftsize+1:fftsize+length(origMix));
toc

rec_err = norm(origMix-InvDubnovSTFT)/norm(origMix);
fprintf(['  Normalized Reconstruction Error:'...
    '   %e \n'],rec_err);

% This implementation does the hop factor normalization correctly, but only
% haphazardly had some form of window correction factor. In the code, this
% has been updated to also work in the 'smooth' case.
% The 'perfect' version uses Hann on analysis and rectangular on synth.
% The 'smooth' version uses Hann on both.
% For this implementation to be a unitary transform, use forward window with
% average value 0.5 ?
% TODO : adapt for general windows ? 



% Now try an implementation from Matthieu Hodgkinson

% ( http://www.cs.nuim.ie/~matthewh/ISTFT.m )
disp('*** STFT, Hodgkinson implementation (periodic hann):') 
tic
[HodgSTFT, indices] = stft_hodg(input, hann(fftsize, 'periodic'), hopsize, fftsize);
InvHodgSTFT = istft_hodg(HodgSTFT, indices, fftsize);
InvHodgSTFT = InvHodgSTFT(fftsize+1:fftsize+length(origMix));
toc

rec_err = norm(origMix-InvHodgSTFT)/norm(origMix);
fprintf(['  Normalized Reconstruction Error:'...
    '   %e \n'],rec_err);

% This implementation needed a hop correction factor added after the 
% ISTFT, to compensate for window application in the forward transform.
% It already does appropriate zero-padding itself, therefore it's the only
% one that implemented a true exact STFT from scratch.
% For this implementation to be a unitary transform, use forward window with
% average value 0.5 ?
% TODO : adapt for general windows ? 


disp('--- Finished ---')
