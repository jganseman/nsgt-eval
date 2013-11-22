% TestNSGTFlinearity.m
% We assumed there was a problem with the NSGTF implementation, 
% not being entirely linear.
% The code underneath was used to verify this.
% Author: J. Ganseman

%% some initialization
clear;

disp('*** Initializing data...');

            % add subdirectories to path
addpath(genpath('../'));

samplerate = 44100;         % precondition: sampling rate of all files must be 44100.

% STFT parameters
fftsize = 1024;             
hopsize = fftsize / 4;
hopsizepct = hopsize/fftsize;
window = 'hann';

% CQT parameters
bins_per_octave = 48;     %  (multiple of 12 for music)
fmax = samplerate/3;     % center frequency of the highest frequency bin .
fmin = fmax/512;        
% at 44100Hz sampling rate, max=14700, min=28.7
% FYI, a piano goes from 27.5 to 4186 Hz.

% Make a sound mixture file
datadir = getDataDirectory();
source1wav = wavread([datadir 'pianofluid.wav']);         % or choose pianowin/violinwin if you like the WinXP MIDI synth better
source2wav = wavread([datadir 'violinfluid.wav']);
source1wav = makeMono(source1wav);              % depends on makeMono.m function
source2wav = makeMono(source2wav);              % at the same time, adapt gain for mixing
sourceMix = (source1wav+source2wav)./2;

signallength = length(sourceMix);


%% Assert that it holds that the sum of two NSGTs equals the NSGT of the sum

disp('*** Running NSGTF linearity test');

% use the standard windows
[windows,shifts,windowlengths] = nsgcqwin(fmin,fmax,bins_per_octave,samplerate,signallength);

% compute dual windows for the inverse transform
dualwindows = nsdual(windows,shifts,windowlengths);

% compute NSGTF representation of the two sources in the mixture
source1nsgt = nsgtf(source1wav,windows,shifts,windowlengths);
source2nsgt = nsgtf(source2wav,windows,shifts,windowlengths);

% sum the two, elementwise
if iscell(source1nsgt)
    for i = 1:size(source1nsgt,1)
        sumnsgts{i} = source1nsgt{i} + source2nsgt{i};
    end
else                % if we have our results in a matrix representation
    sumnsgts = source1nsgt + source2nsgt;
end

%invert the result, see if we get the source mixture back
assertMix = nsigtf(sumnsgts,dualwindows,shifts,signallength);


rec_err = norm(sourceMix-assertMix)/norm(sourceMix);
fprintf(['Inverse of summed NSGTF transforms Relative error:'...
    '   %e \n'],rec_err);

onlyreal = isreal(assertMix);
fprintf(['Is the NGSTF result fully real? %d \n'], onlyreal);
if ~onlyreal
    largestimag = max(abs(imag(assertMix)));
    fprintf(['Largest imaginary number is %e \n'], largestimag);
    realAssertMix = real(assertMix);
    rec_errReal = norm(sourceMix-realAssertMix)/norm(sourceMix);
    fprintf(['Relative error of real part of signal is:'...
    '   %e \n'],rec_errReal);
end

% % PLOT:
%figure;
%plotnsgtf(sumnsgts',shifts,samplerate,2,60);

% % PLAY:
% soundsc(assertMix, 44100);

% When the result contains imaginary numbers, the file will not be played


%% Assert that additivity of the STFT holds

disp('*** Running STFT linearity test...');

% compute transforms
source1stft = stft(source1wav, fftsize, hopsize, 0, window);
source2stft = stft(source2wav, fftsize, hopsize, 0, window);

% sum transforms
sumstft = source1stft + source2stft;

% invert this
assertMixSTFT = stft(sumstft, fftsize, hopsize, 0, window);

% cut to the right size
assertMixSTFT = assertMixSTFT(1:signallength)';

% doublecheck with original mixture
rec_err = norm(sourceMix-assertMixSTFT)/norm(sourceMix);
fprintf(['Inverse of summed STFT transforms Relative error:'...
    '   %e \n'],rec_err);

% check if result contains imaginary elements
onlyrealSTFT = isreal(assertMixSTFT);
fprintf(['Is the STFT result fully real? %d \n'], onlyrealSTFT);
if ~onlyrealSTFT
    largestimagSTFT = max(abs(imag(assertMixSTFT)));
    fprintf(['Largest imaginary number is %e \n'], largestimagSTFT);
    realAssertMixSTFT = real(assertMixSTFT);
    rec_errRealSTFT = norm(sourceMix-realAssertMixSTFT)/norm(sourceMix);
    fprintf(['Relative error of real part of signal is:'...
    '   %e \n'],rec_errRealSTFT);
end

% % PLOT:
%figure;
%plotnsgtf(sumnsgts',shifts,samplerate,2,60);

% % PLAY:
% soundsc(assertMixSTFT, 44100);
