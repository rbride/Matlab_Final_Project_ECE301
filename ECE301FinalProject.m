%Ryan Bride and Steven Larkin

clear all; 
close all;
clc

% Pass bands
wp1 = (0.4*pi);
wp2 = (0.6*pi); 
%Stop Bands
ws1 = (0.2*pi);
ws2 = (0.8*pi);

%converting from rad/s to samples/s
fp1 = (wp1/(2*pi)); % wp1
fp2 = (wp2/(2*pi)); % wp2
fs1 = (ws1/(2*pi)); % ws1
fs2 = (ws2/(2*pi)); % ws2
 
passTol = 0.01;
stopTol = 0.057501127785;  

passTolLSEM = 0.02; % for Least Square Error Minimization 
stopTolLSEM = 0.01;


% -- FIR filters --
% Rectangular window (Pass)
nPass = 69; % number of samples
orderPass = (nPass - 1);

rectWindo = figure(1);
wnPass = [fp1 fp2];
bPass = fir1(orderPass,wnPass, 'bandpass', rectwin(nPass)');
[hPass, wPass] = freqz(bPass);
freqz(bPass, 1, 512)
saveas(gcf, 'RectangularWindow.png');


% Hamming window (Hamm)
nHamm = 33; %number of samples             
orderHamm = (nHamm - 1);

hammingWin = figure(2);
wnHamm = [fp1 fp2];
bHamm = fir1(orderHamm, wnHamm, 'bandpass', hamming(nHamm)');
freqz(bHamm, 1, 512)
saveas(gcf, 'HammingWindow.png');


% Frequency-sampling (fs)
orderfs = 33;

freqSamp = figure(3);
ffs = [0 0.2 0.2 0.7 0.7 1];
mfs = [0 0 1 1 0 0];
bfs = fir2(orderfs, ffs, mfs, hamming(orderfs + 1));
[hfs, wfs, ffs] = freqz(bfs, 1, 512);
freqz(bfs, 1, 512)
saveas(freqSamp, 'FrequencySample.png');


% Parks-McCellan Optimal Equiripple Filter (PMOE)
nPMOE = 32; 

fPMOE = [0 fs1 fp1 fp2 fs2 1];
aPMOE = [0 0.0 1.0 1.0 0.0 0.0];

PMOE = figure(4);
wPMOE = [(stopTol/passTol) (1) (stopTol/passTol)];
bPMOE = firpm(nPMOE, fPMOE, aPMOE, wPMOE);
freqz(bPMOE, 1, 512)
saveas(PMOE, 'Parks-McCellanOptimalEquiripple.png');


% Least Square Error Minimization (LSEM)
nLSEM = 32;

fLSEM = [0 fs1 fp1 fp2 fs2 1];
aLSEM = [0 0.0 1.0 1.0 0.0 0.0];
LSEM = figure(5);

wLSEM = [(stopTolLSEM/passTolLSEM) (1) (stopTolLSEM/passTolLSEM)];
bLSEM = firls(nLSEM, fLSEM, aLSEM, wLSEM);
freqz(bLSEM, 1, 512)   
saveas(LSEM, 'LeastSquareErrorMinimization.png');


% --- IIR filters ---

% Butterworth (Butt)
nButt = 12; % 12 chosen visually as a best fit, Order = N*2.

orderButt = nButt-1;
wnButt = [fp1 fp2];

Butt = figure(6);
[bButt, aButt] = butter(orderButt, wnButt, 'bandpass');
[hButt, wButt] = freqz(bButt, aButt);
freqz(bButt, 1, 512)
saveas(Butt, 'Butterworth.png');



% Chebyshev Type I (Cheb1)
nCheb1 = 14; % 14 chosen visually as a best fit, Order = N/2.

orderCheb1 = nCheb1-1;
wnCheb1 = [fp1 fp2];
RpCheb1 = (2*(40*log10((1+passTol)/(1-passTol)))); % Peak to Peak Pb R ( = 2*val)
Cheb1 = figure(7);

[bCheb1, aCheb1] = cheby1(orderCheb1, RpCheb1, wnCheb1, 'bandpass');
freqz(bCheb1, aCheb1, 512)
saveas(Cheb1, 'ChebyshevType1.png');

 

% Chebyshev Type II (Cheb2)
nCheb2 = 12; % number of samples

orderCheb2 = (nCheb2 - 1);
wsCheb2 = [fp1 fp2];
RsCheb2 = (2 *(-20 * log10(passTol))); % stopband atten.

Cheb2 = figure(8);
[bCheb2, aCheb2] = cheby2(orderCheb2, RsCheb2, wsCheb2, 'bandpass');
freqz(bCheb2, aCheb2, 512)
saveas(Cheb2, 'ChebyshevType2.png');
 

% Elliptic (E)
nE = 12; % 12 chosen visually as a best fit.

orderE = (nE - 1);
wsE = [fp1 fp2];
RpE = (2*(40*log10((1+passTol)/(1-passTol)))); % Peak to Peak Pb R ( = 2*val)
RsE = (2 *(-20 * log10(passTol))); % stopband atten.

E = figure(9);
[bE, aE] = ellip(orderE, RpE, RsE, wsE, 'bandpass');
freqz(bE, aE, 512)
saveas(E, 'Elliptic.png');

