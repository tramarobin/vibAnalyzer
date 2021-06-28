clear
close all
clc

addpath(genpath(pwd));
load ACC


%% Parameters (see funtion descriptions for other parameters)
Fs=2000; % sample frequency
infFreq=6; % lowest frequency analyzed (set the range for cwt and fft, discart some IMF or wavelets from wt)
supFreq=100;  % highest frequency analyzed (set the range for cwt and fft, discart some IMF or wavelets from wt)
plotFig=1; % 0 to not plot
preImpact=0.2; % the impact occured at 0.2 second after the start of the signal (400 points at 2000 Hz)
postImpact=0.25; % 0.25 seconds are considered post impact

% The signal was band-pass filtered between 6 and 100 Hz
% three axes
% acc=ACC;

% vertical axis
acc=ACC(:,3);

% gather all fonctions
accParam=vibAnalyzer(acc,'Fs',Fs,'infFreq',infFreq,'supfreq',supFreq,'preImpact',preImpact,'postImpact',postImpact);

% cwt with morse wavelets
% default parameter (3,60)
cwtParam=cwtAnalysis(acc,'Fs',Fs,'infFreq',infFreq,'supfreq',supFreq,'preImpact',preImpact,'postImpact',postImpact,'plotfig',plotFig,'newfig',1);

% better frequency resolution (3,120)
fb3120=cwtfilterbank('waveletparameter',[3 120]);
cwtParamFreq=cwtAnalysis(acc,'fb',fb3120,'Fs',Fs,'infFreq',infFreq,'supfreq',supFreq,'preImpact',preImpact,'postImpact',postImpact,'plotfig',plotFig,'newfig',1);

% better time resolution (3,30)
fb330=cwtfilterbank('waveletparameter',[3 30]);
cwtParamTime=cwtAnalysis(acc,'fb',fb330,'Fs',Fs,'infFreq',infFreq,'supfreq',supFreq,'preImpact',preImpact,'postImpact',postImpact,'plotfig',plotFig,'newfig',1);

% wt with cauchy wavelets (VVT/Enders)
wtParam=wtAnalysis(acc,'Fs',Fs,'infFreq',infFreq,'supfreq',supFreq,'preImpact',preImpact,'postImpact',postImpact,'plotfig',plotFig,'newfig',1);

% fft analysis
fftParam=fftAnalysis(acc,'Fs',Fs,'infFreq',infFreq,'supfreq',supFreq,'preImpact',preImpact,'postImpact',postImpact,'plotfig',plotFig,'newfig',1);

% emd analysis
emdParam=emdAnalysis(acc,'Fs',Fs,'infFreq',infFreq,'supfreq',supFreq,'preImpact',preImpact,'postImpact',postImpact,'plotfig',plotFig,'newfig',1);

% peemd analysis (Khassetarash)
peemdParam=peemdAnalysis(acc,'Fs',Fs,'infFreq',infFreq,'supfreq',supFreq,'preImpact',preImpact,'postImpact',postImpact,'plotfig',plotFig,'newfig',1);

% temporal anaylis 
temporalParam=temporalAnalysis(acc,'Fs',Fs,'preImpact',preImpact,'postImpact',postImpact,'plotFig',plotFig,'newFig',1); % might consider filtering the signal as no frequency range are set

preImpact=0.238; % time juste before first peak
% estimation
estimationParam=accEstimation(acc,'Fs',Fs,'preImpact',preImpact,'postImpact',postImpact,'plotFig',plotFig,'newFig',1); % might consider filtering the signal as no frequency range are set

% estimation 2f
estimationParam2f=accEstimation2f(acc,'Fs',Fs,'preImpact',preImpact,'postImpact',postImpact,'plotFig',plotFig,'newFig',1); % might consider filtering the signal as no frequency range are set

% estimation 2fEnders
estimationParam2fEnders=accEstimation2fEnders(acc,'Fs',Fs,'preImpact',preImpact,'postImpact',postImpact,'plotFig',plotFig,'newFig',1); % might consider filtering the signal as no frequency range are set

% estimation OptiF
estimationParamOptiF=accEstimationOptiF(acc,'Fs',Fs,'preImpact',preImpact,'postImpact',postImpact,'plotFig',plotFig,'newFig',1); % might consider filtering the signal as no frequency range are set
