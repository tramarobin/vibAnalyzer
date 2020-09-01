% Combine different types of analysis to characterize transient vibratory signals
% time (temporalAnalysis), frequencial (fftAnalysis) and time-frequency (cwtAnalysis,wtAnalysis, emdAnalysis, peemdAnalysis)
% see the details for each funtions
% The plot represents the time and cwt analyses
% /!\ require wavelet (for CWT), signal (for EMD), and optimisation (for damping and acc (model) estimation) toolboxes 

%% INPUTS:
% Acceleration signal in m/s² in with the format(time series,axes)

%% OPTIONS
% see addParameter section

%% OUTPUTS :
% one structure with temporal, frequencial and time-frequency parameters

function accParam=vibAnalyzer(acc,varargin)
p = inputParser;
addParameter(p,'Fs',1000,@isnumeric); % sample frequency (Hz)
addParameter(p,'fb',[]); % filterbank (default: morse (3,60))
addParameter(p,'preImpact',[],@isnumeric); % pre impact time (default = start of the signal)
addParameter(p,'postImpact',[],@isnumeric); % total time analyzed (default = end of signal)
addParameter(p,'infFreq',[],@isnumeric); % Freq min analyzed (default = min of fb)
addParameter(p,'supFreq',[],@isnumeric); % Freq max analyzed (default = max du fb)
addParameter(p,'interpFreq',1,@isnumeric); % change interpolation frequency to reduce map size
addParameter(p,'newFs',[],@isnumeric); % change sample frequency to reduce map size (must be < Fs)
addParameter(p,'colorbarOk',1,@isnumeric); % 0 dont plot colorbar
addParameter(p,'plotFig',0,@isnumeric); % 1 to plot figure
addParameter(p,'padding',2048,@isnumeric); % number of point for padding the FFT
addParameter(p,'reflection',0,@isnumeric); % 1 use reflection at the start of the signal and add 0 padding of 2 seconds on each side centered on heel strike, it improve mode separation and allow to investigate lower frequencies, /!\ the signal analyzed is not the one you measured anymore. Enders et al. 2012; http://dx.doi.org/10.1016/j.jbiomech.2012.08.027

parse(p,varargin{:});
Fs=p.Results.Fs;
fb=p.Results.fb;
preImpact=p.Results.preImpact;
postImpact=p.Results.postImpact;
infFreq=p.Results.infFreq;
supFreq=p.Results.supFreq;
interpFreq=p.Results.interpFreq;
newFs=p.Results.newFs;
plotFig=p.Results.plotFig;
colorbarOk=p.Results.colorbarOk;
padding=p.Results.padding;
reflection=p.Results.reflection;

%% Transpose if not the right dimension
acc=transposeColmunIfNot(acc);

%% TEMPORAL
%Time analyzed
[~,preImpactPoints,postImpactPoints,~,~]=defineTime(acc,Fs,newFs,preImpact,postImpact);
accParam.TEMPORAL=temporalAnalysis(acc(preImpactPoints:preImpactPoints+postImpactPoints-1,:),'Fs',Fs);

%% CWT
accParam.CWT=cwtAnalysis(acc,'fb',fb,'Fs',Fs,'infFreq',infFreq,'supFreq',supFreq,'reflection',reflection,...
    'preImpact',preImpact,'postImpact',postImpact,'interpFreq',interpFreq,'newFs',newFs);

%% FFT
accParam.FFT=fftAnalysis(acc(preImpactPoints:preImpactPoints+postImpactPoints-1,:),'Fs',Fs,'infFreq',infFreq,'supFreq',supFreq,'padding',padding);

%% EMD & PEEMD
accParam.EMD=emdAnalysis(acc(preImpactPoints:preImpactPoints+postImpactPoints-1,:),'Fs',Fs,'infFreq',1,'supFreq',supFreq,'padding',padding,'interpFreq',interpFreq,'newFs',newFs,'reflection',reflection);
accParam.PEEMD=peemdAnalysis(acc(preImpactPoints:preImpactPoints+postImpactPoints-1,:),'Fs',Fs,'infFreq',1,'supFreq',supFreq,'padding',padding,'interpFreq',interpFreq,'newFs',newFs,'reflection',reflection);

%% Estimation
accParam.ESTIM=accEstimation(acc(preImpactPoints:preImpactPoints+postImpactPoints-1,:),'Fs',Fs);

%% WT VVT/Enders
accParam.WT=wtAnalysis(acc(preImpactPoints:preImpactPoints+postImpactPoints-1,:),'Fs',Fs,'infFreq',infFreq,'supFreq',supFreq);

%% PLOT (temportal + CWT)
if plotFig==1
    figure('units','normalized','outerposition',[0 0 1 1],'visible','on')
    
    % Temporal analysis
    subplot(221)
    temporalAnalysis(acc(preImpactPoints:preImpactPoints+postImpactPoints-1,:),'Fs',Fs,'plotFig',1);
    
    % Map of coefficients
    time=1/Fs:1/Fs:postImpactPoints/Fs;
    subplot(223)
    if size(acc,2)==1
        displayCoeffMap(accParam.CWT.map.sep.amplitude{1},accParam.CWT.map.f,time,'colorbarOk',colorbarOk,'newFig',0,'title','Map of coefficients in time-frequency domain')
    else
        displayCoeffMap(accParam.CWT.map.norm.amplitude,accParam.CWT.map.f,time,'colorbarOk',colorbarOk,'newFig',0,'title','Map of coefficients in time-frequency domain')
    end
    
    % Frequencial
    subplot(222)
    frequencyEstimation(accParam.CWT.freq.f,accParam.CWT.freq.sep.amplitude,'plotfig',1,'title','Normalised amplitude in frequency domain');
    
    % Time domain analysis (damping)
    subplot(224)
    dampingEstimation(accParam.CWT.damp.sep.amplitude,'Fs',accParam.CWT.damp.samplingFrequency,'plotfig',1,'title','Normalised amplitude in time domain');
    
end
end
