% Compute the continuous wavelet transform using a filterbank of morse
% (3,60) wavelets with 48 voices per octave
% the signal, if you want avoid reflexion,
% Compute the map of coefficients (.map)
% estimate frequency parameters (.freq) and damping (.damp)
% for one or several axes (.sep(:,axis)) and the norm (.norm)


%% INPUTS:
% Acceleration signal in m/s² in column (time,axes)

%% OPTIONS
% see addParameter section


%% OUTPUTS :
% one structure with coeffcient maps, frequency analysis and damping


function cwtParam= cwtAnalysis(acc,varargin)
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
addParameter(p,'newFig',1,@isnumeric); % create new fig
addParameter(p,'isIMF',0,@isnumeric); % sum the maps instead of norm if IMF
addParameter(p,'reflection',[],@isnumeric); % 1 use reflection at the start of the signal and add 0 padding of 2 seconds on each side centered on heel strike, it improve mode separation and allow to investigate lower frequencies, /!\ the signal analyzed is not the one you measured anymore. Enders et al. 2012; http://dx.doi.org/10.1016/j.jbiomech.2012.08.027
addParameter(p,'damping',1,@isnumeric); % 0 to not perform damping estimation
addParameter(p,'delay2peak',[]); % search maximal peak up to this value in sec, or between 2 values [startSearch endSearch]

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
newFig=p.Results.newFig;
colorbarOk=p.Results.colorbarOk;
isIMF=p.Results.isIMF;
reflection=p.Results.reflection;
damping=p.Results.damping;
delay2peak=p.Results.delay2peak;
%% Transpose if not the right dimension
acc=transposeColmunIfNot(acc);

%% reflection if signal cut at the impact
if isempty(reflection)
    if isempty(preImpact)
        reflection=1;
    else
        reflection=0;
    end
end

%% Time analyzed
[newFs,preImpactPoints,postImpactPoints,newPreImpactPoints,newPostImpactPoints,acc]=defineTime(acc,Fs,newFs,preImpact,postImpact,reflection);

%% Filterbank and frequency range analyzed
[fb,f,infLim,supLim,infFreq,supFreq,intFreq,indIntInfLim,indIntSupLim]=defineFilterBank(acc,Fs,interpFreq,fb,infFreq,supFreq,reflection);

%% Wavelet Transform
coefficients=waveletAnalysis(acc,fb,'reflection',reflection);

%% Map Interpolation
for i=1:size(acc,2) % each axis
    intCoefficients{i}=interpCoefficients(coefficients.axisCoefficients{i}(supLim:infLim,:),f(supLim:infLim,:),intFreq,interpFreq,newFs,Fs,indIntInfLim,indIntSupLim,newPreImpactPoints,newPostImpactPoints);
end
if size(acc,2)>1
    intNormCoefficients=interpCoefficients(coefficients.normCoefficients(supLim:infLim,:),f(supLim:infLim,:),intFreq,interpFreq,newFs,Fs,indIntInfLim,indIntSupLim,newPreImpactPoints,newPostImpactPoints);
end

cwtParam.map.sep.amplitude=intCoefficients;
if size(acc,2)>1
    cwtParam.map.norm.amplitude=intNormCoefficients;
end
cwtParam.map.dimensions=size(intCoefficients{1});
cwtParam.map.fs=newFs;
cwtParam.map.f=intFreq;


%% delete NaN coefficient and associated frequencies
nanVal=find(isnan(mean(intCoefficients{1}')));
intFreq(nanVal)=[];
for i=1:size(acc,2)
    intCoefficients{i}(nanVal,:)=[];
    if i>2
    intNormCoefficients(nanVal,:)=[];
    end
end


%% Integration of time and frequency
for i=1:size(acc,2)
    energy(:,i)=trapz(intCoefficients{i}')/newFs;
end

for i=1:size(acc,2)
    power(:,i)=trapz(intFreq,intCoefficients{i})';
end


%% PLOT and variables extraction
if plotFig==1
    if newFig==1
        figure('units','normalized','outerposition',[0 0 1 1],'visible','on')
    end
end

% Map of coefficients
if plotFig==1
    newTime=1/newFs:1/newFs:newPostImpactPoints/newFs;
    time=1/Fs:1/Fs:postImpactPoints/Fs;
    
    subplot(221)
    plot(time,acc(preImpactPoints:preImpactPoints+postImpactPoints-1,:),'k')
    title('Acceleration signal in time domain')
    xlabel('Time (s)')
    ylabel('Acceleration (m\cdots^-^2)')
    yLimits = get(gca, 'YLim');
    ylim(1.1*yLimits);
    box off
    
    subplot(223)
    if size(acc,2)==1
        displayCoeffMap(intCoefficients{1},intFreq,newTime,'colorbarOk',colorbarOk,'newFig',0,'titleFig','Amplitude in the time-frequency domain')
    else
        displayCoeffMap(intNormCoefficients,intFreq,newTime,'colorbarOk',colorbarOk,'newFig',0,'titleFig','Amplitude in the time-frequency domain')
    end
end

% Frequencial analysis
if plotFig==1
    subplot(222)
end
cwtParam.freq=frequencyEstimation(intFreq,energy,'plotfig',plotFig,'isIMF',isIMF,'titleFig','Amplitude in the frequency domain');

% Time domain analysis (damping)
if plotFig==1
    subplot(224)
end
cwtParam.damp=dampingEstimation(power,'calculDamping',damping,'Fs',newFs,'plotfig',plotFig,'isIMF',isIMF,'delay2peak',delay2peak,'titleFig','Amplitude in the time domain');


end