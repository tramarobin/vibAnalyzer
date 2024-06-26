% compute normalized DFT, DFT and PSD from FFT algorithm from FFT
% if there is serveral signals, the algorithm calcul the spetrums for each
% axis (.sep) independantly, and also calcul the norm of of the spectrums
% (.norm)
% estimates the frequency parameters for each spectrum

%% INPUT
% acceleration signal in m/s² in column (time,axes)

%% OPTIONS
% see addParameter below

%% OUPUTS
% one structure containing normalized FT, FT and PSD for each axis and
% the norm

function fftParam=fftAnalysis(acc,varargin)

p = inputParser;
addParameter(p,'infFreq',[],@isnumeric); % Freq min analyzed
addParameter(p,'supFreq',[],@isnumeric); % Freq max analyzed
addParameter(p,'padding',[],@isnumeric); % number of point for padding
addParameter(p,'Fs',2000,@isnumeric); % samplefrequency
addParameter(p,'plotFig',0,@isnumeric); % if 1, plot figure
addParameter(p,'newFig',1,@isnumeric); % if 1, plot new figure
addParameter(p,'isIMF',0,@isnumeric); % 1 sum the spectrums if IMF
addParameter(p,'preImpact',[],@isnumeric); % pre impact time (default = start of the signal)
addParameter(p,'postImpact',[],@isnumeric); % total time analyzed (default = end of signal)

parse(p,varargin{:});
infFreq=p.Results.infFreq;
supFreq=p.Results.supFreq;
padding=p.Results.padding;
Fs=p.Results.Fs;
plotFig=p.Results.plotFig;
newFig=p.Results.newFig;
isIMF=p.Results.isIMF;
preImpact=p.Results.preImpact;
postImpact=p.Results.postImpact;


acc=transposeColmunIfNot(acc);
[~,preImpactPoints,postImpactPoints,~,~,acc]=defineTime(acc,Fs,Fs,preImpact,postImpact,0);
acc=acc(preImpactPoints:preImpactPoints+postImpactPoints-1,:);

%% padding
if isempty(padding)
    padding=2^(nextpow2(size(acc,1))+2); % pad to next+2 power of 2
end

accPad=[acc ;zeros(padding-size(acc,1),size(acc,2))];

L=size(accPad,1);
Y = fft(accPad);
spectralWidth=Fs/L;
f = Fs*(0:(L/2))/L;

%% find freq analyzed
if isempty(infFreq)
    startF=1;
else
    startF=find(f<infFreq,1,'last');
end
if isempty(supFreq)
    stopF=numel(f);
else
    stopF=find(f>supFreq,1);
end
if isempty(stopF)
    stopF=numel(f);
end
if isempty(startF)
    startF=1;
end

f=f(startF:stopF);

%% FFT
% amplitude spectrum
P2 = abs(Y/L);
P1 = P2(1:L/2+1,:);
P1(2:end-1,:) = 2*P1(2:end-1,:);
FT=P1(startF:stopF,:);
nFT=FT/spectralWidth;

% power spectrum
S2 = (1/(Fs*L)) * abs(Y).^2;
S1 = S2(1:L/2+1,:);
S1(2:end-1) = 2*S1(2:end-1);
PSD=S1(startF:stopF,:);


if plotFig==1
    if newFig==1
        figure('units','normalized','outerposition',[0 0 1 1],'visible','on')
    end
    subplot(121)
    plot(1/Fs:1/Fs:size(acc,1)/Fs,acc,'k')
    xlabel('Time (s)')
    ylabel('Acceleration (m\cdots^-^2)')
    title('Acceleration in time domain')
    box off
    
    subplot(122)
end
fftParam.normalizedFT=frequencyEstimation(f,nFT,'plotFig',plotFig,'titleFig','normalized FT','isIMF',isIMF);


fftParam.PSD=frequencyEstimation(f,PSD,'units','m^2\cdots^-^4/Hz','titleFig','PSD','isIMF',isIMF);

fftParam.FT=frequencyEstimation(f,FT,'isIMF',isIMF);

fftParam.spectralWidth=spectralWidth;

end


