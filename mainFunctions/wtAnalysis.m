% Compute the wavelet transform using a filterbank cauchy wavelet with a
% 0.7 scale parameter :
% Enders et al. (2012) : Analysis of damped tissue vibrations in
% time-frequency space-A wavelet-based approach
% (http://dx.doi.org/10.1016/j.jbiomech.2012.08.027)


% Compute the wavelet power and the sum of wavelet power (.amplitude)
% estimate damping (.damp)
% for one or several axes (.sep{axis}) and the norm (.norm)

%% INPUTS:
% Acceleration signal in m/s² in column (time,axes)

%% OPTIONS
% see addParameter section


%% OUTPUTS :
% one structure with amplitude and damping estimation


function wtParam=wtAnalysis(acc,varargin)
p = inputParser;
addParameter(p,'Fs',1000,@isnumeric); % sample frequency (Hz)
addParameter(p,'preImpact',[],@isnumeric); % pre impact time (default = start of the signal)
addParameter(p,'postImpact',[],@isnumeric); % total time analyzed (default = end of signal)
addParameter(p,'infFreq',[],@isnumeric); % Freq min analyzed (default = min of wavelets)
addParameter(p,'supFreq',[],@isnumeric); % Freq max analyzed (default = max of wavelets)
addParameter(p,'plotFig',0,@isnumeric); % 1 to plot figure
addParameter(p,'newFig',1,@isnumeric); % create new fig
addParameter(p,'damping',1,@isnumeric); % 0 to not perform damping estimation
addParameter(p,'delay2peak',0.1,@isnumeric); % search maximal peak up to this value in sec, or between 2 values [startSearch endSearch]

parse(p,varargin{:});
Fs=p.Results.Fs;
preImpact=p.Results.preImpact;
postImpact=p.Results.postImpact;
infFreq=p.Results.infFreq;
supFreq=p.Results.supFreq;
plotFig=p.Results.plotFig;
newFig=p.Results.newFig;
damping=p.Results.damping;
delay2peak=p.Results.delay2peak;
%% Transpose if not the right dimension
acc=transposeColmunIfNot(acc);

%% Cut and mirror the signal
[~,preImpactPoints,postImpactPoints,~,~]=defineTime(acc,Fs,Fs,preImpact,postImpact);
acc=acc(preImpactPoints:preImpactPoints+postImpactPoints-1,:);
originalSize=size(acc,1);
paddingSize=2^(nextpow2(2*originalSize));
acc=[zeros(paddingSize/2-originalSize,size(acc,2)); acc(end:-1:1,:); acc; zeros(paddingSize/2-originalSize,size(acc,2))]; % reflect the signal

%% Define the wavelets
scale=0.7; % scale parameter
nWavelets=100; % number of wavelets
Le=size(acc,1)/2;
dt=1/Fs; % time resolution
df=1/(Le*dt); % frequency resolution
f= [0:df:Fs]; % frequency vector

%% Central frequencies
j=2:nWavelets+2;
q=1.45;
r=1.959;
cfs=1/scale*(q+j-2).^r; % central frequencies with q=1.45 and r=1.959

% Nyquist theorem
cfs(cfs>Fs/2)=[];
% Frequency limits
if ~isempty(infFreq)
    cfs(cfs<infFreq)=[];
end
if ~isempty(supFreq)
    cfs(cfs>supFreq)=[];
end

%% Wavelets computations
for j=1:numel(cfs)
    for k=1:numel(f)
        waveFreq(k,j)=((f(1,k)/cfs(1,j))^(cfs(1,j)*scale))*exp(-f(1,k)/cfs(1,j)+1)^(cfs(1,j)*scale);
    end
    waveRealTime(:,j)=fft_real_vvt(waveFreq(:,j),1);
    waveImgTime(:,j)=fft_real_vvt(i*waveFreq(:,j),1);
end

% Timeoral representation mirrored
half=length(waveRealTime)/2;
P1=waveRealTime(1:half,:);
P2=waveRealTime(half+1:end,:);
waveRealTimesMir=cat(1,P2,P1);
half=length(waveImgTime)/2;
P1=waveImgTime(1:half,:);
P2=waveImgTime(half+1:end,:);
waveImgTimesMir=cat(1,P2,P1);


% Frequency representation
for j=1:numel(cfs)
    waveRealFreq(:,j)=real(fft_real_vvt(waveRealTime(:,j),0));
    waveImgFreq(:,j)=fft_real_vvt(waveImgTime(:,j),0);
end

%% WT analysis
if plotFig==1
    if newFig==1
        figure('units','normalized','outerposition',[0 0 1 1],'visible','on')
    end
end

power4norm=0;
for k=1:size(acc,2) % each axis
    
    fftSignal=fft_real_vvt(acc(:,k),0);
    
    for j=1:numel(cfs)
        freqConv(:,j)=waveRealFreq(:,j).*fftSignal; %convolution in frequency space
        tempConv(:,j)=fft_real_vvt(freqConv(:,j),1);
        power(:,j)=tempConv(:,j).^2;
    end
    
    %% Damping calculation
    % define enveloppe for damping
    for j=1:numel(cfs)
        accEnv(:,j)=envelope(power(:,j),round(cfs(j)/2),'peak');
    end
    
    power4norm=power4norm+power.^2;
    if size(accEnv,1)/2+originalSize>size(accEnv,1)
        acc4damp=accEnv(size(accEnv,1)/2+1:end,:);
        power=power(size(accEnv,1)/2+1:end,:);
        accPlot=acc(size(accEnv,1)/2+1:end,:);
    else
        acc4damp=accEnv(size(accEnv,1)/2+1:(size(accEnv,1)/2)+originalSize,:);
        power=power(size(accEnv,1)/2+1:(size(accEnv,1)/2)+originalSize,:);
        accPlot=acc(size(accEnv,1)/2+1:(size(accEnv,1)/2)+originalSize,:);
    end
    
    if plotFig==1
        time=1/Fs:1/Fs:size(acc4damp,1)/Fs;
        if size(acc,2)>1
            subplot(3,size(accPlot,2)+1,k)
        else
            subplot(3,size(accPlot,2),k)
        end
        plot(time,accPlot(:,k),'k')
        title(['Acceleration of axis # ' num2str(k)])
        xlabel('Time (s)');
        ylabel('Amplitude (m\cdots^-^2)')
        box off
        
        if size(acc,2)>1
            subplot(3,size(acc,2)+1,size(acc,2)+1+k)
        else
            subplot(3,size(accPlot,2),k+1)
        end
        plot(time,sum(power,2),'k'); hold on
        plot(time,power,'k-.')
        title('Wavelet power')
        xlabel('Time (s)');
        ylabel('Amplitude (a.u)')
        legend('sum of wavelets power','wavelets power','box','off')
        box off
        
        if size(acc,2)>1
            subplot(3,size(acc,2)+1,2*size(acc,2)+2+k)
        else
            subplot(3,size(accPlot,2),k+2)
        end
    end
    
    % damping estimation
    if size(acc,2)>1
        wtParam.sep{k}.amplitude.sep=power;
        wtParam.sep{k}.amplitude.norm=sum(power,2);
        wtParam.sep{k}.cfs=cfs;
        wtParam.sep{k}.damp=dampingEstimation(acc4damp,'Fs',Fs,'plotfig',plotFig,'delay2peak',delay2peak,'isIMF',1,'titleFig','Damping estimation from signal envelop','units','a.u');
        wtParam.sep{k}.freq=fftAnalysis(power,'Fs',Fs,'isIMF',1,'infFreq',infFreq,'supFreq',supFreq);
    else
        wtParam.sep.amplitude.sep=power;
        wtParam.sep.amplitude.norm=sum(power,2);
        wtParam.sep.cfs=cfs;
        wtParam.sep.damp=dampingEstimation(acc4damp,'Fs',Fs,'plotfig',plotFig,'delay2peak',delay2peak,'isIMF',1,'titleFig','Damping estimation from signal envelop','units','a.u');
        wtParam.sep.freq=fftAnalysis(power,'Fs',Fs,'isIMF',1,'infFreq',infFreq,'supFreq',supFreq);
    end
    clear power
end

%% norm
if size(acc,2)>1
    power=sqrt(power4norm);
    
    for j=1:numel(cfs)
        accEnv(:,j)=envelope(power(:,j),round(cfs(j)/2),'peak');
    end
    
    if size(accEnv,1)/2+originalSize>size(accEnv,1)
        acc4damp=accEnv(size(accEnv,1)/2+1:end,:);
        power=power(size(accEnv,1)/2+1:end,:);
    else
        acc4damp=accEnv(size(accEnv,1)/2+1:(size(accEnv,1)/2)+originalSize,:);
        power=power(size(accEnv,1)/2+1:(size(accEnv,1)/2)+originalSize,:);
    end
    
    if plotFig==1
        time=1/Fs:1/Fs:size(acc4damp,1)/Fs;
        subplot(3,size(accPlot,2)+1,size(accPlot,2)+1)
        plot(time,signalNorm(accPlot),'k')
        title('Norm of acceleration')
        xlabel('Time (s)');
        ylabel('Amplitude (m\cdots^-^2)')
        box off
        
        subplot(3,size(acc,2)+1,2*(size(accPlot,2)+1))
        plot(time,sum(power,2),'k'); hold on
        plot(time,power,'k-.')
        title('Norm of wavelets power')
        xlabel('Time (s)');
        ylabel('Amplitude (a.u)')
        legend('sum of wavelets power','wavelets power','box','off')
        box off
        
        subplot(3,size(acc,2)+1,3*(size(accPlot,2)+1))
    end
    
    wtParam.norm.amplitude.sep=power;
    wtParam.norm.amplitude.norm=sum(power,2);
    wtParam.norm.cfs=cfs;
    wtParam.norm.damp=dampingEstimation(acc4damp,'Fs',Fs,'plotfig',plotFig,'delay2peak',delay2peak,'isIMF',1,'titleFig','Damping estimation from signal envelop','units','a.u');
    wtParam.norm.freq=fftAnalysis(power,'Fs',Fs,'isIMF',1,'infFreq',infFreq,'supFreq',supFreq);
    
end
