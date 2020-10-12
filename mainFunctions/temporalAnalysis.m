% Trama Robin, LIMB, Lyon, France
% compute several parameters with a temporal analysis such as peak
% acceleration, time to peak, estimated frequency, damping calculated from
% signal envelope
% if there is serveral signals, the algorithm calcul the parameters for each
% axis independantly, and also for the norm

%% INPUT
% acceleration signal in m/s² in column (time,axes)

%% OPTIONS
% see addParameter below

%% OUPUTS
% one structure containing the parameters for each axis (.sep)
% and the norm
% Contains also an estimation of damping .damp.sep for each axis .damp.norm
% for the norm


function temporalParam=temporalAnalysis(acc,varargin)

p = inputParser;
addParameter(p,'Fs',1000,@isnumeric); % samplefrequency
addParameter(p,'plotFig',0,@isnumeric); % if 1, plot figure
addParameter(p,'newFig',0,@isnumeric); % if 1, plot in a new figure
addParameter(p,'isIMF',0,@isnumeric); % 1 sum the IMFs instead of norm if IMFs (in .norm)
addParameter(p,'preImpact',[],@isnumeric); % pre impact time (default = start of the signal)
addParameter(p,'postImpact',[],@isnumeric); % total time analyzed (default = end of signal)
parse(p,varargin{:});
Fs=p.Results.Fs;
plotFig=p.Results.plotFig;
newFig=p.Results.newFig;
isIMF=p.Results.isIMF;
preImpact=p.Results.preImpact;
postImpact=p.Results.postImpact;

acc=transposeColmunIfNot(acc);
[~,preImpactPoints,postImpactPoints,~,~,acc]=defineTime(acc,Fs,Fs,preImpact,postImpact,0);
acc=acc(preImpactPoints:preImpactPoints+postImpactPoints-1,:);

temporalParam.sampleFrequency=Fs;

%% Separated signals
for i=1:size(acc,2)
    [maxAcc(i),tMaxAcc(i)]=max(acc(:,i));
    [minAcc(i),tMinAcc(i)]=min(acc(:,i));
    [peakAcc(i),time2peak(i)]=max(abs(acc(:,i)));
    peak2peak(i)=maxAcc(i)-minAcc(i);
    
    tPeaks(i,:)=sort([tMaxAcc(i) tMinAcc(i)]);
    tp2p(i)=diff(tPeaks(i,:))/Fs;
    estimatedFrequency(i)=1/(2*tp2p(i)); % Boyer et al. 2004; doi:10.1016/j.jbiomech.2004.01.002
    
    integratedAmplitude(i)=trapz(abs(acc(:,i)))/Fs;
    
    rmsAcc(i)=mean(rmsWindow(signalNorm(acc(:,i)),0.01,Fs)); % Ehrstrom et al. 2018; 10.3389/fphys.2018.01627
    standardDeviation(i)=std(acc(:,i)); % Gellaerts et al. (2017); 10.23736/S0022-4707.16.06721-9
    
end

temporalParam.sep.max=maxAcc;
temporalParam.sep.min=minAcc;
temporalParam.sep.peak=peakAcc;
temporalParam.sep.peak2peak=peak2peak;
temporalParam.sep.timePeak2peak=tp2p;
temporalParam.sep.estimatedFrequency=estimatedFrequency;
temporalParam.sep.integratedAmplitude=integratedAmplitude;
temporalParam.sep.rms=rmsAcc;
temporalParam.sep.standardDeviation=standardDeviation;
temporalParam.sep.time2max=tMaxAcc/Fs;
temporalParam.sep.time2min=tMinAcc/Fs;
temporalParam.sep.time2peak=time2peak/Fs;
temporalParam.sep.ind.max=tMaxAcc;
temporalParam.sep.ind.min=tMinAcc;
temporalParam.sep.ind.peak=time2peak;

%% Norm of signal
if size(acc,2)>1
    if isIMF==0
        nAcc=signalNorm(acc);
    else
        nAcc=sum(acc,2);
    end
    [nPeakAcc,nTime2peak]=max(nAcc);
    nIintegratedAmplitude=trapz(nAcc)/Fs;
    
    nRmsAcc=mean(rmsWindow(nAcc,0.01,Fs)); % Ehrstrom et al. 2018; 10.3389/fphys.2018.01627
    nStandardDeviation=std(nAcc); % Gellaerts et al. (2017); 10.23736/S0022-4707.16.06721-9
    
    
    temporalParam.norm.peak=nPeakAcc;
    temporalParam.norm.time2peak=nTime2peak/Fs;
    temporalParam.norm.integratedAmplitude=nIintegratedAmplitude;
    temporalParam.norm.rms=nRmsAcc;
    temporalParam.norm.standardDeviation=nStandardDeviation;
    temporalParam.norm.ind2peak=nTime2peak;
    
    
end

%% Damping
accEnv=envelope(abs(acc),round(mean(0.5./estimatedFrequency*Fs)),'peak');
temporalParam.damp=dampingEstimation(accEnv(min(tPeaks(:,1)):end,:),'Fs',Fs,'isIMF',1,'delay2peak',0.1);


%% PLOT
if plotFig==1
    if newFig==1
        figure('units','normalized','outerposition',[0 0 1 1],'visible','on')
    end
    
    time=1/Fs:1/Fs:size(acc,1)/Fs;
    
    if size(acc,2)>1 % several signals
        
        plot(time,nAcc,'k'); hold on
        plot(time,acc(:,1),':');
        plot(time,acc(:,2),'k--');
        if size(acc,2)>2
            plot(time,acc(:,3),'k-.');
        end
        
        scatter(nTime2peak/Fs,nPeakAcc,'k+','HandleVisibility','off');
        text(nTime2peak/Fs+0.02*max(time),nPeakAcc,['Peak acceleration = ' sprintf('%0.1f',nPeakAcc) ' m\cdots^-^2 at t = ' sprintf('%0.3f',nTime2peak/Fs) ' s'])
        
        if size(acc,2)==2
            legend({'norm','x','y'},'box','off')
        else
            legend({'norm','x','y','z'},'box','off')
        end
        title('Acceleration signals in temporal domain')
        
    else % one signal : plot the signal analysis
        plot(time,acc,'k','HandleVisibility','on'); hold on
        scatter(tMaxAcc/Fs,maxAcc,'k+','HandleVisibility','off');
        scatter(tMinAcc/Fs,minAcc,'kx','HandleVisibility','off');
        text(tMinAcc/Fs+0.02*max(time),minAcc,['Peak deceleration = ' sprintf('%0.1f',minAcc) ' m\cdots^-^2 at t = ' sprintf('%0.3f',tMinAcc/Fs) ' s'])
        text(tMaxAcc/Fs+0.02*max(time),maxAcc,['Peak acceleration = ' sprintf('%0.1f',maxAcc) ' m\cdots^-^2 at t = ' sprintf('%0.3f',tMaxAcc/Fs) ' s'])
        if tMaxAcc<tMinAcc
            text(1.5*max([tMaxAcc tMinAcc])/Fs,maxAcc-0.1*peak2peak,['Peak to peak amplitude = ' sprintf('%0.1f',peak2peak) ' m\cdots^-^2'])
            text(1.5*max([tMaxAcc tMinAcc])/Fs,maxAcc-0.15*peak2peak,['Time between peaks = ' sprintf('%0.3f',tp2p) ' s'])
            text(1.5*max([tMaxAcc tMinAcc])/Fs,maxAcc-0.2*peak2peak,['Estimated frequency = ' sprintf('%0.1f',estimatedFrequency) ' Hz'])
        else
            text(1.5*max([tMaxAcc tMinAcc])/Fs,minAcc+0.1*peak2peak,['Peak to peak amplitude = ' sprintf('%0.1f',peak2peak) ' m\cdots^-^2'])
            text(1.5*max([tMaxAcc tMinAcc])/Fs,minAcc+0.15*peak2peak,['Time between peaks = ' sprintf('%0.3f',tp2p) ' s'])
            text(1.5*max([tMaxAcc tMinAcc])/Fs,minAcc+0.2*peak2peak,['Estimated frequency = ' sprintf('%0.1f',estimatedFrequency) ' Hz'])
        end
        
        
        title('Acceleration signal in time domain')
        
    end
    
    box off
    xlabel('Time (s)')
    ylabel('Acceleration  (m\cdots^-^2)')
    yLimits = get(gca, 'YLim');
    ylim(1.1*yLimits);
    
end
end









