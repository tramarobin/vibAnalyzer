% calculate the damping coefficient
% the maximal power decrease and the set powerLimits (default =10%)
% using optimisation function


%% INPUT
% Frequency-integrated wavelet coefficients (overall power, normalized time
% power) in column (time,axes)

%% OPTIONS
% see addParameter below

%% OUPUTS
% one structure containing the parameters for each axis and the norm
% damping coefficient, settling time, maximal amplitude

function dampingParam=dampingEstimation(power,varargin)

p = inputParser;
addParameter(p,'Fs',1000,@isnumeric); % samplefrequency
addParameter(p,'powerLimits',0.1,@isnumeric); % Percentage of power at which the LSminimisation stop
addParameter(p,'delay2peak',[],@isnumeric); % search maximal peak up to this value in sec, or between 2 values [startSearch endSearch]
addParameter(p,'plotFig',0,@isnumeric); % if 1, plot figure (for one axis, or the norm)
addParameter(p,'newFig',0,@isnumeric); % if 1, plot new figure
addParameter(p,'isIMF',0,@isnumeric); % 1 sum the spectrums if IMF
addParameter(p,'titleFig',[],@ischar); % title fig
addParameter(p,'units','m\cdots^-^2/s',@ischar); % units

parse(p,varargin{:});
Fs=p.Results.Fs;
powerLimits=p.Results.powerLimits;
delay2peak=p.Results.delay2peak;
plotFig=p.Results.plotFig;
newFig=p.Results.newFig;
isIMF=p.Results.isIMF;
units=p.Results.units;
titleFig=p.Results.titleFig;

powers=transposeColmunIfNot(power);
dampingParam.sep.amplitude=powers;

%% SEP
for i=1:size(power,2)
    
    power=powers(:,i);
    
    if isempty(delay2peak)
        [Maximal_power,maxPos]=max(power);
    elseif numel(delay2peak)==1
        [Maximal_power,maxPos]=max(power(1:delay2peak*Fs));
    elseif numel(delay2peak)==2
        [Maximal_power,maxPos]=max(power(delay2peak1(1)*Fs:delay2peak(2)*Fs));
        maxPos=maxPos+delay2peak1(1)*Fs-1;
    end
    
    Diff_Puissance=diff(power);
    if maxPos+round(50*Fs/1000)<numel(Diff_Puissance)
        [~,Damp_start]=min(Diff_Puissance(maxPos:maxPos+round(50*Fs/1000)));
    else
        [~,Damp_start]=min(Diff_Puissance(maxPos:end));
    end
    Damp_start=Damp_start+maxPos-1;
    Damp_seuil=find(power(Damp_start:end)<powerLimits*Maximal_power,1)+Damp_start-1;
    Damp_diff_d=find(Diff_Puissance(Damp_start:end)>0.001*power(Damp_start));
    Damp_diff_s=find(power(Damp_start:end)<0.5*power(Damp_start));
    Damp_diff=intersect(Damp_diff_d,Damp_diff_s);
    if isempty(Damp_diff)==0
        Damp_diff=Damp_diff(1)+Damp_start-1;
    else
        Damp_diff=size(power,1);
    end
    if isempty(Damp_seuil) && isempty(Damp_diff)
        Damp_end=size(power,1);
    else
        Damp_end=min([Damp_seuil Damp_diff]);
    end
    settlingTime=(Damp_end-maxPos)/Fs;
    time_max=(Damp_end-Damp_start)/Fs;
    t=0:1/Fs:time_max;
    Pdec=power(Damp_start:Damp_end)';
    options = optimoptions('fmincon','Display','off');
    if settlingTime>0
        [Damping_property,Damp_err]=fmincon(@(d)dampingMinimisation(Pdec,t,d),10,[],[],[],[],0,1000,[],options);
    else
        Damping_property=nan; Damp_err=nan; Damp_start=nan; Damp_end=nan;
    end
    dampingParam.sep.dampingCoefficient(i)=Damping_property;
    dampingParam.sep.settlingTime(i)=settlingTime;
    dampingParam.sep.maxAmplitude(i)=Maximal_power;
    dampingParam.sep.totalAmplitude(i)=trapz(power)/Fs;
    dampingParam.sep.indices.maxPos(i)=maxPos;
    dampingParam.sep.indices.dampStart(i)=Damp_start;
    dampingParam.sep.indices.dampEnd(i)=Damp_end;
    dampingParam.sep.temp.maxPos(i)=maxPos/Fs;
    dampingParam.sep.temp.dampStart(i)=Damp_start/Fs;
    dampingParam.sep.temp.dampEnd(i)=Damp_end/Fs;
    dampingParam.sep.totalError(i)=Damp_err;
    
end


%% NORM
if size(powers,2)>1
    if isIMF==0
        power=signalNorm(powers);
    else
        power=sum(powers,2);
    end
    
    if isempty(delay2peak)
        [Maximal_power,maxPos]=max(power);
    elseif numel(delay2peak)==1
        [Maximal_power,maxPos]=max(power(1:delay2peak*Fs));
    elseif numel(delay2peak)==2
        [Maximal_power,maxPos]=max(power(delay2peak1(1)*Fs:delay2peak(2)*Fs));
        maxPos=maxPos+delay2peak1(1)*Fs-1;
    end
    
    Diff_Puissance=diff(power);
    if maxPos+round(50*Fs/1000)<numel(Diff_Puissance)
        [~,Damp_start]=min(Diff_Puissance(maxPos:maxPos+round(50*Fs/1000)));
    else
        [~,Damp_start]=min(Diff_Puissance(maxPos:end));
    end
    Damp_start=Damp_start+maxPos-1;
    Damp_seuil=find(power(Damp_start:end)<powerLimits*Maximal_power,1)+Damp_start-1;
    Damp_diff_d=find(Diff_Puissance(Damp_start:end)>0.001*power(Damp_start));
    Damp_diff_s=find(power(Damp_start:end)<0.5*power(Damp_start));
    Damp_diff=intersect(Damp_diff_d,Damp_diff_s);
    if isempty(Damp_diff)==0
        Damp_diff=Damp_diff(1)+Damp_start-1;
    else
        Damp_diff=size(power,1);
    end
    if isempty(Damp_seuil) && isempty(Damp_diff)
        Damp_end=size(power,1);
    else
        Damp_end=min([Damp_seuil Damp_diff]);
    end
    settlingTime=(Damp_end-maxPos)/Fs;
    time_max=(Damp_end-Damp_start)/Fs;
    t=0:1/Fs:time_max;
    Pdec=power(Damp_start:Damp_end)';
    options = optimoptions('fmincon','Display','off');
    [Damping_property,Damp_err]=fmincon(@(d)dampingMinimisation(Pdec,t,d),10,[],[],[],[],0,1000,[],options);
    
    dampingParam.norm.amplitude=power;
    dampingParam.norm.dampingCoefficient=Damping_property;
    dampingParam.norm.settlingTime=settlingTime;
    dampingParam.norm.maxAmplitude=Maximal_power;
    dampingParam.norm.totalAmplitude=sum(power);
    dampingParam.norm.indices.maxPos=maxPos;
    dampingParam.norm.indices.dampStart=Damp_start;
    dampingParam.norm.indices.dampEnd=Damp_end;
    dampingParam.norm.temp.maxPos=maxPos/Fs;
    dampingParam.norm.temp.dampStart=Damp_start/Fs;
    dampingParam.norm.temp.dampEnd=Damp_end/Fs;
    dampingParam.norm.totalError=Damp_err;
    
end

dampingParam.samplingFrequency=Fs;

%% PLOT
if plotFig==1
    if newFig==1
        figure('visible','on')
    end
    
    time=1/Fs:1/Fs:size(power,1)/Fs;
    if size(powers,2)>1 & size(powers,2)<4
        plot(time,power,'k','HandleVisibility','on'); hold on
        plot(time,powers(:,1),'-.k','HandleVisibility','on');
        plot(time,powers(:,2),'-.k','HandleVisibility','off');
        if size(powers,2)>2
            plot(time,powers(:,3),'-.k','HandleVisibility','off');
        end
    elseif size(powers,2)==1
        plot(time,power,'k','HandleVisibility','off'); hold on
    else
        plot(time,power,'k','HandleVisibility','on'); hold on
        plot(time,powers(:,1),'-.k','HandleVisibility','on');
        plot(time,powers(:,2:end),'-.k','HandleVisibility','off');
    end
    scatter(time(maxPos),Maximal_power,'k+','HandleVisibility','off')
    if Damp_start<numel(time)
        scatter(time(Damp_start),power(Damp_start),'kv')
        if Damp_end<=size(time,2)
            scatter(time(Damp_end),power(Damp_end),'ko')
            plot(time(Damp_start:Damp_end),Pdec(1).*exp(-Damping_property*((0:size(Damp_start:Damp_end,2)-1)/Fs)),'k--','Linewidth',2)
        else
            plot(time(Damp_start:end),Pdec(1).*exp(-Damping_property*((0:size(Damp_start:numel(time),2)-1)/Fs)),'k--','Linewidth',2)
            warning('increase time post impact to see complete damping on the figure')
        end
    end
    
    text(time(maxPos)+0.05*max(time),Maximal_power,['Peak amplitude = ' sprintf('%0.1f',Maximal_power) ' ' units ' at t = ' sprintf('%3.3f',maxPos/Fs) ' s'])
    text(0.02*max(time),0.2*max(power),['Damping coefficient = ' sprintf('%3.1f',Damping_property) ' s^-^1'])
    text(0.02*max(time),0.1*max(power),['Settling time = ' sprintf('%3.3f',settlingTime) ' s'])
    
    xlabel('Time (s)');
    ylabel(['Amplitude (' units ')'])
    title(titleFig)
    if size(powers,2)>1
        if isIMF==0
            legend({'norm of signals','signals from different axes',['Start of damping estimation at t = ' sprintf('%3.3f',Damp_start/Fs) ' s'],['End of damping estimation at t = ' sprintf('%3.3f',Damp_end/Fs) ' s'],'Least-square minimisation'},'box','off')
        else
            legend({'sum of sub signals','sub signals',['Start of damping estimation at t = ' sprintf('%3.3f',Damp_start/Fs) ' s'],['End of damping estimation at t = ' sprintf('%3.3f',Damp_end/Fs) ' s'],'Least-square minimisation'},'box','off')
        end
        
    else
        legend({['Start of damping estimation at t = ' sprintf('%3.3f',Damp_start/Fs) ' s'],['End of damping estimation at t = ' sprintf('%3.3f',Damp_end/Fs) ' s'],'Least-square minimisation'},'box','off')
    end
    legPos=get(legend,'position');
    set(legend,'position',[legPos(1) legPos(2)-0.05 legPos(3:4)])
    box off
    yLimits = get(gca, 'YLim');
    ylim(1.1*yLimits);
end
end
