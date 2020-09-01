% calculate the frequency variables from power spectrums
% if sevral spectrums, calculate varaibles for each axis and the norm

%% INPUTS
% power spectrums (frequency,axes)
% frequencies

%% OPTION
% see addParameter below

%% OUTPUTS
% one structure containing the parameters for each axis and the norm
% Main, median, mean frequencies
% peak and total energy

function frequencyParam=frequencyEstimation(f,energy,varargin)

%% SETUP
p = inputParser;
addParameter(p,'plotFig',0,@isnumeric); % if 1, plot figure
addParameter(p,'newFig',0,@isnumeric); % if 1, plot new figure
addParameter(p,'units','m\cdots^-^2/Hz',@ischar); % units
addParameter(p,'titleFig',[],@ischar); % title fig
addParameter(p,'isIMF',0,@isnumeric); % 1 sum the spectrums if IMF
addParameter(p,'ranges',[],@isnumeric); % range of frequency to calculate an additional totalAmplitude ([linInf 1, limSup 1; liInf 2, limSup 2]
parse(p,varargin{:});
plotFig=p.Results.plotFig;
newFig=p.Results.newFig;
units=p.Results.units;
titleFig=p.Results.titleFig;
isIMF=p.Results.isIMF;

ranges=p.Results.ranges;
% dimensions
energies=transposeColmunIfNot(energy);
f=transposeColmunIfNot(f);

frequencyParam.sep.amplitude=energies;


%% Separated axes
for i=1:size(energies,2)
    energy=energies(:,i);
    mainFrequency=f(energy==max(energy));
    sumEnergy=cumsum(energy,'omitnan'); % cumul de chaque scale au cours du temps    for j=1:size(CUM,2) % pour chaque temps
    medianFrequency=f(find(sumEnergy>sumEnergy(end)/2,1)); % trouver le cumul > moitié de la puissance max
    totalEnergy=trapz(f,energy);
    meanFrequency=trapz(f,energy.*f)/totalEnergy;
    
    frequencyParam.sep.main(i)=mainFrequency;
    frequencyParam.sep.median(i)=medianFrequency;
    frequencyParam.sep.mean(i)=meanFrequency;
    frequencyParam.sep.maxAmplitude(i)=max(energy);
    frequencyParam.sep.totalAmplitude(i)=totalEnergy;
    
    % range
    if ~isempty(ranges)
        indRange=findIndRange(f,ranges);
        frequencyParam.sep.range.ranges=ranges;
        for j=1:numel(indRange)
            energy2range=energy(indRange{j});
            frequencyParam.sep.range.totalAmplitude(j,i)=trapz(f(indRange{j}),energy2range);
        end
    end
    
end
frequencyParam.sep.ratioTotalAmplitude=100*frequencyParam.sep.totalAmplitude/sum(frequencyParam.sep.totalAmplitude);
%% Norm of the signal
if size(energies,2)>1
    if isIMF==0
        energy=signalNorm(energies);
    else
        energy=sum(energies,2);
    end
    mainFrequency=f(energy==max(energy));
    sumEnergy=cumsum(energy,'omitnan'); % cumul de chaque scale au cours du temps    for j=1:size(CUM,2) % pour chaque temps
    medianFrequency=f(find(sumEnergy>sumEnergy(end)/2,1)); % trouver le cumul > moitié de la puissance max
    totalEnergy=trapz(f,energy);
    meanFrequency=trapz(f,energy.*f)/totalEnergy;
    
    frequencyParam.norm.amplitude=energy;
    frequencyParam.norm.main=mainFrequency;
    frequencyParam.norm.median=medianFrequency;
    frequencyParam.norm.mean=meanFrequency;
    frequencyParam.norm.maxAmplitude=max(energy);
    frequencyParam.norm.totalAmplitude=totalEnergy;
    
    % range
    if ~isempty(ranges)
        indRange=findIndRange(f,ranges);
        frequencyParam.norm.range.ranges=ranges;
        
        for j=1:numel(indRange)
            energy2range=energy(indRange{j});
            frequencyParam.norm.range.totalAmplitude(j)=trapz(f(indRange{j}),energy2range);
        end
    end
    
end
frequencyParam.f=f;

%% PLOT
if plotFig==1
    if newFig==1
        figure
    end
    
    % find close points for mean frequency
    [~,meanFrequency4plot]=min(abs(meanFrequency-f));
    meanFrequency4plot=[meanFrequency4plot-1 meanFrequency4plot meanFrequency4plot+1];
    
    if size(energies,2)>1 & size(energies,2)<4
        plot(f,energy,'k','HandleVisibility','on'); hold on
        plot(f,energies(:,1),'k-.','HandleVisibility','on');
        plot(f,energies(:,2),'k-.','HandleVisibility','off');
        if size(energies,2)>2
            plot(f,energies(:,3),'k-.','HandleVisibility','off');
        end
    elseif size(energies,2)==1
        plot(f,energy,'k','HandleVisibility','off'); hold on
    else
        plot(f,energy(:,1),'k','HandleVisibility','on'); hold on
        plot(f,energies(:,1),'-.k','HandleVisibility','on');
        plot(f,energies(:,2:end),'-.k','HandleVisibility','off');
    end
    xlim([min(f) max(f)])
    xlabel('Frequency (Hz)');
    ylabel(['Amplitude (' units ')'])
    scatter(mainFrequency,energy(f==mainFrequency),'k+','HandleVisibility','off')
    text(mainFrequency+0.025*max(f),energy(f==mainFrequency),['Peak amplitude = ' sprintf('%0.1f',energy(f==mainFrequency)) ' ' units ' at ' sprintf('%0.1f',mainFrequency) ' Hz'])
    vline(medianFrequency,'displayLegend',1,'vlimits',[0 energy(f==medianFrequency)]);
    vline(meanFrequency,'displayLegend',1,'linetype',':k','vlimits',[0 max(energy(meanFrequency4plot))])
    if size(energies,2)==1
        legend({['Median frequency = ' sprintf('%0.1f',medianFrequency) ' Hz'],['Mean frequency = ' sprintf('%0.1f',meanFrequency) ' Hz']},'Box','off')
    else
        if isIMF==1
            legend({'sum of sub signals','sub signals',['Median frequency = ' sprintf('%0.2f',medianFrequency) ' Hz'],['Mean frequency = ' sprintf('%0.2f',meanFrequency) ' Hz']},'Box','off')
        else
            legend({'norm of signals','signals from different axes',['Median frequency = ' sprintf('%0.2f',medianFrequency) ' Hz'],['Mean frequency = ' sprintf('%0.2f',meanFrequency) ' Hz']},'Box','off')
        end
    end
    legPos=get(legend,'position');
    set(legend,'position',[legPos(1) legPos(2)-0.1 legPos(3:4)])
    box off
    yLimits = get(gca, 'YLim');
    ylim(1.1*yLimits);
    title(titleFig)
    
end
end
