% compute the intrinsec mode functions (IMF) of the signals with partly emsemble empirical
% mode decomposition
% if there is serveral signals, the algorithm calcul the IMFs for each axis independantly
% compute the FFT and the CWT of each IMFs of each axis

%% INPUT
% acceleration signal in m/s² in column (time,axes)

%% OPTIONS
% see addParameter below

%% OUPUTS
% one structure containing the IMFs per axis (.IMF{axis})
% the temporal analysis (.TEMPORAL{axis}), the damping is computed with the
% decaying envelope
% the FFT transforms (.FF{axis}) and the CWT (.CWT{axis})

function peemdParam=peemdAnalysis(acc,varargin)

p = inputParser;
addParameter(p,'infFreq',6,@isnumeric); % Freq min analyzed
addParameter(p,'supFreq',[],@isnumeric); % Freq max analyzed
addParameter(p,'Fs',1000,@isnumeric); % sample frequency
addParameter(p,'padding',2048,@isnumeric); % number of point for padding the FFT analysis
addParameter(p,'plotFig',0,@isnumeric); % if 1, plot figure
addParameter(p,'newFig',0,@isnumeric); % if 1, plot new figure
addParameter(p,'postImpact',[],@isnumeric); % total time analyzed (default = end of signal)
addParameter(p,'interpFreq',1,@isnumeric); % change interpolation frequency to reduce map size
addParameter(p,'newFs',[],@isnumeric); % change sample frequency to reduce map size (must be < Fs)
addParameter(p,'reflection',0,@isnumeric); % 1 use reflection at the start of the signal and add 0 padding of 2048 points centered on heel strike, it improve mode separation and allow to investigate lower frequencies, /!\ the signal analyzed is not the one you measured anymore. Enders et al. 2012; http://dx.doi.org/10.1016/j.jbiomech.2012.08.027
addParameter(p,'noiseAmplitude',0.2,@isnumeric); % 1 use reflection at the start of the signal and add 0 padding of 2048 points centered on heel strike, it improve mode separation and allow to investigate lower frequencies, /!\ the signal analyzed is not the one you measured anymore. Enders et al. 2012; http://dx.doi.org/10.1016/j.jbiomech.2012.08.027
addParameter(p,'nRealisations',100,@isnumeric); % 1 use reflection at the start of the signal and add 0 padding of 2048 points centered on heel strike, it improve mode separation and allow to investigate lower frequencies, /!\ the signal analyzed is not the one you measured anymore. Enders et al. 2012; http://dx.doi.org/10.1016/j.jbiomech.2012.08.027
addParameter(p,'maxIterations',500,@isnumeric); % 1 use reflection at the start of the signal and add 0 padding of 2048 points centered on heel strike, it improve mode separation and allow to investigate lower frequencies, /!\ the signal analyzed is not the one you measured anymore. Enders et al. 2012; http://dx.doi.org/10.1016/j.jbiomech.2012.08.027

parse(p,varargin{:});
infFreq=p.Results.infFreq;
supFreq=p.Results.supFreq;
Fs=p.Results.Fs;
padding=p.Results.padding;
plotFig=p.Results.plotFig;
newFig=p.Results.newFig;
postImpact=p.Results.postImpact;
interpFreq=p.Results.interpFreq;
newFs=p.Results.newFs;
reflection=p.Results.reflection;
noiseAmplitude=p.Results.noiseAmplitude;
nRealisations=p.Results.nRealisations;
maxIterations=p.Results.maxIterations;

acc=transposeColmunIfNot(acc);

%% IMF and FFT
for i=1:size(acc,2)
    peemdParam.IMF{i}=peemd(acc(:,i),noiseAmplitude,nRealisations,maxIterations)';
    fftParam=fftAnalysis(peemdParam.IMF{i},'padding',padding,'Fs',Fs,'infFreq',1);
    if ~isempty(infFreq) || ~isempty(supFreq)
        if ~isempty(infFreq)
            peemdParam.IMF{i}=peemdParam.IMF{i}(:,fftParam.FT.sep.main>infFreq);
            fftParam.FT.sep.main=fftParam.FT.sep.main(fftParam.FT.sep.main>infFreq);
            fftParam.FT.sep.ratioTotalAmplitude=fftParam.FT.sep.ratioTotalAmplitude(fftParam.FT.sep.main>infFreq);
        end
        if ~isempty(supFreq)
            peemdParam.IMF{i}=peemdParam.IMF{i}(:,fftParam.FT.sep.main<supFreq);
            fftParam.FT.sep.main=fftParam.FT.sep.main(fftParam.FT.sep.main<supFreq);
            fftParam.FT.sep.ratioTotalAmplitude=fftParam.FT.sep.ratioTotalAmplitude(fftParam.FT.sep.main<supFreq);
            
        end
    end
    
    peemdParam.IMF{i}=peemdParam.IMF{i}(:,fftParam.FT.sep.ratioTotalAmplitude>5); % only if IMF represent 5% of the signal energy
    %% FFT and CWT
    if ~isempty(peemdParam.IMF{i})
        peemdParam.FFT{i}=fftAnalysis(peemdParam.IMF{i},'padding',padding,'Fs',Fs,'infFreq',infFreq,'supFreq',supFreq,'isIMF',1);
        peemdParam.CWT{i}=cwtAnalysis(peemdParam.IMF{i},'Fs',Fs,'infFreq',infFreq,'supFreq',supFreq,'isIMF',1,'interpFreq',interpFreq,'newFs',newFs,'reflection',1,'delay2peak',0.1);
    else
        warning('No IMF found in the frequency range');
    end
    
    peemdParam.TEMPORAL{i}=temporalAnalysis(peemdParam.IMF{i},'isIMF',1);
    
end


%% PLOT and damp
if plotFig==1
    if newFig==1
        figure('units','normalized','outerposition',[0 0 1 1],'visible','on')
    end
    time=1/Fs:1/Fs:size(acc,1)/Fs;
    for i=1:numel(peemdParam.FFT)
        subplot(2,numel(peemdParam.FFT),i)
        plot(time,acc(:,i),'k--'); hold on
        plot(time,sum(peemdParam.IMF{i},2),'k');
        plot(time,peemdParam.IMF{i},'k-.');
        legend({'Original signal','sum of IMFs','IMFs'},'box','off')
        title(['IMF reconstruction for axe #' num2str(i)])
        subplot(2,numel(peemdParam.FFT),numel(peemdParam.FFT)+i)
        
        f=peemdParam.FFT{i}.normalizedFT.f;
        energy=peemdParam.FFT{i}.normalizedFT.norm.amplitude;
        energies=peemdParam.FFT{i}.normalizedFT.sep.amplitude;
        mainFrequency=peemdParam.FFT{i}.normalizedFT.norm.main;
        meanFrequency=peemdParam.FFT{i}.normalizedFT.norm.mean;
        medianFrequency=peemdParam.FFT{i}.normalizedFT.norm.median;
        
        % find close points for mean frequency
        [~,meanFrequency4plot]=min(abs(meanFrequency-f));
        meanFrequency4plot=[meanFrequency4plot-1 meanFrequency4plot meanFrequency4plot+1];
        
        if size(energies,2)>1
            plot(f,energy,'k','HandleVisibility','on'); hold on
            plot(f,energies(:,1),'k-.','HandleVisibility','on');
            plot(f,energies(:,2:end),'k-.','HandleVisibility','off');
        else
            plot(f,energy,'k','HandleVisibility','off'); hold on
        end
        xlim([min(f) max(f)])
        xlabel('Frequency (Hz)');
        ylabel('Amplitude (m\cdots^-^2/Hz)')
        scatter(mainFrequency,energy(f==mainFrequency),'k+','HandleVisibility','off')
        text(mainFrequency+0.025*max(f),energy(f==mainFrequency),['Peak amplitude = ' sprintf('%0.2f',energy(f==mainFrequency)) ' m\cdots^-^2/Hz at ' sprintf('%0.2f',mainFrequency) ' Hz'])
        vline(medianFrequency,'displayLegend',1);
        vline(meanFrequency,'displayLegend',1,'linetype',':k')
        if size(energies,2)==1
            legend({['Median frequency = ' sprintf('%0.2f',medianFrequency) ' Hz'],['Mean frequency = ' sprintf('%0.2f',meanFrequency) ' Hz']},'Box','off')
        else
            legend({'sum of IMFs','IMFs',['Median frequency = ' sprintf('%0.2f',medianFrequency) ' Hz'],['Mean frequency = ' sprintf('%0.2f',meanFrequency) ' Hz']},'Box','off')
        end
        legPos=get(legend,'position');
        set(legend,'position',[legPos(1) legPos(2)-0.1 legPos(3:4)])
        box off
        title(['FT of IMFs for axe #' num2str(i)])
    end
end

end