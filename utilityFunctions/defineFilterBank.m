function [fb,f,infLim,supLim,infFreq,supFreq,intFreq,indIntInfLim,indIntSupLim]=defineFilterBank(acc,Fs,interpFreq,fb,infFreq,supFreq,reflection)

%% FB
if reflection==0
    if isempty (fb)
        [minfreq,maxfreq] = cwtfreqbounds(size(acc,1),Fs,...
            'VoicesPerOctave',48);
    else
        [minfreq,maxfreq] = cwtfreqbounds(size(acc,1),Fs,...
            'WaveletParameters',fb.WaveletParameters,...
            'VoicesPerOctave',fb.VoicesPerOctave);
    end
else
    if isempty (fb)
        [minfreq,maxfreq] = cwtfreqbounds(4.096*Fs,Fs,...
            'VoicesPerOctave',48);
    else % double padded at 2048 = 4.096*Fs points
        [minfreq,maxfreq] = cwtfreqbounds(4.096*Fs,Fs,...
            'WaveletParameters',fb.WaveletParameters,...
            'VoicesPerOctave',fb.VoicesPerOctave);
    end
end

if isempty(infFreq)
    if minfreq>1
        infFreq=ceil(minfreq);
    else
        infFreq=minfreq;
    end
end
if isempty(supFreq)
    if maxfreq>1
        supFreq=floor(maxfreq);
    else
        supFreq=maxfreq;
    end
end

if ~isempty(infFreq)
    if infFreq<minfreq
        if minfreq>1
            infFreq=ceil(minfreq);
        else
            infFreq=minfreq;
        end
        warning('foo:bar',['low frequency too low, increase signal length\n new low frequency limit = ' num2str(infFreq) ' Hz'])
    end
end
if ~isempty(supFreq)
    if supFreq>maxfreq
        if maxfreq>1
            supFreq=floor(maxfreq);
        else
            supFreq=maxfreq;
        end
        warning('foo:bar',['high frequency too high, increase sample frequency\n new high frequency = ' num2str(supFreq) ' Hz'])
    end
end

if reflection==0
    if isempty(fb)
        fb=cwtfilterbank(...
            'SamplingFrequency',Fs,...
            'SignalLength',size(acc,1),...
            'FrequencyLimits',[minfreq,maxfreq],...
            'VoicesPerOctave',48);
    elseif fb.SignalLength~=size(acc,1)
        [minfreq,maxfreq] = cwtfreqbounds(size(acc,1),Fs,...
            'WaveletParameters',fb.WaveletParameters,...
            'VoicesPerOctave',fb.VoicesPerOctave);
        fb2=cwtfilterbank(...
            'SamplingFrequency',Fs,...
            'SignalLength',size(acc,1),...
            'WaveletParameters',fb.WaveletParameters,...
            'VoicesPerOctave',fb.VoicesPerOctave,...
            'FrequencyLimits',[minfreq,maxfreq]);
        fb=fb2;
    end
else
    if isempty(fb)
        fb=cwtfilterbank(...
            'SamplingFrequency',Fs,...
            'SignalLength',4.096*Fs,...
            'FrequencyLimits',[minfreq,maxfreq],...
            'VoicesPerOctave',48);
    elseif fb.SignalLength~=4.096*Fs
        [minfreq,maxfreq] = cwtfreqbounds(4.096*Fs,Fs,...
            'WaveletParameters',fb.WaveletParameters,...
            'VoicesPerOctave',fb.VoicesPerOctave);
        fb2=cwtfilterbank(...
            'SamplingFrequency',Fs,...
            'SignalLength',4.096*Fs,...
            'WaveletParameters',fb.WaveletParameters,...
            'VoicesPerOctave',fb.VoicesPerOctave,...
            'FrequencyLimits',[minfreq,maxfreq]);
        fb=fb2;
    end
end


%% FREQ

f=centerFrequencies(fb);

infLim=find(f<=infFreq-interpFreq);
if isempty(infLim)
    infLim=numel(f);
elseif minfreq==infFreq
    infLim=infLim(1);
else
    infLim=infLim(1)+1;
end
if infLim>numel(f)
    infLim=numel(f);
end

supLim=find(f>=supFreq+interpFreq);
if isempty(supLim)
    supLim=1;
elseif maxfreq==supFreq
    supLim=supLim(end);
else
    if supLim>1
        supLim=supLim(end)-1;
    else
        supLim=supLim(end);
    end
end
 
% infFreq=ceil(f(infLim)/interpFreq)*interpFreq;
% supFreq=floor(f(supLim)/interpFreq)*interpFreq;

%% Limits for interpolated maps
intFreq=transpose(infFreq:interpFreq:supFreq);

indIntInfLim=find(intFreq>=supFreq,1,'last');
if isempty(indIntInfLim)
    indIntInfLim=numel(intFreq);
end
indIntSupLim=find(intFreq<=infFreq,1);
if isempty(indIntSupLim)
    indIntSupLim=1;
end


end