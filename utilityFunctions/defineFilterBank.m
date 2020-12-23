function [fb,f,infLim,supLim,infFreq,supFreq,intFreq,indIntInfLim,indIntSupLim]=defineFilterBank(acc,Fs,interpFreq,fb,infFreq,supFreq,reflection)

if reflection==0
    if size(acc,1)>round(2.048*Fs)
        signalLength=size(acc,1);
    else
        signalLength=round(2.048*Fs);
    end
else
    if size(acc,1)>round(2.048*Fs)
        signalLength=2*size(acc,1);
    else
        signalLength=2*round(2.048*Fs);
    end
end


%% FB
if isempty(fb)
    fb=cwtfilterbank(...
        'SamplingFrequency',Fs,...
        'SignalLength',signalLength,...
        'VoicesPerOctave',48);
else
    if fb.SignalLength~=signalLength
        if ~isempty(fb.WaveletParameters)
            fb2=cwtfilterbank(...
                'SamplingFrequency',Fs,...
                'SignalLength',signalLength,...
                'VoicesPerOctave',fb.VoicesPerOctave,...
                'Wavelet',fb.Wavelet,...
                'WaveletParameters',fb.WaveletParameters);
        else
            fb2=cwtfilterbank(...
                'SamplingFrequency',Fs,...
                'SignalLength',signalLength,...
                'VoicesPerOctave',fb.VoicesPerOctave,...
                'Wavelet',fb.Wavelet);
        end
        fb=fb2; clear fb2

    end
end

%% FREQ
f=centerFrequencies(fb);
[minfreq maxfreq]=cwtfreqbounds(fb.SignalLength,fb.SamplingFrequency,'Wavelet',fb.Wavelet);

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