function [coefficients]=waveletAnalysis(acc,fb,varargin)
p = inputParser;
addParameter(p,'isIMF',0,@isnumeric); % sum the map if IMF
addParameter(p,'reflection',0,@isnumeric); % 1 use reflection at the start of the signal and add 0 padding of 2048 points centered on heel strike, it improve mode separation and allow to investigate lower frequencies, /!\ the signal analyzed is not the one you measured anymore. Enders et al. 2012; http://dx.doi.org/10.1016/j.jbiomech.2012.08.027
parse(p,varargin{:});
isIMF=p.Results.isIMF;
reflection=p.Results.reflection;

%% each axis
if reflection==0
    for i=1:size(acc,2)
        cfs=wt(fb,acc(:,i));
        axisCoefficients{i}=abs(cfs);
    end
    
else
    Fs=fb.SamplingFrequency;
    originalSize=size(acc,1);
    paddingSize=2.048*Fs-size(acc,1);
    acc=[zeros(paddingSize,size(acc,2)); acc(end:-1:1,:); acc; zeros(paddingSize,size(acc,2))]; % reflect the signal
    for i=1:size(acc,2)
        cfs=wt(fb,acc(:,i));
        axisCoefficients{i}=abs(cfs(:,2.048*Fs+1:2.048*Fs+originalSize));
    end
    
end

coefficients.axisCoefficients=axisCoefficients;

%% norm
if size(acc,2)>1
    Coeffs=0;
    if isIMF==0
        for i=1:size(acc,2)
            Coeffs=Coeffs+axisCoefficients{i}.^2;
        end
        coefficients.normCoefficients=sqrt(Coeffs);
        
    else
        for i=1:size(acc,2)
            Coeffs=Coeffs+axisCoefficients{i};
        end
    end
end


end