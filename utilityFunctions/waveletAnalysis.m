function [coefficients]=waveletAnalysis(acc,fb,varargin)
p = inputParser;
addParameter(p,'isIMF',0,@isnumeric); % sum the map if IMF
addParameter(p,'reflection',0,@isnumeric); % 1 use reflection at the start of the signal and add 0 padding of 2048 points centered on heel strike, it improve mode separation and allow to investigate lower frequencies, /!\ the signal analyzed is not the one you measured anymore. Enders et al. 2012; http://dx.doi.org/10.1016/j.jbiomech.2012.08.027
parse(p,varargin{:});
isIMF=p.Results.isIMF;
reflection=p.Results.reflection;

%% padding
originalSize=size(acc,1);
if reflection==0
    paddingSize=fb.SignalLength-originalSize;
    acc=[acc; zeros(paddingSize,size(acc,2))];
else
    paddingSize=(fb.SignalLength-2*originalSize)/2;
    acc=[zeros(paddingSize,size(acc,2)); acc(end:-1:1,:); acc; zeros(paddingSize,size(acc,2))]; % reflect the signal
end
%% each axis
if reflection==0
    for i=1:size(acc,2)
        cfs=wt(fb,acc(:,i));
        axisCoefficients{i}=abs(cfs(:,1:originalSize));
    end
    
else
 
    for i=1:size(acc,2)
        cfs=wt(fb,acc(:,i));
        axisCoefficients{i}=abs(cfs(:,fb.SignalLength/2+1:fb.SignalLength/2+1+originalSize));
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