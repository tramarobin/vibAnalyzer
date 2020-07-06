function modelParam=accEstimation(acc,varargin)
p = inputParser;
addParameter(p,'Fs',1000,@isnumeric); % samplefrequency
addParameter(p,'plotFig',0,@isnumeric); % if 1, plot figure (for one axis, or the norm)
addParameter(p,'newFig',0,@isnumeric); % if 1, plot new figure
parse(p,varargin{:});
Fs=p.Results.Fs;
plotFig=p.Results.plotFig;
newFig=p.Results.newFig;


acc=transposeColmunIfNot(acc);
options = optimoptions('fmincon','Display','off');

for i=1:size(acc,2)
    [optiParam, err]=fmincon(@(optiP)accMinimisation(acc(:,i),Fs,optiP),...
        [0 10 15 0],[],[],[],[],[],[],[],options);
    
    t=transpose(1/Fs:1/Fs:size(acc,1)/Fs);
    modelAcc=optiParam(1)*exp(-optiParam(2).*t).*sin(2*pi*optiParam(3).*t+optiParam(4));
    
    modelParam.measuredAcc(:,i)=acc(:,i);
    modelParam.modelAcc(:,i)=modelAcc;
    modelParam.amplitude(i)=optiParam(1);
    modelParam.damping(i)=optiParam(2);
    modelParam.frequency(i)=optiParam(3);
    modelParam.phase(i)=optiParam(4);
    modelParam.error(i)=err;
    
end


%% PLOT

if plotFig==1
    if newFig==1
        figure
    end
    
    time=1/Fs:1/Fs:size(acc,1)/Fs;
    
    for i=1:size(acc,2)
        subplot(size(acc,2),1,i)
        plot(time,[acc(:,i),modelParam.modelAcc(:,i)])
        legend({'Original signal','Modeled signal'},'box','off')
        box off
    end
    
end



end
