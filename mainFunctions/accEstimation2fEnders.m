% estimates the characteristics of the vibrations by fitting a model and return these characteristics as output
% model = amplitude*exp(-damping.*t).*sin(2*pi*frequency.*t+phase);

function modelParam=accEstimation2fEnders(acc,varargin)
p = inputParser;
addParameter(p,'Fs',1000,@isnumeric); % samplefrequency
addParameter(p,'plotFig',0,@isnumeric); % if 1, plot figure (for one axis, or the norm)
addParameter(p,'newFig',0,@isnumeric); % if 1, plot new figure
addParameter(p,'preImpact',[],@isnumeric); % pre impact time (default = start of the signal)
addParameter(p,'postImpact',[],@isnumeric); % total time analyzed (default = end of signal)
parse(p,varargin{:});
Fs=p.Results.Fs;
plotFig=p.Results.plotFig;
newFig=p.Results.newFig;
preImpact=p.Results.preImpact;
postImpact=p.Results.postImpact;

acc=transposeColmunIfNot(acc);

[~,preImpactPoints,postImpactPoints,~,~,acc]=defineTime(acc,Fs,Fs,preImpact,postImpact,0);
acc=acc(preImpactPoints:preImpactPoints+postImpactPoints-1,:);

options = optimoptions('fmincon','Display','off');
for i=1:size(acc,2)
    
    [optiParam, err]=fmincon(@(optiP)accMinimisation2fEnders(acc(:,i),Fs,optiP),...
        [0 0 0 0],[],[],[],[],[-inf 0 0 0],[],[],options);
    
    f=mean([optiParam(3) optiParam(4)]);
    df=diff([optiParam(3) optiParam(4)]);
    
    t=transpose(1/Fs:1/Fs:size(acc,1)/Fs);
    modelAcc=(2*optiParam(1).*sin((2*pi*f.*t)).*cos((2*pi*(df/2).*t)).*exp(-optiParam(2).*t));
    
    modelParam.measuredAcc(:,i)=acc(:,i);
    modelParam.modelAcc(:,i)=modelAcc;
    modelParam.amplitude(i)=optiParam(1);
    modelParam.damping(i)=optiParam(2);
    modelParam.frequency(i,:)=optiParam([3 4]);
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
        xlabel('Time (s)')
        ylabel('Acceleration (m\cdots^-^2)')
        
    end
    
end



end
