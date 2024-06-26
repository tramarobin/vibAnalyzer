% estimates the characteristics of the vibrations by fitting a model and return these characteristics as output
% model = amplitude*exp(-damping.*t).*sin(2*pi*frequency.*t+phase);

function modelParam=accEstimation2f(acc,varargin)
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
    [optiParam, err]=fmincon(@(optiP)accMinimisation2f(acc(:,i),Fs,optiP),...
        [0 10 10 0 0 30 30 0],[],[],[],[],[-inf 0 0 0 -inf 0 0 0],[],[],options);
    
    t=transpose(1/Fs:1/Fs:size(acc,1)/Fs);
    subModelAcc(:,1)=(optiParam(1).*exp(-optiParam(2).*t)).*sin((2.*pi.*optiParam(3).*t)+optiParam(4));
    subModelAcc(:,2)=(optiParam(5).*exp(-optiParam(6).*t)).*sin((2.*pi.*optiParam(7).*t)+optiParam(8));
modelAcc=sum(subModelAcc,2);
    
modelParam.measuredAcc(:,i)=acc(:,i);
modelParam.modelAcc(:,i)=modelAcc;
modelParam.amplitude(i,:)=optiParam([1 5]);
modelParam.damping(i,:)=optiParam([2 6]);
modelParam.dampingRatio(i,1)=optiParam(2)/(optiParam(3)*sqrt(4*pi.^pi+optiParam(2).^2/optiParam(3).^2));
modelParam.dampingRatio(i,2)=optiParam(6)/(optiParam(7)*sqrt(4*pi.^pi+optiParam(6).^2/optiParam(7).^2));
modelParam.frequency(i,:)=optiParam([3 7]);
modelParam.phase(i,:)=optiParam([4 8]);
modelParam.error(i)=err;
modelParam.r(i)=corr(modelAcc,acc(:,i));
end


%% PLOT

if plotFig==1
    if newFig==1
        figure
    end
    
    time=1/Fs:1/Fs:size(acc,1)/Fs;
    
    for i=1:size(acc,2)
        subplot(size(acc,2),1,i)
        plot(time,[acc(:,i),modelParam.modelAcc(:,i)],'linewidth',1.2); hold on
        plot(time, subModelAcc,'-.k')
        legend({'Original signal','Modeled signal','subModel sinuses'},'box','off')
        box off
        xlabel('Time (s)')
        ylabel('Acceleration (m\cdots^-^2)')
        
    end
    
end



end
