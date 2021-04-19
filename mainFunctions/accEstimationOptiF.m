% estimates the characteristics of the vibrations by fitting a model and return these characteristics as output
% model = amplitude*exp(-damping.*t).*sin(2*pi*frequency.*t+phase);

function modelParam=accEstimationOptiF(acc,varargin)
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
    LB=[-inf 0 0 0];
    for k=1:10
        [optiParam{k}, err(k)]=fmincon(@(optiP)accMinimisationOptiF(acc(:,i),Fs,k,optiP),...
            [zeros(1,4*k)],[],[],[],[],repmat(LB,1,k),[],[],options);
    end
    errP=err.*(1:10);
    [errP,nf]=min(errP);
    optiParam=optiParam{nf};
    err=err(nf);
    
    t=transpose(1/Fs:1/Fs:size(acc,1)/Fs);
    for j=1:numel(optiParam)/4
        subModelAcc(:,j)=(optiParam(4*(j-1)+1).*exp(-optiParam(4*(j-1)+2).*t)).*sin((2.*pi.*optiParam(4*(j-1)+3).*t)+optiParam(4*(j-1)+4));
    end
    
    modelAcc=sum(subModelAcc,2);
    
    modelParam.measuredAcc(:,i)=acc(:,i);
    modelParam.modelAcc(:,i)=modelAcc;
    modelParam.subModelAcc{i}=subModelAcc;
    
    for j=1:numel(optiParam)/4
        modelParam.amplitude(i,j)=optiParam(4*(j-1)+1);
        modelParam.damping(i,j)=optiParam(4*(j-1)+2);
        modelParam.frequency(i,j)=optiParam(4*(j-1)+3);
        modelParam.dampingRatio(i,j)=modelParam.damping(i,j)/(modelParam.frequency(i,j)*sqrt(4*pi.^pi+modelParam.damping(i,j).^2/modelParam.frequency(i,j).^2));
        modelParam.phase(i,j)=optiParam(4*(j-1)+4);
    end
    
    modelParam.error(i)=err;
    modelParam.r(i)=corr(modelAcc,acc);
    
    clear optiParam err subModelAcc
    
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
        plot(time, modelParam.subModelAcc{i},'-.k')
        legend({'Original signal','Modeled signal','subModel sinuses'},'box','off')
        box off
        xlabel('Time (s)')
        ylabel('Acceleration (m\cdots^-^2)')
        
    end
    
end



end
