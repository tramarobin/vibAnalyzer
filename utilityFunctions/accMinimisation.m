function err=accMinimisation(measuredAcc,Fs,optiP)


t=transpose(0:1/Fs:(numel(measuredAcc)-1)/Fs);
modelAcc=(optiP(1).*exp(-optiP(2).*t)).*sin((2.*pi.*optiP(3).*t)+optiP(4));

err=sum((measuredAcc-modelAcc).^2);


end