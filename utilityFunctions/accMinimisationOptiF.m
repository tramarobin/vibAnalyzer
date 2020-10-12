function err=accMinimisationOptiF(measuredAcc,Fs,nf,optiP)

t=transpose(0:1/Fs:(numel(measuredAcc)-1)/Fs);

for i=1:nf
subModelAcc(:,i)=(optiP(4*(i-1)+1).*exp(-optiP(4*(i-1)+2).*t)).*sin((2.*pi.*optiP(4*(i-1)+3).*t)+optiP(4*(i-1)+4));
end

modelAcc=sum(subModelAcc,2);
err=sum((measuredAcc-modelAcc).^2);


end