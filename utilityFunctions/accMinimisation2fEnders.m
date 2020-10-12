function err=accMinimisation2fEnders(measuredAcc,Fs,optiP)


t=transpose(0:1/Fs:(numel(measuredAcc)-1)/Fs);

f=mean([optiP(3) optiP(4)]);
df=diff([optiP(3) optiP(4)]);

modelAcc=(2*optiP(1).*sin((2*pi*f.*t)).*cos((2*pi*(df/2).*t)).*exp(-optiP(2).*t));
err=sum((measuredAcc-modelAcc).^2);


end