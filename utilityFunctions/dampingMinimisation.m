function Err=dampingMinimisation(Pdec,t,d)

Err=sum((Pdec(1).*exp(-d*t)-Pdec).^2);

end