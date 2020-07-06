function Norm=signalNorm(Signal)

Norm=zeros(size(Signal,1),1);
for i=1:size(Signal,2)
    
    
    Norm=Norm+Signal(:,i).^2;
    
end

Norm=sqrt(Norm);





end