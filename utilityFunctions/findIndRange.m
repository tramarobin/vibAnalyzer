function indRange=findIndRange(f,ranges)

for i=1:size(ranges,1)
    indRange{i}=find(f>ranges(i,1) & f<ranges(i,2));    
end


end