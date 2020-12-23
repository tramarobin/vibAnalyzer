function [newFs,preImpactPoints,postImpactPoints,newPreImpactPoints,newPostImpactPoints,acc]=defineTime(acc,Fs,newFs,preImpact,postImpact,reflection)


if isempty(newFs)
    newFs=Fs;
end

if ~isempty(preImpact)
    preImpactPoints=preImpact*Fs;
    newPreImpactPoints=preImpactPoints*newFs/Fs;
else
    preImpactPoints=1;
    newPreImpactPoints=1;
end

if isempty(postImpact)
    if preImpactPoints==1
        postImpactPoints=size(acc,1)-preImpactPoints+1;
    else
        postImpactPoints=size(acc,1)-preImpactPoints;
    end
    newPostImpactPoints=postImpactPoints*newFs/Fs;
else
    postImpactPoints=postImpact*Fs;
    newPostImpactPoints=postImpact*newFs;
end

if reflection==1
    if preImpactPoints>1
    acc=acc(preImpactPoints+1:end,:);
    end
    preImpactPoints=1;
    newPreImpactPoints=1;
end
% newPreImpactPoints=ceil(newPreImpactPoints);
% newPostImpactPoints=ceil(newPostImpactPoints);
end