function [out]=valnorm(input,kval);

sz=size(input);
for i=1:sz(2)
    out(:,i)=input(:,i)./kval(:,i);
end
end
