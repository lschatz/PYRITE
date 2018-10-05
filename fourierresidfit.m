function [trival, quadval]=fourierresidfit(quadresidual, triresidual, fbs,kval)

%%
sz=size(triresidual);
sz2=size(fbs);
for i=1:sz(3)
    for j=1:sz2(3)
    trivals(j,i)=sum(sum(triresidual(:,:,i).*fbs(:,:,j)));
    quadvals(j,i)=sum(sum(quadresidual(:,:,i).*fbs(:,:,j)));
    end
end

% %% Normalize by the kolmogorov phase screen values
% 
sz=size(trivals);
for i=1:sz(2)
    trival(:,i)=trivals(:,i)./kval(:,i);
    quadval(:,i)=quadvals(:,i)./kval(:,i);

end


%% Normalize
% sz=size(kval);
% for i=1:sz(2)
%  Rwavefront(:,i)=Rwavefront(:,i)./kval(:,i);
% end

% norm=rms(kval,2);
% 
% trival=trival./norm;
% quadval=quadval./norm;

end

