function [val]=fourierfit(basis,input)

%% Fit fourier modes to kolmogorov phase screen

sz=size(input);
sz2=size(basis);

for i=1:sz(3)
    for j=1:sz2(3)
    val(j,i)=sum(sum(input(:,:,i).*basis(:,:,j)));
    end
end





end
