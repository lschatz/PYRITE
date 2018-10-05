function [trispatialfreq, quadspatialfreq, xp,yp,cp]= fourierbasis(npix,nmodes,triresidual,quadresidual)

%The number of modes returned in fbs is (2n+1)^2-1.
[fourier,xp,yp,cp]=makefbs(npix,nmodes);
%%
sz=size(triresidual);
sz2=size(fourier);

%TRI

    for j=1:sz2(3)
    freq=triresidual.*fourier(:,:,j);
    freq=sum(sum(freq));
    trispatialfreq(:,j)=freq;
    end


%QUAD
sz=size(quadresidual);
sz2=size(fourier);


    for j=1:sz2(3)
    freq=quadresidual.*fourier(:,:,j);
    freq=sum(sum(freq));
    quadspatialfreq(:,j)=freq;
    end


end