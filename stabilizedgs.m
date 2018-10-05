function [out]=stabilizedgs(in)

%% reshape array so that each fourier mode is a row vector
sz1=size(in);
for i=1:sz1(3)
A=reshape(in(:,:,i),1,[]);
v(i,:,:)=A;
end
v=squeeze(v);



%% Stabilized Gram Schmidt
% v= input matrix
% u= output matrix


sz=size(v);
u=zeros(sz(1),sz(2));

   u(1,:)=v(1,:)/norm(v(1,:)); % normalize
    
   for i=2:sz(1)
     u(i,:)=v(i,:);
     for j=1:i-1
     u(i,:)=u(i,:)-dot(u(i,:),u(j,:)).*u(j,:); % subtract off the projection
     end
     u(i,:)=u(i,:)/norm(u(i,:)); % normalize
   end

dot(u(1,:),u(2,:)) %check that it is close to zero
%% Reshape back
for i=1:sz1(3)
    B=reshape(u(i,:,:),[sz1(1),sz1(2)]);
    out(:,:,i)=B;
end
