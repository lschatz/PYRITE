function [imageout]=addnoise(imagein, rdns)

%imagein=PupilOne;
root=imagein.^0.5;
imageout=normrnd(imagein, root);
imageout=normrnd(imageout, rdns);

end
