clear all; clc; close all;
addpath('/Users/mahadevan-group/Documents/MATLAB/matlab2tikz/src');
n = 1;
m = 3;
imagefiles = dir('*.png');
nfiles = length(imagefiles);    % Number of files found
for ind = 1:5
    fname = imagefiles(ind).name;
    cimg = imread(fname);
    cimg = cimg(:,:,1);
    %plot(cimg(1:xcen,ycen),'o-'); hold on;
    rProf = [];
    k = 1;
    angLst = 0:5:360;
    %subplot(n,m,1);
    %imshow(cimg);
    for ang = angLst
        rot = imrotate(cimg,-ang);
        
        cen = floor(size(rot)/2);
        xcen = cen(1);
        ycen = cen(2);
        
        rval = 255-rot(xcen,ycen:end);
        rval(find(rval==255)) = 0.0;
        %subplot(n,m,2);
        %plot(rval,'o-'); hold on;
        %axis square;
        
        rProf(k) = trapz(rval);
        k = k+1;
    end
    
    %subplot(n,m,3);
    finLst = rProf./trapz(angLst*2*pi/360,rProf);
    finAng =  angLst*2*pi/360 - 3.491;
    %plot(finAng,finLst,'-'); hold on;
    plot(finAng,smoothdata(finLst,'gaussian'),'-'); hold on;
    axis square;
end
matlab2tikz('radialDensityUpdated.tex');
%subplot(n,m,3);
%imshow(rot);