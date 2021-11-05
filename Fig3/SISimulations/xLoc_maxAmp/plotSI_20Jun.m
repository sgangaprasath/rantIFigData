clear all; close all; clc;
addpath('/Users/mahadevan-group/Documents/MATLAB/matlab2tikz/src');
%addpath('/home/kalyani/Documents/MATLAB/matlab2tikz-matlab2tikz-56c94a1/src');
%datMat=load('tunneling.txt');

dinfo(1) = dir('partialtunneling_19Jun.txt');
dinfo(2) = dir('partialtunnelingII_19Jun.txt');
dinfo(3) = dir('tunneling_19Jun.txt');
%dinfo = dir('p*.txt');

nr = 1;
nc = 3;

normFn = fittype('exp((x-mu).^2./(2*sig.^2))');

for K = 1 : length(dinfo)
    thisfilename = dinfo(K).name;
    datMat = load(thisfilename);
    kind = 1;
    len = size(datMat);
    x = datMat(:,1);
    %x = spread_phase(1,:)-5*pi;
    maxEvol = [];
    xilLoc = [];
    xLEvol = [];
    for ind = 1:(len(2)/3-1)
        %         if ind == 1
        %             subplot(nr,nc,4);
        %             plot(x, exp(-(x-0.4).^2/0.05),'-'); hold on;
        %         end
        
        rhos = datMat(:,kind);
        cp = datMat(:,kind+1);
        rhoa = datMat(:,kind+2);
        maxEvol(ind,:) = [0.05*kind/3 max(rhoa) max(cp)];
        
        ixLoc = find(rhoa == max(rhoa));
        %widLoc = find((rhoa >= 0.5*max(rhoa)-0.1) & (rhoa <= 0.5*max(rhoa)+0.5));
        xiL(ind) = x(ixLoc(1));
        %widL(ind) = abs(x(max(widLoc)) - x(min(widLoc)));
        widL(ind) = var(x,abs(rhoa));
        %f0 = fit(x,rhoa,normFn);
        %widL(ind) = f0.sig;
        xLEvol(ind,:) = [0.05*kind/3 xiL(ind) sqrt(widL(ind))];
        
%                 if rem(ind,2) == 0
%                     subplot(nr,nc,4);
%                     plot(x+K*0.25,rhoa,'Linewidth',1.2); hold on;
%                     xlabel('x');
%                     ylabel('\rho_a(x)');
%                     pbaspect([2 1 1]);
%                     xlim([0.25 2]);
%                     ylim([0 1]);
%         
%                     subplot(nr,nc,5);
%                     plot(x+K*0.25,rhos,'Linewidth',1.2); hold on;
%                     xlabel('x');
%                     ylabel('\rho_s(x)');
%                     pbaspect([2 1 1]);
%                     xlim([0.25 2]);
%                     ylim([0 1]);
%         
%                     subplot(nr,nc,6);
%                     plot(x+K*0.25,cp); hold on;
%                     xlabel('x');
%                     ylabel('c(x)');
%                     pbaspect([2 1 1]);
%                     xlim([0.25 2]);
%                     ylim([0 1]);
%                 end
        
        kind = kind + 3;
    end
    
    subplot(nr,nc,1);
    %semilogy(kind,xiL,'ko'); hold on;
    plot(xLEvol(:,1),xLEvol(:,2)); hold on;
    xlabel('t');
    ylabel('xLoc');
    xlim([0 7]);
    pbaspect([2 1 1]);
    
    subplot(nr,nc,2);
    plot(xLEvol(:,1),xLEvol(:,3)); hold on;
    xlabel('t');
    xlim([0 7]);
    pbaspect([2 1 1]);
    ylabel('width');
    
    
    
    subplot(nr,nc,3);
    plot(maxEvol(:,1),maxEvol(:,2)); hold on;
    plot(maxEvol(:,1),maxEvol(:,3));
    xlim([0 7]);
    pbaspect([2 1 1]);
    xlabel('t');
    %ylabel('max(c)');
    
    %    sName = sprintf('maxEvol%d.dat',K);
    %    writematrix(maxEvol,sName,'Delimiter','\t');
    
    %    sName = sprintf('xLEvol%d.dat',K);
    %    writematrix(xLEvol,sName,'Delimiter','\t');
end
matlab2tikz('xlocPlot.tex')