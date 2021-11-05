clear all; close all; clc;
%datMat=load('tunneling.txt');
datMat=load('partial_tunneling.txt');
%addpath('/Users/mahadevan-group/Documents/MATLAB/matlab2tikz/src');
addpath('/home/kalyani/Documents/MATLAB/matlab2tikz-matlab2tikz-56c94a1/src');

len = size(datMat);
x = datMat(:,1);
%x = spread_phase(1,:)-5*pi;
kind = 2;
nr = 3;
nc = 1;

maxEvol = [];
xilLoc = [];
tabX = [];
for ind = 1:(len(2)/3-1)
    cp = datMat(:,kind);
    rhoa = datMat(:,kind+1);
    rhos = datMat(:,kind+2);
    
    ixLoc = find(rhoa == max(rhoa));
    %ixLoc = find(abs(rhos - 0.5) <= 0.05);
    xiL(ind) = x(ixLoc(1));
    
    subplot(nr,nc,1);
    %semilogy(kind,xiL,'ko'); hold on;
    plot(ind/(len(2)/3-1),xiL(ind),'ko'); hold on;
    tabX(ind, :) = [ind/(len(2)/3-1) xiL(ind)];
    pbaspect([1 1 1]);
    xlabel('t');
    ylabel('x');
    
    %     ixLoc = find(cp == max(cp));
    %     xiL = x(ixLoc);
    %
    %     subplot(nr,nc,5);
    %     loglog(kind,-xiL,'ko'); hold on;
    %     xilLoc(ind,:) = [kind -xiL];
    %     pbaspect([1 1 1]);
    %     xlabel('t');
    %     ylabel('x');
    %
    %     subplot(nr,nc,6);
    %     plot(kind,max(cp),'ko'); hold on;
    %     plot(kind,max(rhob),'ro');
    %     maxEvol(ind,:) = [kind max(cp) max(rhob)];
    %     pbaspect([1 1 1]);
    %     xlabel('t');
    %     ylabel('Amplitude');
    if rem(ind,5) == 0
        
        %         subplot(nr,nc,1);
        %         %plot(x,rhob./max(rhob)); hold on;
        %         plot(x,rhoa); hold on;
        %         xlabel('x');
        %         ylabel('c_p(x)');
        %         pbaspect([2 1 1]);
        
        subplot(nr,nc,2);
        %     plot(x,cp./max(cp)); hold on;
        plot(x,cp); hold on;
        xlabel('x');
        ylabel('\rho_a(x)');
        pbaspect([2 1 1]);
        
        subplot(nr,nc,3);
        %     plot(x,cp./max(cp)); hold on;
        plot(x,rhos); hold on;
        xlabel('x');
        ylabel('\rho_s(x)');
        pbaspect([2 1 1]);
        
        
        %         subplot(nr,nc,4);
        %         plot(-x./xiL,cp./max(cp)); hold on;
        % %         xiLFull = [xiLFull -x'./xiL];
        %         cpScFull = [cpScFull -x'./xiL cp'./max(cp)];
        %         xlim([-2.5 2.5]);
        %         pbaspect([2 1 1]);
        %         xlabel('\tilde{x}');
        %         ylabel('\varrho(x)/max(\varrho)');
        %
        %
        %         %     subplot(nr,nc,5);
        %         ixLoc = find(rhob >= 1e-2);
        %         xiL = min(x(ixLoc));
        %         %     plot(kind,xiL,'ko'); hold on;
        %         %     pbaspect([1 1 1]);
        %         subplot(nr,nc,3);
        %         plot(-x./xiL,rhob./max(rhob)); hold on;
        % %         xiLRhoFull = [xiLRhoFull -x'./xiL];
        %         rhoScFull = [rhoScFull -x'./xiL rhob'./max(rhob)];
        %
        %         xlim([-2.5 2.5]);
        %         pbaspect([2 1 1]);
        %         xlabel('\tilde{x}');
        %         ylabel('c_p(x)/max(c_p)');
        
    end
    kind = kind + 3;
end
writematrix(tabX,'xLocPartialTunnel.txt','Delimiter','\t');
% writematrix(cpFull(1:2:end,:),'cpFull.txt','Delimiter','\t');
% writematrix(rhoFull(1:2:end,:),'rhoFull.txt','Delimiter','\t');
% writematrix(cpScFull(1:2:end,:),'cpScFull.txt','Delimiter','\t');
% writematrix(rhoScFull(1:2:end,:),'rhoScFull.txt','Delimiter','\t');
% writematrix(maxEvol,'maxEvol.txt','Delimiter','\t');
% writematrix(xilLoc,'xilLoc.txt','Delimiter','\t');
%  cleanfigure; matlab2tikz('1D_growth.tex');