clear all; close all; clc;
%addpath('/Users/mahadevan-group/Documents/MATLAB/matlab2tikz/src');
addpath('/home/kalyani/Documents/MATLAB/matlab2tikz-matlab2tikz-56c94a1/src');

tunl=load('tunneling.txt');
partunl=load('partial_tunneling.txt');
jammed=load('jammed.txt');

len = size(tunl);
x = tunl(:,1);
kind = 2;
nr = 1;
nc = 3;

maxEvol = [];
xilLoc = [];
tabX = [];
indx = 1:(len(2)/3-1);
xiL = zeros(1,max(indx));
mcp = zeros(1,max(indx));
mra = zeros(1,max(indx));
for ind = indx
    cp = tunl(:,kind);
    rhoa = tunl(:,kind+1);
    rhos = tunl(:,kind+2);
    
    ixLoc = find(rhoa == max(rhoa));
    xiL(ind) = x(ixLoc(1));
    
    subplot(nr,nc,1);
    %semilogy(kind,xiL,'ko'); hold on;
    plot(ind/(len(2)/3-1),xiL(ind),'ko'); hold on;
%    tabX(ind, :) = [ind/(len(2)/3-1) xiL(ind)];
    pbaspect([1 1 1]);
    xlabel('t');
    ylabel('x');
    
    subplot(nr,nc,2);
    mcp(ind) = max(cp);
    plot(ind/(len(2)/3-1),mcp(ind),'ko'); hold on;
    xlabel('t');
    ylabel('max(c_p)');
    pbaspect([1 1 1]);
    
    subplot(nr,nc,3);
    mra(ind) = max(rhoa);
    plot(ind/(len(2)/3-1),mra(ind),'ko'); hold on;
    xlabel('t');
    ylabel('\rho_a(x)');
    pbaspect([1 1 1]);
    
    kind = kind + 3;
end
tabX = [indx'/max(ind), mcp', mra', xiL'];
writematrix(tabX,'Tunnel.txt','Delimiter','\t');

len = size(partunl);
x = partunl(:,1);
kind = 2;
indx = 1:(len(2)/3-1);
xiL = zeros(1,max(indx));
mcp = zeros(1,max(indx));
mra = zeros(1,max(indx));
for ind = indx
    cp = partunl(:,kind);
    rhoa = partunl(:,kind+1);
    rhos = partunl(:,kind+2);
    
    ixLoc = find(rhoa == max(rhoa));
    xiL(ind) = x(ixLoc(1));
    
    subplot(nr,nc,1);
    xiL(ind) = x(ixLoc(1));
    %semilogy(kind,xiL,'ko'); hold on;
    plot(ind/(len(2)/3-1),xiL(ind),'ro'); hold on;
%    tabX(ind, :) = [ind/(len(2)/3-1) xiL(ind)];
    pbaspect([1 1 1]);
    xlabel('t');
    ylabel('x');
    
    
    subplot(nr,nc,2);
    mcp(ind) = max(cp);
    plot(ind/(len(2)/3-1),max(cp),'ro'); hold on;
    xlabel('t');
    ylabel('max(c_p)');
    pbaspect([1 1 1]);
    
    subplot(nr,nc,3);
    mra(ind) = max(rhoa);
    plot(ind/(len(2)/3-1),max(rhoa),'ro'); hold on;
    xlabel('t');
    ylabel('\rho_a(x)');
    pbaspect([1 1 1]);
    
    kind = kind + 3;
end
tabX = [indx'/max(ind), mcp', mra', xiL'];
writematrix(tabX,'PartialTunnel.txt','Delimiter','\t');

len = size(jammed);
x = jammed(:,1);
kind = 2;
for ind = 1:(len(2)/3-1)
    cp = jammed(:,kind);
    rhoa = jammed(:,kind+1);
    rhos = jammed(:,kind+2);
    
    ixLoc = find(rhoa == max(rhoa));
    xiL(ind) = x(ixLoc(1));
    
    subplot(nr,nc,1);
    %semilogy(kind,xiL,'ko'); hold on;
    plot(ind/(len(2)/3-1),xiL(ind),'bo'); hold on;
%    tabX(ind, :) = [ind/(len(2)/3-1) xiL(ind)];
    pbaspect([1 1 1]);
    xlabel('t');
    ylabel('x');
    
    
    subplot(nr,nc,2);
    plot(ind/(len(2)/3-1),max(cp),'bo'); hold on;
    xlabel('t');
    ylabel('max(c_p)');
    pbaspect([1 1 1]);
    
    subplot(nr,nc,3);
    plot(ind/(len(2)/3-1),max(rhoa),'bo'); hold on;
    xlabel('t');
    ylabel('max(\rho_a)');
    pbaspect([1 1 1]);
    
    kind = kind + 3;
end

%writematrix(tabX,'xLocPartialTunnel.txt','Delimiter','\t');