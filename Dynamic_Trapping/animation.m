function [] = animation(T,X,s,p)

figure
photo = imshow(zeros(p.ny,p.nx),'InitialMagnification','fit');
hold on
%rants = line('color','r','Marker','o','LineStyle','None','MarkerSize',20,'LineWidth',5);
for i=1:p.nr
    u = ['r',num2str(i)];
    ra.(u) = line('color','r','LineWidth',2);
end
sz = 1/2*p.sr;
arrw = [-sz,0;sz,0;sz/2,sz/2;sz,0;sz/2,-sz/2];

if p.vidFlag==1
    if ~exist('../video', 'dir')
        mkdir('../video');
    end
    v = VideoWriter(p.vidDir);
    open(v)
end

tic
while toc<p.tE/s
    tQ = toc*s;
    
    [~,it] = min(abs(T-tQ));
    r = X.r(it,:);
    th = X.th(it,:);
    c = interp1(T,X.c,tQ);
    
    r = reshape(r,[p.nr 2]);
    c = reshape(c,[p.ny p.nx]);
    
    set(photo,'CData',c)
    %set(rants,'XData',r(:,1)/p.L*p.nx,'YData',(1-r(:,2)/p.L)*p.ny)
    
    for i=1:p.nr
        u = ['r',num2str(i)];
        A = [cos(th(i)),-sin(th(i));sin(th(i)),cos(th(i))];
        a = (A*arrw')';
        set(ra.(u),'XData',(r(i,1)+a(:,1))/p.L*p.nx,'YData',(1-(r(i,2)+a(:,2))/p.L)*p.ny);
    end
    
    drawnow;
    
    if p.vidFlag==1
        frame = getframe(gcf);
        writeVideo(v,frame);
    end
end
if p.vidFlag==1
    close(v)
end
end