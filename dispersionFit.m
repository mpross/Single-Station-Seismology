function [bestPar,bestDispers,frames]=dispersionFit(obsFreq,obsDispers,layers)
% Iteratively fits observed Rayleigh dispersion curves to a set of P and S wave
% velocities and depths for two or three layer models. Also outputs frames
% of the fit to be used to watch how well the fitting worked using
% movie(frames).
%
% [bestPar,bestDispers,frames]=dispersionFit(obsFreq,obsDispers,layers)

% Initial variable set up
t0=cputime;
bestPar=[];
bestDispers=zeros(size(obsDispers));
N=1e3;
iter=2e1;
nDispers=zeros([N,length(obsFreq)]);
errSeries=[];
err=zeros([N,1]);
frames=[];

% Initialize randomly sampled parameter space
if(layers==2)
    %First layer
    vP1=rand([N,1])*3e3;
    vS1=rand([N,1])*3e3;
    d1=rand([N,1])*5e3;

    %Second layer
    vP2=rand([N,1])*3e3;
    vS2=rand([N,1])*3e3;
    d2=zeros([N,1]);

    %Third layer
    vP3=zeros([N,1]);
    vS3=zeros([N,1]);
elseif(layers==3)
    %First layer
    vP1=rand([N,1])*3e3;
    vS1=rand([N,1])*3e3; 
    d1=rand([N,1])*5e3;

    %Second layer
    vP2=rand([N,1])*3e3;
    vS2=rand([N,1])*3e3;
    d2=rand([N,1])*5e3;

    %Third layer
    vP3=rand([N,1])*3e3;
    vS3=rand([N,1])*3e3;
end

% Create plots
% % % 
% figure(15)
% ax1=subplot(2,2,1);
% plot1=scatter3(vP1,vS1,err,16,err,'filled');
% ylabel('Layer 1 S-wave Velocity')
% xlabel('Layer 1 P-wave Velocity')
% zlabel('Error')
% ax2=subplot(2,2,2);
% plot2=scatter3(vP2,vS2,err,16,err,'filled');
% ylabel('Layer 2 S-wave Velocity')
% xlabel('Layer 2 P-wave Velocity')
% zlabel('Error')
% ax3=subplot(2,2,3);
% plot3=scatter3(vP3,vS3,err,16,err,'filled');
% ylabel('Layer 3 S-wave Velocity')
% xlabel('Layer 3 S-wave Velocity')
% zlabel('Error')
% ax4=subplot(2,2,4);
% plot4=scatter3(d1,d2,err,16,err,'filled');
% ylabel('Layer 1 Thickness')
% xlabel('Layer 2 Thickness')
% zlabel('Error')
% 
% figure(19)
% plot5=plot(obsFreq,bestDispers);
% ylabel('Velocity (m/s)')
% xlabel('Frequency (Hz)')
% set(plot2,'LineWidth',1.5);
% set(gca,'FontSize',16);
% set(plot5,'MarkerSize',16);
% 
% figure(20)
% plot6=semilogy(errSeries);
% ylabel('Minimum Error')
% xlabel('Iteration')
% set(plot2,'LineWidth',1.5);
% set(gca,'FontSize',16);
% set(plot5,'MarkerSize',16);

%% Main loop
for k=(1:iter)
    for n=(1:N)
        fitDispers=dispersionCalc(vP1(n),vS1(n),d1(n),vP2(n),vS2(n),d2(n),vP3(n),vS3(n),obsFreq,layers);
        nDispers(n,:)=fitDispers;        
        err(n)=sum((fitDispers-obsDispers).^2);        
    end
    w1=0.5;
    w2=0.01;
    sigma=(iter-k)*10;
    
    for n=(1:N)
        dist=sqrt((vP1(n)-vP1).^2+(vS1(n)-vS1).^2+(d1(n)-d1).^2 ...
             +(vP2(n)-vP2).^2+(vS2(n)-vS2).^2+(d2(n)-d2).^2 ...
             +(vP3(n)-vP3).^2+(vS3(n)-vS3).^2);
         if not(isempty(find(and(dist>0,err<err(n)))))
            nearestN=min(find(dist==min(dist(find(and(dist>0,err<err(n)))))));
            globalN=min(find(err==min(err)));            

            vP1(n)=vP1(n)+w1*(vP1(nearestN)-vP1(n))+w2*(vP1(globalN)-vP1(n))+sigma*rand;
            vP2(n)=vP2(n)+w1*(vP2(nearestN)-vP2(n))+w2*(vP2(globalN)-vP2(n))+sigma*rand;
            vP3(n)=vP3(n)+w1*(vP3(nearestN)-vP3(n))+w2*(vP3(globalN)-vP3(n))+sigma*rand;
            vS1(n)=vS1(n)+w1*(vS1(nearestN)-vS1(n))+w2*(vS1(globalN)-vS1(n))+sigma*rand;
            vS2(n)=vS2(n)+w1*(vS2(nearestN)-vS2(n))+w2*(vS2(globalN)-vS2(n))+sigma*rand;
            vS3(n)=vS3(n)+w1*(vS3(nearestN)-vS3(n))+w2*(vS3(globalN)-vS3(n))+sigma*rand;
            d1(n)=d1(n)+w1*(d1(nearestN)-d1(n))+w2*(d1(globalN)-d1(n))+sigma*rand;
            d2(n)=d2(n)+w1*(d2(nearestN)-d2(n))+w2*(d2(globalN)-d2(n))+sigma*rand;
         end
    end
    
    bestN=min(find(err==min(err)));
    bestDispers=nDispers(bestN,:);
    errSeries=[errSeries; err(bestN)];
    bestPar=[vP1(bestN) vS1(bestN) d1(bestN) vP2(bestN) vS2(bestN) d2(bestN) vP3(bestN) vS3(bestN)];
%     if mod(k,1)==0
%         
%         set(plot1,'XData',vP1)
%         set(plot1,'YData',vS1)
%         set(plot1,'ZData',err)
%         set(plot1,'CData',err)
%         
%         set(plot2,'XData',vP2)
%         set(plot2,'YData',vS2)
%         set(plot2,'ZData',err)
%         set(plot2,'CData',err)
%         
%         set(plot3,'XData',vP3)
%         set(plot3,'YData',vS3)
%         set(plot3,'ZData',err)
%         set(plot3,'CData',err)
%         
%         set(plot4,'XData',d1)
%         set(plot4,'YData',d2)
%         set(plot4,'ZData',err)
%         set(plot4,'CData',err)
%         
%         figure(15)
%         axis([ax1 ax2 ax3 ax4],[0 5e3 0 5e3 0 1e8])
%         
%         set(plot5,'XData',obsFreq)
%         set(plot5,'YData',bestDispers)     
%         
%         if(k==1)
%             figure(20)
%             plot6=plot(errSeries);
%             ylabel('Minimum Error')
%             xlabel('Iteration')
%             set(plot2,'LineWidth',1.5);
%             set(gca,'FontSize',16);
%             set(plot5,'MarkerSize',16);
%         end
%         set(plot6,'YData',errSeries)	
%         figure(20); 
%         set(gca,'YScale','log')
%         
%         refreshdata
%         drawnow
%         
%         figure(15)
%         frames(k)=getframe(gcf);
%         
%         disp([num2str(k/iter*100) ' % done'])
%         disp([num2str((iter-k)/k*(cputime-t0)) ' s left'])
%     end
end
% v = VideoWriter('dispersionSwarm.avi');
% v.FrameRate=2;
% open(v);
% writeVideo(v,frames);
% close(v);

if(isempty(bestDispers))
        'Can not find solution'
        bestPar=[0 0 0 0 0 0 0 0];
        bestDispers=nan*obsFreq;
end
