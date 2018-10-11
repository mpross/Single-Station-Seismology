function [bestPar,bestDispers]=dispersionFit(obsFreq,obsDispers,layers)
%% Fit
t0=cputime;
bestPar=[];
bestDispers=zeros(size(obsDispers));
N=1e4;
iter=1e2;
nDispers=zeros([N,length(obsFreq)]);
errSeries=[];
%%
err=zeros([N,1]);
%Assumes two layers and depth is equal to wavelength   
if(layers==2)
    %First layer
    vP1=rand([N,1])*3e3;
    vS1=rand([N,1])*3e3; %must be between 0 and /srt(2) Landau
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
    vS1=rand([N,1])*3e3; %must be between 0 and /srt(2) Landau
    d1=rand([N,1])*5e3;

    %Second layer
    vP2=rand([N,1])*3e3;
    vS2=rand([N,1])*3e3;
    d2=rand([N,1])*5e3;

    %Third layer
    vP3=rand([N,1])*3e3;
    vS3=rand([N,1])*3e3;
end
figure(15)
subplot(2,2,1)
plot1=scatter3(vP1,vS1,err,16,err,'filled');
ylabel('Layer 1 S-wave Velocity')
xlabel('Layer 1 P-wave Velocity')
zlabel('Error')
subplot(2,2,2)
plot2=scatter3(vP2,vS2,err,16,err,'filled');
ylabel('Layer 2 S-wave Velocity')
xlabel('Layer 2 P-wave Velocity')
zlabel('Error')
subplot(2,2,3)
plot3=scatter3(vP3,vS3,err,16,err,'filled');
ylabel('Layer 3 S-wave Velocity')
xlabel('Layer 3 S-wave Velocity')
zlabel('Error')
subplot(2,2,4)
plot4=scatter3(d1,d2,err,16,err,'filled');
ylabel('Layer 1 Thickness')
xlabel('Layer 2 Thickness')
zlabel('Error')
%%

figure(19)
plot5=plot(obsFreq,bestDispers);
ylabel('Velocity (m/s)')
xlabel('Frequency (Hz)')
set(plot2,'LineWidth',1.5);
set(gca,'FontSize',16);
set(plot5,'MarkerSize',16);

figure(20)
plot6=semilogy(errSeries);
ylabel('Minimum Error')
xlabel('Iteration')
set(plot2,'LineWidth',1.5);
set(gca,'FontSize',16);
set(plot5,'MarkerSize',16);

for k=(1:iter)
    for n=(1:N)

        coeffs=[1/vS1(n)^6 0 -8/vS1(n)^4 0 8/vS1(n)^2*(3-2*vS1(n)^2/vP1(n)^2) 0 -16*(1-vS1(n)^2/vP1(n)^2)];
        rots=roots(coeffs);
        fitV1=max(rots(find(and(and(and(imag(rots)==0,rots>0),rots<vP1(n)),vS3(n)/vP3(n)<1/sqrt(2)))));

        coeffs=[1/vS2(n)^6 0 -8/vS2(n)^4 0 8/vS2(n)^2*(3-2*vS2(n)^2/vP2(n)^2) 0 -16*(1-vS2(n)^2/vP2(n)^2)];
        rots=roots(coeffs);
        fitV2=max(rots(find(and(and(and(imag(rots)==0,rots>0),rots<vP2(n)),vS3(n)/vP3(n)<1/sqrt(2)))));
        if(layers==3)
            coeffs=[1/vS3(n)^6 0 -8/vS3(n)^4 0 8/vS3(n)^2*(3-2*vS3(n)^2/vP3(n)^2) 0 -16*(1-vS3(n)^2/vP3(n)^2)];
            rots=roots(coeffs);
            fitV3=max(rots(find(and(and(and(imag(rots)==0,rots>0),rots<vP3(n)),vS3(n)/vP3(n)<1/sqrt(2)))));
        end

        fitDispers=[];
        for j=(1:length(obsFreq))
            if(layers==2)
                % Two layers square attenuation
                coeffs=[1 -fitV2 (fitV2-fitV1)*d1(n)*obsFreq(j)];
                % Two layers linear attenuation
    %             coeffs=[1 -1/2*d1^2*obsFreq(j)^2*(fitV1-fitV2) -obsFreq(j)^3*d1(fitV1-fitV2) obsFreq(j)^4*fitV2/2];
            elseif(layers==3)
                % Three layers
                coeffs=[1 -fitV3 (((fitV2-fitV1)*d1(n)+(fitV3-fitV2)*d2(n))*obsFreq(j))];
            else
                'Invalid layer number.'
                break
            end
            rots=roots(coeffs);        
            if isempty(max(rots(find(and(rots>0,imag(rots)==0)))))
                coeffs=[1 -fitV2 (fitV2-fitV1)*d1(n)*obsFreq(j)];
                rots=roots(coeffs);
                if isempty(max(rots(find(and(rots>0,imag(rots)==0)))))
                    fitDispers=[fitDispers; nan];
                else
                    fitDispers=[fitDispers; max(rots(find(and(rots>0,imag(rots)==0))))];
                end
            else
                fitDispers=[fitDispers; max(rots(find(and(rots>0,imag(rots)==0))))];
            end
        end        
        nDispers(n,:)=fitDispers;
%         fitDispers=fitDispers(find(not(isnan(fitDispers))));
%         cutDispers=obsDispers(find(not(isnan(fitDispers))));
        err(n)=sum((fitDispers-obsDispers).^2);        
    end
    w1=0;
    w2=0.05;
    sigma=iter/k*10;
%     
%     vP1=randomizeNans(vP1);
%     vS1=randomizeNans(vS1);
%     d1=randomizeNans(d1);
%     vP2=randomizeNans(vP2);
%     vS2=randomizeNans(vS2);
%     d2=randomizeNans(d2);
%     vP3=randomizeNans(vP3);
%     vS3=randomizeNans(vS3);
    
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
%         else
%             vP1(n)=vP1(n)+sigma*rand;
%             vP2(n)=vP2(n)+sigma*rand;
%             vP3(n)=vP3(n)+sigma*rand;
%             vS1(n)=vS1(n)+sigma*rand;
%             vS2(n)=vS2(n)+sigma*rand;
%             vS3(n)=vS3(n)+sigma*rand;
%             d1(n)=d1(n)+sigma*rand;
%             d2(n)=d2(n)+sigma*rand;
        end
    end
    
    bestN=min(find(err==min(err)));
    bestDispers=nDispers(bestN,:);
    errSeries=[errSeries; err(bestN)];
    bestPar=[vP1(bestN) vS1(bestN) d1(bestN) vP2(bestN) vS2(bestN) d2(bestN) vP3(bestN) vS3(bestN)];
    if mod(k,1)==0
        
        set(plot1,'XData',vP1)
        set(plot1,'YData',vS1)
        set(plot1,'ZData',err)
        set(plot1,'CData',err)
        
        set(plot2,'XData',vP2)
        set(plot2,'YData',vS2)
        set(plot2,'ZData',err)
        set(plot2,'CData',err)
        
        set(plot3,'XData',vP3)
        set(plot3,'YData',vS3)
        set(plot3,'ZData',err)
        set(plot3,'CData',err)
        
        set(plot4,'XData',d1)
        set(plot4,'YData',d2)
        set(plot4,'ZData',err)
        set(plot4,'CData',err)
        
        set(plot5,'XData',obsFreq)
        set(plot5,'YData',bestDispers)     
        
        if(k==1)
            figure(20)
            plot6=plot(errSeries);
            ylabel('Minimum Error')
            xlabel('Iteration')
            set(plot2,'LineWidth',1.5);
            set(gca,'FontSize',16);
            set(plot5,'MarkerSize',16);
        end
        set(plot6,'YData',errSeries)	
        figure(20); 
        set(gca,'YScale','log')
        
        refreshdata
        drawnow
        
        figure(15)
        F(k)=getframe(gcf);
        
        disp([num2str(k/iter*100) ' % done'])
        disp([num2str((iter-k)/k*(cputime-t0)) ' s left'])
    end
end
v = VideoWriter('dispersionSwarm.avi');
v.FrameRate=2;
open(v);
writeVideo(v,F);
close(v);

if(isempty(bestDispers))
        'Can not find solution'
        bestPar=[0 0 0 0 0 0 0 0];
        bestDispers=nan*obsFreq;
end
