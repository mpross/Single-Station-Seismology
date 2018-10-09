function [bestPar,bestDispers]=dispersionFit(obsFreq,obsDispers,layers)
%% Fit
bestPar=[];
bestDispers=zeros(size(obsDispers));
N=1e3;
iter=1e3;
nDispers=zeros([N,length(obsFreq)]);
%%
err=zeros([N,1]);
%Assumes two layers and depth is equal to wavelength   
if(layers==2)
    %First layer
    vP1=rand([N,1])*3e3;
    vS1=rand([N,1])*3e3; %must be between 0 and /srt(2) Landau
    d1=rand([N,1])*10e3;

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
    d1=rand([N,1])*10e3;

    %Second layer
    vP2=rand([N,1])*3e3;
    vS2=rand([N,1])*3e3;
    d2=rand([N,1])*10e3;

    %Third layer
    vP3=rand([N,1])*3e3;
    vS3=rand([N,1])*3e3;
end
figure(15)
plot1=scatter3(vP1,vS1,err,16,err,'filled')
figure(16)
plot2=scatter3(vP2,vS2,err,16,err,'filled')
figure(17)
plot3=scatter3(vP3,vS3,err,16,err,'filled')
figure(18)
plot4=scatter3(d1,d2,err,16,err,'filled')
%%

% figure(14)
% plot2=plot(obsFreq,bestDispers);
% ylabel('Velocity (m/s)')
% xlabel('Frequency (Hz)')
% set(plot2,'LineWidth',1.5);
% set(gca,'FontSize',16);
% set(plot2,'MarkerSize',16);

for k=(0:iter)
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
        err(n)=sum((fitDispers-obsDispers).^2);        
    end
    w=0.5;
    for n=(1:N)
        dist=sqrt((vP1(n)-vP1).^2+(vS1(n)-vS1).^2+(d1(n)-d1).^2 ...
            +(vP2(n)-vP2).^2+(vS2(n)-vS2).^2+(d2(n)-d2).^2 ...
            +(vP3(n)-vP3).^2+(vS3(n)-vS3).^2);
        if not(isempty(find(and(dist>0,err<err(n)))))
            nearestN=min(find(dist==min(dist(find(and(dist>0,err<err(n)))))));
            vP1(n)=vP1(n)+w*(vP1(nearestN)-vP1(n));
            vP2(n)=vP2(n)+w*(vP2(nearestN)-vP2(n));
            vP3(n)=vP3(n)+w*(vP3(nearestN)-vP3(n));
            vS1(n)=vS1(n)+w*(vS1(nearestN)-vS1(n));
            vS2(n)=vS2(n)+w*(vS2(nearestN)-vS2(n));
            vS3(n)=vS3(n)+w*(vS3(nearestN)-vS3(n));
            d1(n)=d1(n)+w*(d1(nearestN)-d1(n));
            d2(n)=d2(n)+w*(d2(nearestN)-d2(n));
        end
    end
    
    bestN=min(find(err==min(err)));
    bestDispers=nDispers(bestN,:);
    bestErr=err(bestN);
    bestPar=[vP1(bestN) vS1(bestN) d1(bestN) vP2(bestN) vS2(bestN) d2(bestN) vP3(bestN) vS3(bestN)];
    if mod(k,1)==0
        set(plot2,'XData',obsFreq)
        set(plot2,'YData',bestDispers)
        
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
        
        refreshdata
        drawnow
    end
end
if(isempty(bestDispers))
        'Can not find solution'
        bestPar=[0 0 0 0 0 0 0 0];
        bestDispers=nan*obsFreq;
end