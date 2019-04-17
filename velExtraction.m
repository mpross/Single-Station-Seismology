function [V, errV, F]=velExtraction(z, x, y, rx,ry,sampf, inFreq, angOffset)
% Extracts velocity from z, rx, and ry
%
% [V, err, F]=ampExtraction(z,rx,ry,sampF)

startFreq=0.02;
freqStep=0.02;
endFreq=1;

iter=floor((endFreq-startFreq)/freqStep);

z=z-mean(z);
x=x-mean(x);
y=y-mean(y);
rx=rx-mean(rx);
ry=ry-mean(ry);

F=[];
V=[];
A=[];
errV=[];
errA=[];
angHist=[];
deltaXHist=[];
deltaYHist=[];

for a=0:iter
%     if max(a==inFreq)==1
        %% Data Crunching  
        % Bandpass filtering to get data into frequency bins
        freq=(startFreq+a*freqStep);
        F=[F; freq];

        [bb,aa]= butter(3,[2*((a-1/2)*freqStep+startFreq)/sampf 2*((a+1/2)*freqStep+startFreq)/sampf],'bandpass');

        filtZ=filter(bb,aa,z);
        filtX=filter(bb,aa,x);
        filtY=filter(bb,aa,y);
        filtRX=filter(bb,aa,rx);
        filtRY=filter(bb,aa,ry);
        
        hilRX=hilbert(filtRX(100*sampf:end));
        hilRY=hilbert(filtRY(100*sampf:end));
        hilX=hilbert(filtX(100*sampf:end));
        hilY=hilbert(filtY(100*sampf:end));
        hilZ=hilbert(filtZ(100*sampf:end));       
        
        
%         tempV=abs(hilZ)./sqrt(abs(hilRX).^2+abs(hilRY).^2);
%         tempA=atan2(abs(hilRX),abs(hilRY));
        
%         cutLevel=0.1*max(abs(hilZ));
%         cutIndex=find(abs(hilZ)<cutLevel);
        
%         tempV(cutIndex)=0;
%         tempA(cutIndex)=0;
        
%         tim1=(1:length(hilZ))/sampf;
        
        
        
%         tempV=tempV(cutIndex);
%         tempA=tempA(cutIndex);
%         V=[V mean(nonzeros(tempV))];
%         errV=[errV std(nonzeros(tempV))/sqrt(length(nonzeros(tempV)))];

%         cutIndex=find(abs(hilZ)>cutLevel);
% %         
%         hilZ=hilZ(cutIndex);
%         hilX=hilX(cutIndex);
%         hilY=hilY(cutIndex);
%         hilRY=hilRY(cutIndex);
%         hilRX=hilRX(cutIndex);
%         
%         filtZ=filtZ(cutIndex);
%         filtX=filtX(cutIndex);
%         filtY=filtY(cutIndex);
%         filtRY=filtRY(cutIndex);
%         filtRX=filtRX(cutIndex);
        
%         phiRX=((unwrap(angle(hilRX)))-(unwrap(angle(hilZ))));
%         phiRY=((unwrap(angle(hilRY)))-(unwrap(angle(hilZ))));
%         rotRX=real(hilRX).*cos(phiRX)+imag(hilRX).*sin(phiRX)...
%                 +i.*(real(hilRX).*-sin(phiRX)+imag(hilRX).*cos(phiRX));
%         rotRY=real(hilRY).*cos(phiRY)+imag(hilRY).*sin(phiRY)...
%             +i.*(real(hilRY).*-sin(phiRY)+imag(hilRY).*cos(phiRY));
%         
%         tempV=abs(hilZ)./sqrt(abs(hilRX).^2+abs(hilRY).^2);
%         tempA=unwrap(atan2(sign(cos(phiRX)).*abs(hilRX),sign(cos(phiRY)).*abs(hilRY)));
%         
%         tim2=(1:length(hilZ))/sampf;
       
%         figure(100)
%         plot(tim2,tempV)
%         xlabel('Time')
%         ylabel('Velocity (m/s)')
%         
%         figure(101)
%         plot(tim2,tempA*180/pi)        
%         xlabel('Time')
%         ylabel('Azimuth')
%         
%         figure(102)
%         plot(tim2,sign(cos(phiRX)).*abs(hilRX),tim2,sign(cos(phiRY)).*abs(hilRY),tim2,abs(hilZ)/1000,tim2,tim2*0+cutLevel/1000)
%         xlabel('Time')     
%                 
%         figure(104)         
%         plot(tim2,((unwrap(angle(hilRX)))-(unwrap(angle(hilZ))))*180/pi) 
%         hold on
%         plot(tim2,((unwrap(angle(hilRY)))-(unwrap(angle(hilZ))))*180/pi) 
%         hold off
% 
        hilX=detrend(abs(hilX),'constant');
        hilY=detrend(abs(hilY),'constant');
        hilRX=detrend(abs(hilRX),'constant');
        hilRY=detrend(abs(hilRY),'constant');
        hilZ=detrend(abs(hilZ),'constant');
        
        fitX=[hilRX, hilX]';
        fitY=[hilRY, hilY]';
       
        wX=inv(fitX*fitX')*fitX*hilZ;
        wY=inv(fitY*fitY')*fitY*hilZ; 
        
%         figure(100)
%         plot(tim2, filtZ, tim2, wX'*fitX)%, tim2, 4e3*filtRX, tim2, filtX
        
        % Translational coupling
        deltaXHist=[deltaXHist wX(2)/wX(1)];
        deltaYHist=[deltaYHist wY(2)/wY(1)];
        
        vx=wX(1);
        vy=wY(1);
        
%         vx=inv(hilRX'*hilRX)*hilZ'*filtRX;
%         vy=inv(filtRY'*filtRY)*filtZ'*filtRY;
%         sqrt(vx^2+vy^2)
%         mean(abs(hilZ))./sqrt(mean(abs(hilRX)).^2+mean(abs(hilRY)).^2)  
        
%         V=[V mean(abs(hilZ))./sqrt(mean(abs(hilRX)).^2+mean(abs(hilRY)).^2)];
        errV=[errV 0];
        V=[V sqrt(vx^2+vy^2)];
        angHist=[angHist atan2(vy,vx)];
%         angHist=[angHist tempA'];
%         errV=[errV sqrt( (std(abs(hilZ))^2*(mean(abs(hilRX)).^2+mean(abs(hilRY)).^2).^2+...
%             (std(abs(hilRX))^2)*mean(abs(hilRX)).^2+...
%             (std(abs(hilRY))^2)*mean(abs(hilRY)).^2)./(mean(abs(hilRX)).^2+mean(abs(hilRY)).^2)^3 )/sqrt(length(hilZ))];
%     
%     end
    
end

fig8=figure(8);
polarhistogram(angHist+225/180*pi,20,'Normalization','probability')   
legend('','Mexico','Fiji','Venezula','Peru','NewZealand','Canada','Iceland','Peru','Peru')

figure(9)
histogram(log10(deltaXHist),20)