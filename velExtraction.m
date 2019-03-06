function [V, errV, A, errA, F]=velExtraction(z, x, y, rx,ry,sampf, inFreq, angOffset)
% Extracts velocity from z, rx, and ry
%
% [V, err, F]=ampExtraction(z,rx,ry,sampF)

startFreq=0.01;
freqStep=0.01;
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

for a=0:iter
    if max(a==inFreq)==1
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
        
        hilRX=hilbert(filtRX);
        hilRY=hilbert(filtRY);
        hilX=hilbert(filtX);
        hilY=hilbert(filtY);
        hilZ=hilbert(filtZ);
        
        tempV=abs(hilZ)./sqrt(abs(hilRX).^2+abs(hilRY).^2);
        tempA=atan2(abs(hilRY),abs(hilRX));
        
        cutIndex=find(abs(hilZ)>0.1*max(abs(hilZ)));
        tempA=tempA(cutIndex);
        tempV=tempV(cutIndex);
    
        angHist=[angHist; tempA];
        V=[V mean(nonzeros(tempV))];
        errV=[errV std(nonzeros(tempV))/sqrt(length(nonzeros(tempV)))];
        A=[A mean(nonzeros(tempA))];
        errA=[errA std(nonzeros(tempA))/sqrt(length(nonzeros(tempA)))];
    
    end
    
end

fig8=figure(8);
polarhistogram(nonzeros(angHist)+202.5*pi/180,'Normalization','probability')
legend('','Mexico','Fiji','Venezula','Peru','NewZealand','Canada','Iceland')
