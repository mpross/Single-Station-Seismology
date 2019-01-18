function [V, err, F]=velExtraction(z,rx,ry,sampf)
% Extracts velocity from z, rx, and ry
%
% [V, err, F]=ampExtraction(z,rx,ry,sampF)

startFreq=0.01;
freqStep=.005;
endFreq=1;

iter=floor((endFreq-startFreq)/freqStep);

z=z-mean(z);
rx=rx-mean(z);
ry=ry-mean(z);
F=(startFreq:freqStep:endFreq);

V=[];
err=[];

for a=0:iter
    
    %% Data Crunching  
    % Bandpass filtering to get data into frequency bins
    freq=(startFreq+a*freqStep);
    [bb,aa]= butter(3,[2*((a-1/2)*freqStep+startFreq)/sampf 2*((a+1/2)*freqStep+startFreq)/sampf],'bandpass');

    filtZ=filter(bb,aa,z);
    filtRX=filter(bb,aa,rx);
    filtRY=filter(bb,aa,ry);

    fitLength=floor(1/(freq/sampf)/4);
        
    tempZ=zeros(floor(length(z)/fitLength)-2,1);
    tempRX=zeros(floor(length(rx)/fitLength)-2,1);
    tempRY=zeros(floor(length(ry)/fitLength)-2,1);
    tempV=zeros(floor(length(ry)/fitLength)-2,1);
    
    parfor j=1:floor(length(filtZ)/fitLength)-2

        tim=(j*fitLength:(j+1)*fitLength)'./sampf;
        cutZ=filtZ(j*fitLength:(j+1)*fitLength);
        cutRX=filtRX(j*fitLength:(j+1)*fitLength);
        cutRY=filtRY(j*fitLength:(j+1)*fitLength);
        
        x=[sin(2*pi*freq*tim), cos(2*pi*freq*tim)];
        
        wZ=cutZ'*x/(x'*x);
        wRX=cutRX'*x/(x'*x);
        wRY=cutRY'*x/(x'*x);
        
        tempZ(j)=abs(wZ(2))+abs(wZ(1))*1i;
        tempRX(j)=abs(wRX(2))+abs(wRX(1))*1i;
        tempRY(j)=abs(wRY(2))+abs(wRY(1))*1i;
        if or(abs(angle(tempZ(j))-angle(tempRX(j)))<10,abs(angle(tempZ(j))-angle(tempRY(j)))<10)
            
            tempV(j)=abs(tempZ(j))./sqrt(abs(tempRX(j)).^2+abs(tempRY(j)).^2);
            
        end
    end
    
    V=[V mean(tempV)];
    err=[err std(tempV)];
    
end
