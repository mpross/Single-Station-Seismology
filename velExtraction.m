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

        fitLength=floor(1/(freq/sampf));

        tempZ=zeros(floor(length(z)/fitLength)-2,1);
        tempX=zeros(floor(length(x)/fitLength)-2,1);
        tempY=zeros(floor(length(y)/fitLength)-2,1);
        tempRX=zeros(floor(length(rx)/fitLength)-2,1);
        tempRY=zeros(floor(length(ry)/fitLength)-2,1);

        tempV=zeros(floor(length(ry)/fitLength)-2,1);
        tempA=zeros(floor(length(ry)/fitLength)-2,1);

        for j=1:floor(length(filtZ)/fitLength)-2

            tim=(j*fitLength:(j+1)*fitLength)'./sampf;
            cutZ=filtZ(j*fitLength:(j+1)*fitLength);
            cutX=filtX(j*fitLength:(j+1)*fitLength);
            cutY=filtY(j*fitLength:(j+1)*fitLength);
            cutRX=filtRX(j*fitLength:(j+1)*fitLength);
            cutRY=filtRY(j*fitLength:(j+1)*fitLength);

            if max(abs(cutZ))>=0.5*max(abs(filtZ))
                fitx=[sin(2*pi*freq*tim), cos(2*pi*freq*tim)];

                wZ=cutZ'*fitx/(fitx'*fitx);
                wX=cutX'*fitx/(fitx'*fitx);
                wY=cutY'*fitx/(fitx'*fitx);
                wRX=cutRX'*fitx/(fitx'*fitx);
                wRY=cutRY'*fitx/(fitx'*fitx);

                tempZ(j)=wZ(2)+wZ(1)*1i;
                tempX(j)=wX(2)+wX(1)*1i;
                tempY(j)=wY(2)+wY(1)*1i;

                % Extra term is translational coupling subtraction
                tempRX(j)=wRX(2)+wRX(1)*1i+3.8e-4*tempY(j);
                tempRY(j)=wRY(2)+wRY(1)*1i+1.5e-4*tempX(j);

                tempV(j)=abs(tempZ(j))./sqrt(abs(tempRX(j)).^2+abs(tempRY(j)).^2);
                tempA(j)=atan2(abs(tempRY(j)),abs(tempRX(j)));

            end
        end    
    
        angHist=[angHist; tempA];
        V=[V mean(nonzeros(tempV))];
        errV=[errV std(nonzeros(tempV))/sqrt(length(nonzeros(tempV)))];
        A=[A mean(nonzeros(tempA))];
        errA=[errA std(nonzeros(tempA))/sqrt(length(nonzeros(tempA)))];
    
    end
    
end

fig8=figure(8);
polarhistogram(nonzeros(angHist)+angOffset*pi/180,'Normalization','probability')   
legend('','Mexico','Fiji','Venezula','Peru','NewZealand','Canada','Iceland')
