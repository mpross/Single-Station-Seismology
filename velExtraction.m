function [V, errV, A, errA, F]=velExtraction(z, x, y, rx,ry,sampf, inFreq)
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
F=(startFreq:freqStep:endFreq);

V=[];
A=[];
errV=[];
errA=[];
angHist=[];

for a=0:iter
    
    %% Data Crunching  
    % Bandpass filtering to get data into frequency bins
    freq=(startFreq+a*freqStep);
    
    [bb,aa]= butter(3,[2*((a-1/2)*freqStep+startFreq)/sampf 2*((a+1/2)*freqStep+startFreq)/sampf],'bandpass');

    filtZ=filter(bb,aa,z);
    filtX=filter(bb,aa,x);
    filtY=filter(bb,aa,y);
    filtRX=filter(bb,aa,rx);
    filtRY=filter(bb,aa,ry);

%     figure(150)
%     plot((1:length(filtZ)),filtZ/1e3,(1:length(filtZ)),filtX/1e3,(1:length(filtZ)),filtY/1e3)
%     figure(151)
%     plot((1:length(filtZ)),filtRX+3.8e-4*filtY,(1:length(filtZ)),filtRY+1.5e-4*filtX,...
%         (1:length(filtZ)),filtRX,(1:length(filtZ)),filtRY)

    fitLength=floor(1/(freq/sampf));

    tempZ=zeros(floor(length(z)/fitLength)-2,1);
    tempX=zeros(floor(length(x)/fitLength)-2,1);
    tempY=zeros(floor(length(y)/fitLength)-2,1);
    tempRX=zeros(floor(length(rx)/fitLength)-2,1);
    tempRY=zeros(floor(length(ry)/fitLength)-2,1);

    tempV=zeros(floor(length(ry)/fitLength)-2,1);
    tempA=zeros(floor(length(ry)/fitLength)-2,1);
%     
%     crossRX=zeros(floor(length(ry)/fitLength)-2,1);
%     crossRY=zeros(floor(length(ry)/fitLength)-2,1);

    for j=1:floor(length(filtZ)/fitLength)-2

        tim=(j*fitLength:(j+1)*fitLength)'./sampf;
        cutZ=filtZ(j*fitLength:(j+1)*fitLength);
        cutX=filtX(j*fitLength:(j+1)*fitLength);
        cutY=filtY(j*fitLength:(j+1)*fitLength);
        cutRX=filtRX(j*fitLength:(j+1)*fitLength);
        cutRY=filtRY(j*fitLength:(j+1)*fitLength);
%         
%         crossRX(j)=sum(cutRX.*cutZ)/sum(cutRX.*cutRX);
%         crossRY(j)=sum(cutRY.*cutZ)/sum(cutRY.*cutRY);
        if max(abs(cutZ))>=0.9*max(abs(filtZ))
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

            if and(abs(mod(angle(tempZ(j))-angle(tempRX(j)),pi))<45*pi/180,abs(mod(angle(tempZ(j))-angle(tempRY(j)),pi))<45*pi/180)
                tempV(j)=abs(tempZ(j))./sqrt(abs(tempRX(j)).^2+abs(tempRY(j)).^2);

                phi=angle(tempZ(j));
                rotRX=real(tempRX(j))*cos(phi)+imag(tempRX(j))*sin(phi)...
                    +i*(real(tempRX(j))*-sin(phi)+imag(tempRX(j))*cos(phi));
                rotRY=real(tempRY(j))*cos(phi)+imag(tempRY(j))*sin(phi)...
                    +i*(real(tempRY(j))*-sin(phi)+imag(tempRY(j))*cos(phi));

                tempA(j)=atan2(sign(real(rotRY))*abs(tempRY(j)),sign(real(rotRX))*abs(tempRX(j)));
%               if or(crossRY>=0.9*tempV(j),crossRX>=0.9*tempV(j))
%                 tempA(j)=atan2(sign(crossRY)*abs(tempRY(j)),sign(crossRX)*abs(tempRX(j)));
            end
%             end
        end
    end

%     figure(10)
%     plot(1:floor(length(filtZ)/fitLength)-2, angle(tempZ)*180/pi,...
%         1:floor(length(filtZ)/fitLength)-2, angle(tempX)*180/pi,...
%         1:floor(length(filtZ)/fitLength)-2, angle(tempY)*180/pi,...
%         1:floor(length(filtZ)/fitLength)-2, angle(tempRX)*180/pi,...
%         1:floor(length(filtZ)/fitLength)-2, angle(tempRY)*180/pi)
%     
%     figure(11)
%     plot(1:floor(length(filtZ)/fitLength)-2, abs(tempZ)/1e3,...
%         1:floor(length(filtZ)/fitLength)-2, abs(tempX)/1e3,...
%         1:floor(length(filtZ)/fitLength)-2, abs(tempY)/1e3,...
%         1:floor(length(filtZ)/fitLength)-2, abs(tempRX),...
%         1:floor(length(filtZ)/fitLength)-2, abs(tempRY));
%     ,...
%         1:floor(length(filtZ)/fitLength)-2, abs(tempRX+2.91e-4*tempY),...
%         1:floor(length(filtZ)/fitLength)-2, abs(tempRY+1.27e-4*tempX))
%     %
    
    
    if max(a==inFreq)==1
        angHist=[angHist; tempA];
%         
    end

    V=[V mean(nonzeros(tempV))];
    errV=[errV std(nonzeros(tempV))/sqrt(length(nonzeros(tempV)))];
    A=[A mean(nonzeros(tempA))];
    errA=[errA std(nonzeros(tempA))/sqrt(length(nonzeros(tempA)))];
    
    
end

fig8=figure(8);
polarhistogram(nonzeros(angHist)-225*pi/180,20,'Normalization','probability')   
legend('','Mexico','Fiji','Venezula','Peru','NewZealand','Canada','Iceland')
