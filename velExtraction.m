function [V, errV, F]=velExtraction(z, x, y, rx, ry, sampf, inFreq, freqSpace)
% Extracts velocity from z, rx, and ry
%
% [V, errV, F]=velExtraction(z, x, y, rx, ry, sampf, inFreq)

iter=length(freqSpace);

z=z-mean(z);
x=x-mean(x);
y=y-mean(y);
rx=rx-mean(rx);
ry=ry-mean(ry);

F=[];
V=[];
errV=[];
deltaXHist=[];
deltaYHist=[];

for a=1:iter-1
%     if max(a==inFreq)==1
        %% Data Crunching  
        % Bandpass filtering to get data into frequency bins
        freq=freqSpace(a);
        freqStep=(freqSpace(a+1)-freqSpace(a))/2;
        F=[F; freq];

        [bb,aa]= butter(3,[2*(freq-freqStep)/sampf 2*(freq+freqStep)/sampf],'bandpass');

        filtZ=filter(bb,aa,z);
        filtX=filter(bb,aa,x);
        filtY=filter(bb,aa,y);
        filtRX=filter(bb,aa,rx);
        filtRY=filter(bb,aa,ry);
        
        fitErr=[];
        shiftMax=ceil(sampf/freq);
        minFitErr=inf;
        minShift=0;
        for shift=(-shiftMax+1:shiftMax-1)
        
            hilRX=hilbert(filtRX(shiftMax:end-shiftMax));
            hilRY=hilbert(filtRY(shiftMax:end-shiftMax));
            hilX=hilbert(filtX(shiftMax+shift:end-shiftMax+shift));
            hilY=hilbert(filtY(shiftMax+shift:end-shiftMax+shift));
            hilZ=hilbert(filtZ(shiftMax:end-shiftMax));       

            cutLevel=0.1*max(abs(hilZ));

            cutIndex=find(abs(hilZ)>cutLevel);

            hilZ=hilZ(cutIndex);
            hilX=hilX(cutIndex);
            hilY=hilY(cutIndex);
            hilRY=hilRY(cutIndex);
            hilRX=hilRX(cutIndex);

            hilX=detrend(abs(hilX),'constant');
            hilY=detrend(abs(hilY),'constant');
            hilRX=detrend(abs(hilRX),'constant');
            hilRY=detrend(abs(hilRY),'constant');
            hilZ=detrend(abs(hilZ),'constant');

            fitX=[hilRX, hilY]';
            fitY=[hilRY, hilX]';
            
            wX=inv(fitX*fitX')*fitX*hilZ;
            wY=inv(fitY*fitY')*fitY*hilZ; 
            
            fitErr=[fitErr; sum((hilZ-(wX(1).*hilRX+wX(2).*hilX)).^2)];
            if sum((hilZ-(wX(1).*hilRX+wX(2).*hilX)).^2) < minFitErr
                minFitErr=sum((hilZ-(wX(1).*hilRX+wX(2).*hilX)).^2);
                minShift=shift;
            end
        end
        
        hilRX=hilbert(filtRX(shiftMax:end-shiftMax));
        hilRY=hilbert(filtRY(shiftMax:end-shiftMax));
        hilX=hilbert(filtX(shiftMax+minShift:end-shiftMax+minShift));
        hilY=hilbert(filtY(shiftMax+minShift:end-shiftMax+minShift));
        hilZ=hilbert(filtZ(shiftMax:end-shiftMax));       

        cutLevel=0.01*max(abs(hilZ));

        cutIndex=find(abs(hilZ)>cutLevel);

        hilZ=hilZ(cutIndex);
        hilX=hilX(cutIndex);
        hilY=hilY(cutIndex);
        hilRY=hilRY(cutIndex);
        hilRX=hilRX(cutIndex);

        hilX=detrend(abs(hilX),'constant');
        hilY=detrend(abs(hilY),'constant');
        hilRX=detrend(abs(hilRX),'constant');
        hilRY=detrend(abs(hilRY),'constant');
        hilZ=detrend(abs(hilZ),'constant');

        fitX=[hilRX, hilY]';
        fitY=[hilRY, hilX]';

        wX=inv(fitX*fitX')*fitX*hilZ;
        wY=inv(fitY*fitY')*fitY*hilZ; 
                
        % Translational coupling
        deltaXHist=[deltaXHist wX(2)/wX(1)];
        deltaYHist=[deltaYHist wY(2)/wY(1)];
        
        vx=wX(1);
        vy=wY(1);
        
        err=sqrt(abs(inv(fitX*fitX')*[1e-9 1e-9]'));
        
        errV=[errV err(1)];
        V=[V sqrt(vx^2+vy^2)];
   
%     end
    
end