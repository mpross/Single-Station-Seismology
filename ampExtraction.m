function [A, err, F]=ampExtraction(signal,sampf)
% Extracts complex amplitudes of signal
%
% [A, err, F]=ampExtraction(signal,sampF)

startFreq=0.01;
freqStep=.005;
endFreq=1;

iter=floor((endFreq-startFreq)/freqStep);

signal=signal-mean(signal);
F=(startFreq:freqStep:endFreq);

A=[];
err=[];

for a=0:iter
    
    temp=zeros(floor(length(filtSignal)/fitLength)-2,1);
    
    %% Data Crunching  
    % Bandpass filtering to get data into frequency bins
    freq=(startFreq+a*freqStep);
    [bb,aa]= butter(3,[2*((a-1/2)*freqStep+startFreq)/sampf 2*((a+1/2)*freqStep+startFreq)/sampf],'bandpass');

    filtSignal=filter(bb,aa,signal);

    fitLength=floor(1/(freq/sampf)/4);
    parfor j=1:floor(length(filtSignal)/fitLength)-2

        tim=(j*fitLength:(j+1)*fitLength)'./sampf;
        cut=filtSignal(j*fitLength:(j+1)*fitLength);
        
        x=[sin(2*pi*freq*tim), cos(2*pi*freq*tim)];
        
        w=cut'*x/(x'*x);
        
        temp(j)=abs(w(2))+abs(w(1))*1i
        
    end
    
    A=[A mean(temp)];
    err=[err std(temp)/sqrt(length(temp))];
    
end
