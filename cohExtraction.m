function [C, err, F]=cohExtraction(signal1, signal2, sampf)
% Extracts coherence between two signals
%
% [C, err, F]=cohExtraction(signal1, signal2, sampf)

startFreq=0.01;
freqStep=.005;
endFreq=1;

iter=floor((endFreq-startFreq)/freqStep);

signal1=signal1-mean(signal1);
signal2=signal2-mean(signal2);
F=(startFreq:freqStep:endFreq);

C=[];
err=[];

for a=0:iter
    
    %% Data Crunching  
    % Bandpass filtering to get data into frequency bins
    freq=(startFreq+a*freqStep);
    [bb,aa]= butter(3,[2*((a-1/2)*freqStep+startFreq)/sampf 2*((a+1/2)*freqStep+startFreq)/sampf],'bandpass');

    filtSignal1=filter(bb,aa,signal1);
    filtSignal2=filter(bb,aa,signal2);

    fitLength=floor(2/(freq/sampf));
    
        
    temp=zeros(floor(length(filtSignal1)/fitLength)-2,1);
    
    parfor j=1:floor(length(filtSignal1)/fitLength)-2

        cut1=filtSignal1(j*fitLength:(j+1)*fitLength);
        cut2=filtSignal2(j*fitLength:(j+1)*fitLength);
  
        crossCor=sqrt(max(abs(xcorr(cut1,cut2))).^2./(max(abs(xcorr(cut1))).*max(abs(xcorr(cut2)))));
        
        temp(j)=crossCor;

    end
        
    C=[C mean(temp)];
    err=[err std(temp)/sqrt(length(temp))];
    
end
