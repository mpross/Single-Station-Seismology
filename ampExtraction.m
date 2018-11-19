function [A,F]=ampExtraction(signal,sampf)
% Extracts complex amplitudes of signal
%
% [F,A]=ampExtraction(signal,sampF)

startFreq=0.025;
freqStep=.005;
endFreq=0.1;

iter=floor((endFreq-startFreq)/freqStep);

signal=signal-mean(signal);
F=(startFreq:freqStep:endFreq);

for a=0:iter
    
    A=[];
    r2=[];
    %% Data Crunching  
    % Bandpass filtering to get data into frequency bins
    freq=(startFreq+a*freqStep);
    [bb,aa]= butter(3,[2*((a-1/2)*freqStep+startFreq)/sampf 2*((a+1/2)*freqStep+startFreq)/sampf],'bandpass');

    filtSignal=filter(bb,aa,signal);

    fitLength=floor(1/(freq/sampf)/4);

    for j=1:floor(length(filtSignal)/fitLength)-2

        tim=(j*fitLength:(j+1)*fitLength)'./sampf;
        cut=1e9*filtSignal(j*fitLength:(j+1)*fitLength);
        g = fittype( @(a,b,cen_fr,x) a*sin(2*pi*cen_fr*x)+b*cos(2*pi*cen_fr*x), 'problem', 'cen_fr' );
        
        fo = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',[0,0],...
               'Upper',[Inf,max(tim)],...
               'StartPoint',[1 1]);
           
        [myfit,st] = fit(tim, cut, g, 'problem', freq, fo);

        r2=[r2;st.rsquare];
        a=abs(coeffvalues(myfit));
        
        A=[A abs(a1(2))+abs(a1(1))*i];
        
    end
end

% 
% figure(8)
% hold on
% plot(F,abs(mag));
% fill([F, fliplr(F)], [err1, fliplr(err2)],'b','LineStyle','None');
% alpha(0.1)
% hold off
% set(gca,'YScale','log')
% set(gca,'XScale','log')
% 
% figure(9)
% X=[real(A);imag(A)]';
% hist3(X,[200 200],'CdataMode','auto','LineStyle','none')
% colorbar
% view(2)
% 
% figure(10)
% histogram(real(A),100)
% set(gca,'YScale','log')