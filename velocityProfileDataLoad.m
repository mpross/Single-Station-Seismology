%% Set-up
% close all

earthquakes=["Mexico_5_9" "Oklahoma_4_4" "Indonesia_6_9" "Fiji_8_2" "CostaRica_6_1" ...
    "Fiji_6_8" "Oregon_6_2" "Venezuela_7_3" "Peru_7_1" "Fiji_7_8" "NewZealand_6_9" "Canada_6_6" "Iceland_6_8"];
timeStamp=[1214366228 1212587999 1218725806 1218673195 1218583362 ...
    1218688157 1218965525 1218922324 1219136664 1220284172 1220588360 1224221998 1225763398];

vel=[];
vFreq=[];
sampF=8;

%% Data pull and decimate
j=length(earthquakes);
earthquakes(j);
filename=strcat('/home/michael/Google Drive/Seismology/Data/GPS',num2str(timeStamp(j)),'_',earthquakes(j),'.mat');

rawData=load(filename);

rawRY=rawData.rawData(1);
inRY=decimate(rawRY.data,rawRY.rate/8)*1e-9;

rawRX=rawData.rawData(2);
inRX=decimate(rawRX.data,rawRX.rate/8)*1e-9;

rawX=rawData.rawData(3); 
inX=decimate(rawX.data,rawX.rate/8)*1e-9;

rawY=rawData.rawData(4); 
inY=decimate(rawY.data,rawY.rate/8)*1e-9;

rawZ=rawData.rawData(5); 
inZ=decimate(rawZ.data,rawZ.rate/8)*1e-9;

inRY=inRY(1:1e5);
inRX=inRX(1:1e5);
inX=inX(1:1e5);
inY=inY(1:1e5);
inZ=inZ(1:1e5);

%% Time cut
timeThreshold=3e-6;

timeCut=find(abs(inZ-mean(inZ))>timeThreshold);
startTime=(timeCut(1)-1000);

timeCut=find(abs(fliplr(inZ')-mean(inZ))>timeThreshold);
endTime=(length(inZ)-(timeCut(1)-2000));
if(startTime<0)
    startTime=1;
end

if(endTime>length(inZ))
    endTime=length(inZ);
end
% startTime=1;
% endTime=length(inZ)/2;
%% Inversion and filtering

STSInvertFilt = zpk(-2*pi*[pairQ(8.33e-3,0.707)],-2*pi*[0 0],1);
STSInvertFilt = 1*STSInvertFilt/abs(freqresp(STSInvertFilt,2*pi*100));

time=(startTime:endTime)/sampF;

[b,a]=butter(3,0.02*2/sampF,'high');

RY=inRY(startTime:endTime);
RX=inRX(startTime:endTime);

X=lsim(STSInvertFilt,inX(startTime:endTime), time-startTime/sampF);
Y=lsim(STSInvertFilt,inY(startTime:endTime), time-startTime/sampF);
Z=lsim(STSInvertFilt,inZ(startTime:endTime), time-startTime/sampF);

X=filter(b,a,X);
Y=filter(b,a,Y);
Z=filter(b,a,Z);
RY=filter(b,a,RY);
RX=filter(b,a,RX);

X=X(250*sampF:end);
Y=Y(250*sampF:end);
Z=Z(250*sampF:end);
RY=RY(250*sampF:end);
RX=RX(250*sampF:end);
time=time(250*sampF:end);

% figure(1)
% plot1=plot((1:length(inZ))/8, inZ,(1:length(inZ))/8, (1:length(inZ))*0+2^15*1e-9, 'r',...
%     (1:length(inZ))/8, (1:length(inZ))*0-2^15*1e-9,'r')
% set(plot1,'LineWidth',1.5);
% set(gca,'FontSize',16);
% ylabel('Vertical Motion (m/s)')
% xlabel('Time (s)')
% legend('Seismometer','Limit')
% grid on
%     
% figure(2)
% spectrogram(Z,1024,512,1024,8,'yaxis')
% ax = gca;
% ax.YScale = 'log';
% ylim([1e-2 0.5])
% title('Z Velocity')
% 
% figure(3)
% spectrogram(RX,1024,512,1024,8,'yaxis')
% ax = gca;
% ax.YScale = 'log';
% ylim([1e-2 0.5])
% title('RX')
% 
% figure(4)
% spectrogram(RY,1024,512,1024,8,'yaxis')
% ax = gca;
% ax.YScale = 'log';
% ylim([1e-2 0.5])
% title('RY')
% 
% fig=figure(1)
% set(fig,'Units','Inches');
% pos = get(fig,'Position');
% set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print(fig,strcat('/home/michael/Google Drive/Seismology/EventPlots/',earthquakes(j),'_ClipCheck.pdf'),'-dpdf','-r1200')
% 
% 
% fig=figure(2)
% set(fig,'Units','Inches');
% pos = get(fig,'Position');
% set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print(fig,strcat('/home/michael/Google Drive/Seismology/EventPlots/', earthquakes(j), '_ZSpec.pdf'),'-dpdf','-r1200')
% 
% 
% fig=figure(3)
% set(fig,'Units','Inches');
% pos = get(fig,'Position');
% set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print(fig,strcat('/home/michael/Google Drive/Seismology/EventPlots/', earthquakes(j), '_RXSpec.pdf'),'-dpdf','-r1200')
% 
% 
% fig=figure(4)
% set(fig,'Units','Inches');
% pos = get(fig,'Position');
% set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print(fig,strcat('/home/michael/Google Drive/Seismology/EventPlots/', earthquakes(j), '_RYSpec.pdf'),'-dpdf','-r1200')