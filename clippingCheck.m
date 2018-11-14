%% Set-up
% close all

earthquakes=["Mexico_5_9" "Oklahoma_4_4" "Indonesia_6_9" "Fiji_8_2" "CostaRica_6_1" ...
    "Fiji_6_8" "Oregon_6_2" "Venezuela_7_3" "Peru_7_1" "Fiji_7_8" "NewZealand_6_9" "Canada_6_6" "Iceland_6_8"];
timeStamp=[1214366228 1212587999 1218725806 1218673195 1218583362 ...
    1218688157 1218965525 1218922324 1219136664 1220284172 1220588360 1224221998 1225763398];

vel=[];
vFreq=[];
sampF=8;

pass=[];

%% Data pull and decimate
for j=(1:length(earthquakes))
    earthquakes(j)
    filename=strcat('/home/michael/Google Drive/Seismology/Data/GPS',num2str(timeStamp(j)),'_',earthquakes(j),'.mat');

    rawData=load(filename);

    rawBRSX=rawData.rawData(1);
    inBRSX=decimate(rawBRSX.data,rawBRSX.rate/8)*1e-9;

    rawBRSY=rawData.rawData(2);
    inBRSY=decimate(rawBRSY.data,rawBRSY.rate/8)*1e-9;

    rawSTSX=rawData.rawData(3); 
    inSTSX=decimate(rawSTSX.data,rawSTSX.rate/8)*1e-9;

    rawSTSY=rawData.rawData(4); 
    inSTSY=decimate(rawSTSY.data,rawSTSY.rate/8)*1e-9;

    rawSTSZ=rawData.rawData(5); 
    inSTSZ=decimate(rawSTSZ.data,rawSTSZ.rate/8)*1e-9;

    inBRSX=inBRSX(1:1e5);
    inBRSY=inBRSY(1:1e5);
    inSTSX=inSTSX(1:1e5);
    inSTSY=inSTSY(1:1e5);
    inSTSZ=inSTSZ(1:1e5);

    %% Time cut
    timeThreshold=3e-6;
    % 
    % timeCut=find(abs(inSTSZ-mean(inSTSZ))>timeThreshold);
    % startTime=(timeCut(1)-1000);
    % 
    % timeCut=find(abs(fliplr(inSTSZ')-mean(inSTSZ))>timeThreshold);
    % endTime=(length(inSTSZ)-(timeCut(1)-2000));
    % if(startTime<0)
    %     startTime=1;
    % end
    % 
    % if(endTime>length(inSTSZ))
    %     endTime=length(inSTSZ);
    % end
    startTime=1;
    endTime=length(inSTSZ)/2;
    %% Inversion and filtering

    STSInvertFilt = zpk(-2*pi*[pairQ(8.33e-3,0.707)],-2*pi*[0 0],1);
    STSInvertFilt = 1*STSInvertFilt/abs(freqresp(STSInvertFilt,2*pi*100));

    time=(startTime:endTime)/sampF;

    [b,a]=butter(3,0.02*2/sampF,'high');

    BRSX=inBRSX(startTime:endTime);
    BRSY=inBRSY(startTime:endTime);

    STSX=lsim(STSInvertFilt,inSTSX(startTime:endTime), time-startTime/sampF);
    STSY=lsim(STSInvertFilt,inSTSY(startTime:endTime), time-startTime/sampF);
    STSZ=lsim(STSInvertFilt,inSTSZ(startTime:endTime), time-startTime/sampF);

    STSX=filter(b,a,STSX);
    STSY=filter(b,a,STSY);
    STSZ=filter(b,a,STSZ);
    BRSX=filter(b,a,BRSX);
    BRSY=filter(b,a,BRSY);

    STSX=STSX(250*sampF:end);
    STSY=STSY(250*sampF:end);
    STSZ=STSZ(250*sampF:end);
    BRSX=BRSX(250*sampF:end);
    BRSY=BRSY(250*sampF:end);
    time=time(250*sampF:end);
    
    pass=[pass; max(abs(inSTSZ))<=3.333e-4];

    figure(1)
    plot1=plot((1:length(inSTSZ))/8, inSTSZ,(1:length(inSTSZ))/8, (1:length(inSTSZ))*0+3.333e-4, 'r',...
        (1:length(inSTSZ))/8, (1:length(inSTSZ))*0-3.333e-4,'r')
    set(plot1,'LineWidth',1.5);
    set(gca,'FontSize',16);
    ylabel('Vertical Motion (m/s)')
    xlabel('Time (s)')
    legend('Seismometer','Limit')
    grid on

    figure(2)
    spectrogram(STSZ,1024,512,1024,8,'yaxis')
    ax = gca;
    ax.YScale = 'log';
    ylim([1e-2 0.5])
    title('Z Velocity')

    figure(3)
    spectrogram(BRSY,1024,512,1024,8,'yaxis')
    ax = gca;
    ax.YScale = 'log';
    ylim([1e-2 0.5])
    title('RX')

    figure(4)
    spectrogram(BRSX,1024,512,1024,8,'yaxis')
    ax = gca;
    ax.YScale = 'log';
    ylim([1e-2 0.5])
    title('RY')

    fig=figure(1)
    set(fig,'Units','Inches');
    pos = get(fig,'Position');
    set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(fig,strcat('/home/michael/Google Drive/Seismology/EventPlots/',earthquakes(j),'_ClipCheck.jpg'), '-djpeg','-r0')


    fig=figure(2)
    set(fig,'Units','Inches');
    pos = get(fig,'Position');
    set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(fig,strcat('/home/michael/Google Drive/Seismology/EventPlots/', earthquakes(j), '_ZSpec.jpg'), '-djpeg','-r0')


    fig=figure(3)
    set(fig,'Units','Inches');
    pos = get(fig,'Position');
    set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(fig,strcat('/home/michael/Google Drive/Seismology/EventPlots/', earthquakes(j), '_RXSpec.jpg'), '-djpeg','-r0')


    fig=figure(4)
    set(fig,'Units','Inches');
    pos = get(fig,'Position');
    set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(fig, strcat('/home/michael/Google Drive/Seismology/EventPlots/', earthquakes(j), '_RYSpec.jpg'), '-djpeg','-r0')
end
pass