props = java.lang.System.getProperties;
props.setProperty('mail.smtp.port', '587');
props.setProperty('mail.smtp.auth','true');
props.setProperty('mail.smtp.starttls.enable','true');

setpref('Internet','E_mail','mprossmatlab@gmail.com');
setpref('Internet','SMTP_Server','smtp.gmail.com');
setpref('Internet','SMTP_Username','mprossmatlab');
setpref('Internet','SMTP_Password','matlab123');

setPar=[1.5 0.5 1 3 1 2 3.5 2.25]*1e3;
t0=cputime;
testFreq=(5:35)*1e-2;

%First layer
vP1=rand*1e1+setPar(1);
vS1=rand*1e1+setPar(2); %must be between 0 and /srt(2) Landau
d1=rand*1e1+setPar(3);

%Second layer
vP2=rand*1e1+setPar(4);
vS2=rand*1e1+setPar(5);
d2=rand*1e1+setPar(6);

%Third layer
vP3=rand*1e1+setPar(7);
vS3=rand*1e1+setPar(8);

testDispers=dispersionCalc(setPar(1),setPar(2),setPar(3),setPar(4),setPar(5),setPar(6),setPar(7),setPar(8),testFreq,3);

[bestPar,bestDispers,frames]=dispersionFit(testFreq,testDispers,3);

fig1=figure(4);
plot2=plot(testFreq,testDispers,'.',testFreq,bestDispers);
ylabel('Velocity (m/s)')
xlabel('Frequency (Hz)')
set(plot2,'LineWidth',1.5);
set(gca,'FontSize',16);
set(plot2,'MarkerSize',16);

t=(cputime-t0)/3600;
disp(["Completion Time: "+num2str(t)+" hours"+newline+...
   newline+"Fit Parameters: "+num2str(bestPar,4)+newline+"Set Parameters: "+num2str(setPar)+...
   newline+"Difference: "+num2str(bestPar-setPar)])
print(fig1,'-dpng','Test_Dispersion.png');
sendmail('mpross2@uw.edu','Fit Test Complete',"Completion Time: "+num2str(t)+" hours"+newline+...
   newline+"Fit Parameters: "+num2str(bestPar)+newline+"Set Parameters: "+num2str(setPar),{'Test_Dispersion.png','dispersionSwarm.avi'});
