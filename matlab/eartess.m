%Feature Extraction for a sliding win
%configuration for phone on table : 9.65, 10
theVal=0; 

filename='experiment_5_sample/tessfullset.csv';
fnameString="test"; 
fnameNeumericStart=1; 

num = csvread(filename) ;
num=num(2293000:end,:);
[r,c] = size(num) ;
timeValue=num(:,1); 
startVal=timeValue(1);
endVal=timeValue(end);
theval=endVal-0;
disp(theVal)
endVal=(endVal-10)*10;
disp(startVal);
disp(endVal);

%Delete unecessary information
num(2:2:end,:) = [] ;


regionCount=0;
consecutiveChecker=0;
consecutiveStart=0;
consecutiveEnd=0;
savedEnd=0;

%Fs = 1/mean(diff(num(:,1)));  
%y_highpass=highpass(num(:,4),20,Fs);
%num(:,4)=y_highpass;
%high pass filter
Fs = 1/mean(diff(num(:,1)));  
y_highpass=highpass(num(:,4),18,Fs);
num(:,4)=y_highpass;


calculate=num;

 %Delete rows for specific condition
  clowIndices = find(calculate(:,1)<10);
  calculate(clowIndices,:) = []; 

  chighIndices = find(calculate(:,1)>80);
  calculate(chighIndices,:) = [];
  
  mainZ=calculate(:,4) ;
  meanV=mean(mainZ);
  
  %disp(mainZ)



notInside=0;
startObserver=0;
startConsecutive=0;
appStart=0;
endConsecutive=0;
tempEnd=0;
appEnd=0;
endObserver=0;
regionCount=0;
errorCount=0;
errorState=0;  
largeCount=0;
smallCount=0;





start=5447;
for x=1:10000 
    iterate=num;
    %disp(iterate);
    start=start+0.1;
    initValue=start;
    final=start+0.1;

    %Delete rows for specific condition
    lowIndices = find(iterate(:,1)<initValue);
    iterate(lowIndices,:) = [];

    highIndices = find(iterate(:,1)>final);
     iterate(highIndices,:) = [];

    %disp(lowIndices);

    



    %Extract all axes
    ax = iterate(:,2) ;
    ay = iterate(:,3) ;
    az = iterate(:,4) ;
    %disp(az);

    meanX=mean(az);
    %disp(meanX);
    minX=min(az);
    maxX=max(az);
    
    %disp(start);
    %disp(min(az));
    %disp(max(az));
    
   
    
 if min(az)<=-0.005 || max(az)>=0.003
        %disp(start);
        if(notInside==0 && startObserver==0)
            startObserver=1;
            tempStart=start;
            startConsecutive=1;        
        elseif(notInside==0 && startObserver==1)
                startConsecutive=startConsecutive+1;
                if(startConsecutive>=2)
                    notInside=1;
                    startObserver=0;
                    appStart=tempStart;
                    %fprintf('The starting Point is %d\n',appStart);
                    if(appStart-savedEnd)>1.5 && (appStart-errorState)>2.0
                        errorCount=errorCount+1;
                        %disp(appStart);
                        errorState=appStart;
                        fprintf('Possible Error before that %d\n',appStart);
                        
                    end
                end
        elseif(notInside==1 && endObserver==1)
            endObserver=0;
            tempEnd=0;
        end
         
            
    else
        if (notInside==1 && endObserver==0)
            endConsecutive=1;
            endObserver=1;
            tempEnd=final;
            %disp(tempEnd);
            %disp(tempEnd-appStart);
            %if(tempEnd-appStart>1.1) 
               %disp(tempEnd);
               %appEnd=tempEnd;
                %notInside=0; 
            %end
            
        elseif (notInside==1 && endObserver==1)
            endConsecutive=endConsecutive+1;
            if(endConsecutive>=1)
                if(tempEnd-appStart>1.5)
                    appEnd=tempEnd;
                    notInside=0;
                end
                
            end
            if(endConsecutive>=4)
                appEnd=tempEnd;
                notInside=0;              
            end
      
        else
            if(notInside==0 && startObserver==1)
                startObserver=0;
                tempStart=0;
            end
        end  
        
    end 
    
    
    
    %Now printing start and end
    
    if(appStart~=0 && appEnd~=0 && appStart<appEnd && appEnd-appStart>=0.8 && appEnd-appStart<=2.8)
        disp(appStart);
        disp(appEnd);
        disp("==========");
        
        if(appEnd-appStart)>=2.3
            largeCount=largeCount+1;
        elseif(appEnd-appStart)<=1.0
            smallCount=smallCount+1;
        end
        regionCount=regionCount+1;
        savedStart=appStart;
        savedEnd=appEnd;
        appStart=0;
        appEnd=0;
        
        
        class="sad";
        
 
        %my_xls(filename,savedStart,savedEnd,endVal,class);
    end
    
end
fprintf('Total word region found %d\n',regionCount);
fprintf('Total possible error found %d\n',errorCount);
fprintf('Total oversized word regions %d\n',largeCount);
fprintf('Total undersized word regions %d\n',smallCount);

function my_xls(filename,start,funcend,endVal,class)


startValue=start;
endValue=funcend;
%watchCompareStart=20;
%watchCompareEnd=34;
phoneCompareStart=10;
phoneCompareEnd=endVal+10;
%disp(watchStartValue);
%disp(endValue);
%disp(watchEndValue);

%now starting the smartphonw calculations

numPhone = csvread(filename) ;
[r1,c1] = size(numPhone) ;

numPhone(2:2:end,:) = [] ;



comparePhone=numPhone;

%Delete rows for phone compare values
comparePLow = find(comparePhone(:,1)<phoneCompareStart);
comparePhone(comparePLow,:) = [];

comparePHigh = find(comparePhone(:,1)>phoneCompareEnd);
comparePhone(comparePHigh,:) = [];




%Delete rows for specific condition in Phone
lowIndicesp = find(numPhone(:,1)<startValue);
numPhone(lowIndicesp,:) = [];

highIndicesp = find(numPhone(:,1)>endValue);
numPhone(highIndicesp,:) = [];

Fsp = 1/mean(diff(numPhone(:,1)));  
Fn=Fsp/2;
%y_highpass=highpass(numPhone(:,4),20,Fsp);
%numPhone(:,4)=y_highpass;

phoneRms = numPhone(:,5) ;

phoneZ=numPhone(:,4);
comparePZ=comparePhone(:,4);

%Calculating Frequency domain features

Tr = linspace(numPhone(1,1), numPhone(1,end), size(numPhone,1));  
Dr = resample(phoneZ, Tr); 
Dr_mc  = Dr - mean(Dr,1); 


FDr_mc = fft(Dr_mc, [], 1);
Fv = linspace(0, 1, fix(size(FDr_mc,1)/2)+1)*Fn; 

Iv = 1:numel(Fv); 
amplitude=abs(FDr_mc(Iv,:))*2;

upperPart=Fv*amplitude;
ampSum=sum(amplitude);

specCentroid=upperPart/ampSum;
%disp(specCentroid); 

FvSqr=Fv.^2;
stdDevupper=FvSqr*amplitude;
specStdDev=sqrt(stdDevupper/ampSum);
specCrest=max(amplitude)/specCentroid;


specSkewness=(((Fv-specCentroid).^3)*amplitude)/(specStdDev)^3;

specKurt=(sum((((amplitude-specCentroid).^4).*amplitude))/(specStdDev)^4)-3 ;
maxFreq=max(Fv);
maxMagx=max(phoneZ);







meanP=mean(phoneZ);
minP=min(phoneZ);
maxP=max(phoneZ);
meanPZ=mean(comparePZ);
%gradientZ=mean(gradient(phoneZ));
%disp(meanPZ); 
irrk=irregularityk(phoneZ);
irrj=irregularityj(phoneZ);
sharp=sharpness(phoneZ);
smooth=smoothness(phoneZ);

%now adding frequency domain things:





%disp(meanP);
%disp(minP);
%disp(maxP);

meanCrossingP=phoneZ > meanPZ;
numberCrossingP=sum(meanCrossingP(:) == 1);
meanCrossingRateP=numberCrossingP/numel(phoneZ);
%disp(meanCrossingRateP);



%Extracting frequency domain values:

Fp = fft(phoneZ,1024);
FFTCoEffp=Fp/length(phoneZ);
powp = Fp.*conj(Fp);
total_powp = sum(powp);
%disp(total_powp);



Fsp = 1/mean(diff(numPhone(:,1)));

penp=pentropy(phoneZ,Fsp);
sumPenp=sum(penp);

%disp(sumPenp);

%centroid=spectralCentroid(phoneZ,Fsp);
%disp(centroid)

%sharpness = acousticSharpness(phoneZ,Fsp);
%disp(sharpness);



hdp = dfilt.fftfir(phoneZ,1024);
cp=fftcoeffs(hdp);
ampp = 2*abs(cp)/length(phoneZ);
phasep=angle(cp);
magnitudep=abs(ampp);

highestMagp=max(magnitudep);
sumMagp=sum(magnitudep);

frequency_ratiop=highestMagp/sumMagp;



statX=[mean(phoneZ) max(phoneZ) min(phoneZ) std(phoneZ) var(phoneZ) range(phoneZ) (std(phoneZ)/mean(phoneZ))*100 skewness(phoneZ) kurtosis(phoneZ) quantile(phoneZ,[0.25,0.50,0.75]) meanCrossingRateP total_powp sumPenp frequency_ratiop irrk irrj sharp smooth specCentroid specStdDev specCrest specSkewness specKurt maxFreq maxMagx class];

t=array2table(statX);

[numbers, strings, raw] = xlsread('experimentphone.xls');
lastRow = size(raw, 1);
nextRow = lastRow + 1;
cellReference = sprintf('A%d', nextRow);
xlswrite('experimentphone.xls', statX, 'Sheet1', cellReference);










end


function ikSum=irregularityk(phonez)
%disp(phonez);
N=[10 100 1000];

ikSum=0;
for i=1:length(phonez)-2
   ik=(phonez(i+1)-(phonez(i)+phonez(i+1)+phonez(i+2)/3));
   ikSum=ikSum+ik;
    
end
end

function ijSum=irregularityj(phonez)


ijSum=0;
for i=1:length(phonez)-1
   ij1=(phonez(i)-phonez(i+1))^2;
   ij2=phonez(i)^2;
   ij=ij1/ij2;
   ijSum=ijSum+ij;
    
end
end

function finalsharp=sharpness(phonez)
sharpn=0;
tempi=0;
for i=1:length(phonez)
    if(i<15)
        tempi=real(i*phonez(i)^0.23);
        %disp(tempi);
    else
        tempi=real(0.066*exp(0.171*i)*i*phonez(i)^0.23);
    end
    
    sharpn=sharpn+tempi;
end

finalsharp=(0.11*sharpn)/length(phonez);
end

function smoothSum=smoothness(phonez)

smoothSum=0;
for i=1:length(phonez)-2
   
    ismooth=real((20*log(phonez(i))-(20*log(phonez(i))+20*log(phonez(i+1))+20*log(phonez(i+2))))/3);
    
    smoothSum=smoothSum+ismooth;
end


end








