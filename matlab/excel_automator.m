%Feature Extraction for a sliding window

%configuration for phone on table : 9.65, 10.00

filename='experiment_5_sample/sad.csv';
fnameString="test";
fnameNeumericStart=1;

num = csvread(filename) ;
[r,c] = size(num) ;

%Delete unecessary information
num(2:2:end,:) = [] ;


regionCount=0;
consecutiveChecker=0;
consecutiveStart=0;
consecutiveEnd=0;

%Fs = 1/mean(diff(num(:,1)));  
%y_highpass=highpass(num(:,4),20,Fs);
%num(:,4)=y_highpass;


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





start=10;
for x=1:12000 
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
    
     %disp(x);
    %disp(minX)
    %disp(maxX)
    
    
    
    if min(az)<1.000 || max(az)>1.02
       
        if(notInside==0 && startObserver==0)
            startObserver=1;
            tempStart=start;
            startConsecutive=1;  
           
            
        elseif(notInside==0 && startObserver==1)
                startConsecutive=startConsecutive+1;
                if(startConsecutive>=3)
                    notInside=1;
                    startObserver=0;
                    appStart=tempStart;
                     
                    
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
           
            
        elseif (notInside==1 && endObserver==1)
            endConsecutive=endConsecutive+1;
            if(endConsecutive>=3)
                appEnd=tempEnd;
                notInside=0;
                tempEnd=0;
            end
 
        else
            if(notInside==0 && startObserver==1)
                startObserver=0;
                tempStart=0;
            end
        end
        
    end 
    
    
    %Now printing start and end
    
    if(appStart~=0 && appEnd~=0 && appStart<appEnd && appEnd-appStart>=1)
        disp(appStart);
        disp(appEnd);
        disp("==========");
        regionCount=regionCount+1;
        savedStart=appStart;
        savedEnd=appEnd;
        appStart=0;
        appEnd=0;
        
       
        
 
        my_xls(filename,savedStart,savedEnd);
    end
    
end
fprintf('Total word region found %d\n',regionCount);

function my_xls(filename,start,funcend)


class="sad";

startValue=start;
endValue=funcend;
%watchCompareStart=20;
%watchCompareEnd=34;
phoneCompareStart=10;
phoneCompareEnd=80;
%disp(watchStartValue);
%disp(endValue);
%disp(watchEndValue);

%now starting the smartphonw calculations

numPhone = csvread(filename) ;
[r1,c1] = size(numPhone) ;

numPhone(2:2:end,:) = [] ;

Fsp = 1/mean(diff(numPhone(:,1)));  
y_highpass=highpass(numPhone(:,4),20,Fsp);
numPhone(:,4)=y_highpass;

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

phoneRms = numPhone(:,5) ;

phoneZ=numPhone(:,4);

comparePZ=comparePhone(:,4);




meanP=mean(phoneZ);
minP=min(phoneZ);
maxP=max(phoneZ);
meanPZ=mean(comparePZ);

%disp(meanP);
%disp(minP);
%disp(maxP);

meanCrossingP=phoneZ > meanPZ;
numberCrossingP=sum(meanCrossingP(:) == 1);
meanCrossingRateP=numberCrossingP/numel(phoneZ);



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





hdp = dfilt.fftfir(phoneZ,1024);
cp=fftcoeffs(hdp);
ampp = 2*abs(cp)/length(phoneZ);
phasep=angle(cp);
magnitudep=abs(ampp);

highestMagp=max(magnitudep);
sumMagp=sum(magnitudep);

frequency_ratiop=highestMagp/sumMagp;



statX=[mean(phoneZ) max(phoneZ) min(phoneZ) std(phoneZ) var(phoneZ) range(phoneZ) (std(phoneZ)/mean(phoneZ))*100 skewness(phoneZ) kurtosis(phoneZ) quantile(phoneZ,[0.25,0.50,0.75]) meanCrossingRateP total_powp sumPenp frequency_ratiop class];

t=array2table(statX);

[numbers, strings, raw] = xlsread('experimentphone.xls');
lastRow = size(raw, 1);
nextRow = lastRow + 1;
cellReference = sprintf('A%d', nextRow);
xlswrite('experimentphone.xls', statX, 'Sheet1', cellReference);










end


