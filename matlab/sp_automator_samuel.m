   %Feature Extraction for a sliding window

%configuration for phone on table : 9.65, 10.00

filename='experiment_5_sample/crema_sad.csv';
fnameString="sad";
fnameNeumericStart=1;
     
num = csvread('experiment_5_sample/crema_sad.csv') ;
[r,c] = size(num) ;
disp("I am no one");

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

  chighIndices = find(calculate(:,1)>60);
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
largeCount=0;
     




start=10  ;
for x=1:37000               
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
    
    if min(az)<=1.007 || max(az)>=1.017
        %disp(start);
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
            %disp(tempEnd);
            %disp(tempEnd-appStart);
            %if(tempEnd-appStart>1.1) 
               %disp(tempEnd);
               %appEnd=tempEnd;
                %notInside=0; 
            %end
            
        elseif (notInside==1 && endObserver==1)
            endConsecutive=endConsecutive+1;
            if(endConsecutive>=2)
                if(tempEnd-appStart>2.5)
                    appEnd=tempEnd;
                    notInside=0;
                end
                
            end
            if(endConsecutive>=3)
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
    
    
    if(appStart~=0 && appEnd~=0 && appStart<appEnd && appEnd-appStart>=1)
        disp(appStart);
        disp(appEnd);
        disp("=================");
        regionCount=regionCount+1;
        savedStart=appStart;
        appStart=0;
        appEnd=0;
        
       
        
        fnameNeumericStart=fnameNeumericStart+1;
        convertedNeumeric=string(fnameNeumericStart);
        output_file=strcat(fnameString,convertedNeumeric,".png");
        draw_spectrogram(filename,savedStart,output_file)
    end
    
    
    

    %disp(meanX);
    %disp(maxX);
    %disp(minX);



    %disp(statX);



    meanCrossing=az > meanX;
    numberCrossing=sum(meanCrossing(:) == 1);
    meanCrossingRate=numberCrossing/numel(az);
    %disp(meanCrossingRate)

    %absoluteArea=sum(ax(ax<abs(ax)));
    %disp(meanabs(az));



    %Extracting frequency domain values:

    F = fft(az,1024);
    FFTCoEff=F/length(az);
    pow = F.*conj(F);
    total_pow = sum(pow);
    %disp(total_pow);





    Fs = 1/mean(diff(num(:,1)));
    ktest=diff((num(:,1)));
    %disp(Fs);

    pen=pentropy(az ,Fs);
    sumPen=sum(pen);
    %disp(start);
    %disp(final);
    %disp(sumPen);
    
    


    %coeffs=fftcoeffs(ax);

    hd = dfilt.fftfir(az,1024);
    c=fftcoeffs(hd);
    amp = 2*abs(c)/length(az);
    phase=angle(c);
    magnitude=abs(amp);

    highestMag=max(magnitude);
    sumMag=sum(magnitude);

    frequency_ratio=highestMag/sumMag;
    %disp(start);
    %disp(final);

    %disp(frequency_ratio);
    
 
    decision=[start final min(az) max(az)];
    %disp(decision);

    %statX=[mean(az) max(az) min(az) std(az) var(az) range(az) (std(az)/mean(az))*100 skewness(az) kurtosis(az) quantile(az,[0.25,0.50,0.75]) meanCrossingRate total_pow sumPen frequency_ratio];
    %disp(min(az))
    %t=array2table(statX);
    %xlswrite('experiment.xls',statX);
    %xlsappend('experiment.xls',statX);
    %xlsappend( 'experiment.xlsx', { 'test1' , 'test two' } )

    %xlswrite('experiment.xls',statX,'Sheet1','A7');
    %xlswrite('experiment.xls',statX,'Sheet1','A8');
    %xlsappend('experiment.xls',statX);
end

fprintf('Total word region found %d\n',regionCount);
fprintf('Total large region found %d\n',largeCount);



function draw_spectrogram(filename,start_point,output_file)
num = csvread(filename) ;
[r,c] = size(num) ;
%Delete unecessary information    
num(2:2:end,:) = []    ;

end_point=start_point+3;
%15.8-18.8
lowIndices = find(num(:,1)<start_point);
num(lowIndices,:) = [];

highIndices = find(num(:,1)>end_point);
num(highIndices,:) = [];

%Extract all axes
ax = num(:,2) ;
ay = num(:,3) ;
az = num(:,4) ;
%Compute the cumulative/magnitude vector
acc = sqrt(ax.^2 + ay.^2 + az.^2); %possibly in G's
acc_cmpersec = acc.*980; %converts to cm/s^2

num(:,5)= acc_cmpersec;
Fs = 1/mean(diff(num(:,1)));  


figure;
%myMap = [linspace(0,1,256)', zeros(256,2)];
var0=[zeros(170,3)];
var1 = [zeros(30,2),linspace(0,1,30)'];
var2 = [linspace(0,1,56)',linspace(0,1,56)',ones(56,1)];
%hello=ones(256,2);
myMap=vertcat(var0,var1,var2);

%disp(size(myMap))



colormap(myMap)
%colormap(colorcube)
spectrogram(num(:,5),128,120,128,Fs,'yaxis');
%imwrite(imresize(s,[227,227]),"filename.png")
set(gca, 'Visible', 'off')
set(colorbar, 'Visible', 'off')
set(groot,'defaultFigureVisible','off')
saveas(gcf,output_file);
end








