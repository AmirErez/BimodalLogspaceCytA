
clear all
close all

'Check file number & column order & check ppERk gate'

NumberOfFiles=8;
NumberOfGates=30;
A=dir;
concentrations=logspace(-8,-8-NumberOfFiles+1,NumberOfFiles);

step=0.1;
x=-1:step:1;
y=x;
[XI,YI]=meshgrid(x,y);


for i=1:NumberOfFiles
    clear tempdata
    tempdata=[]
    for j=1:size(A,1)
        if findstr(['_00',num2str(i)], A(j).name)>0
        A(j).name
        Atemp=load(A(j).name);
        tempdata=[tempdata' Atemp']';
        end
    end
   
    
   length=size(tempdata,1);     
   tempdata2=tempdata(1+0*floor(length/3):3*floor(length/3),:);
    
   MeanData=mean(log10(tempdata2));
   
   data{i}(:,1)=log10(tempdata2(:,3))-MeanData(3);
   data{i}(:,2)=log10(tempdata2(:,4))-MeanData(4);
   data{i}(:,3)=tempdata2(:,6)>tempdata2(:,4)/1.5;
    
 
   SortedData{i}=sortrows(data{i},2);
  
   EventsPerGate=floor(size(data{i},1)/NumberOfGates);
    
    for j=1:NumberOfGates
            ppERK(i,j)=100*mean(SortedData{i}((j-1)*EventsPerGate+1:j*EventsPerGate,3),1);
            FoldChange(j)=mean(SortedData{i}((j-1)*EventsPerGate+1:j*EventsPerGate,2),1);
    end        
    
    [XI,YI,ZI{i}]=griddata(data{i}(:,2),data{i}(:,1),data{i}(:,3),XI,YI);
    
    subplot(3,3,i)
    image(x,y,ZI{i}*64)
    axis xy
    xlabel('ERK-1')
    ylabel('CD8\alpha')
    zlabel('ppERK (red=no response / blue = response)')
    title(['[OVA] = 10^{',num2str(-7-i),'} Mol']);
    
    
end

concentrations2=logspace(-8,-8-NumberOfFiles+1,21);


                   
figure
for j=1:NumberOfGates
    [estimates, model,sse] = Hill_fit(concentrations,ppERK(:,j)') 
       
    Amplitude2(j)=estimates(1);
    EC50_2(j)=estimates(2);
    Baseline2(j)=estimates(3);
    
    if estimates(2)>0
    hold on
    semilogx(concentrations,ppERK(:,j),'o','MarkerEdgeColor','k',...
                'MarkerFaceColor',[(j-1)/(NumberOfGates-1) 0 (NumberOfGates-j)/(NumberOfGates-1)]);
    hold off
   
    hold on
semilogx(concentrations2, Amplitude2(j)*concentrations2(:)./(concentrations2(:)+EC50_2(j))+Baseline2(j),'-k','LineWidth',1.3)
set(findobj(gca,'Type','line','Color',[0 0 0]),...
    'Color',[(j-1)/(NumberOfGates-1) 0 (NumberOfGates-j)/(NumberOfGates-1)])
    hold off
    end
    
end
set(gca,'XScale','log')
%title({'Differential response for different levels of ERK-1','blue=low levels - red=high levels','[072407]'})
%xlabel('[peptide] (Mol)')
ylabel('%ppERK^{+}','FontSize',24)

results=[(10.^FoldChange)' (Amplitude2)' Baseline2' EC50_2' (Amplitude2+Baseline2)']
save 'ERK1*ppERK_all.txt' results -ascii -tabs
save 'ppERK(diffERK1).txt' ppERK -ascii -tabs
                
for i=1:size(XI,1)
    for j=1:size(YI,1)
        for k=1:NumberOfFiles
            response(k)=ZI{k}(i,j);           
        end


        if sum(isnan(response))==0
            
            [estimates, model,sse] = Hill_fit2(concentrations,response);
            if estimates(2)>0 && estimates(2)<1e-3 %&& estimates(1)<2 && estimates(1)>0 
            
                Amplitude(i,j)=estimates(1)+estimates(3);
                EC50(i,j)=estimates(2);
                Background(i,j)=estimates(3);
                FitQuality(i,j)=sse;
            else
                EC50(i,j)=NaN;
            end
        else EC50(i,j)=NaN;
        end
    end
end

MidPoint=1-(14+log10(EC50))/7;
jet2=jet;
for i=1:5
    jet2(i,:)=[1 1 1];
end

figure
colormap(jet2);
image(10.^x,10.^y,MidPoint*64)
axis xy
xlabel('Normalized [ERK-1]')
ylabel('Normalized [CD8\alpha]')
hcb=colorbar('YTickLabel',{'>100nMol','10nMol','1nMol','100pMol','10pMol','<1pMol'})
set(hcb,'YTickMode','manual')

set(gca,'YTickLabel','manual')
set(gca,'ytick',[0 2.5 5.0 7.5 10])
set(gca,'YTickLabel',{'0.1','0.3','1','3','10'})

set(gca,'XTickLabel','manual')
set(gca,'xtick',[0 2.5 5.0 7.5 10])
set(gca,'XTickLabel',{'0.1','0.3','1','3','10'})


title({'EC_{50} for different levels of ERK-1 and CD8\alpha',...
    'Data [071807]'})

