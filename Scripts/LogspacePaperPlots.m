% Makes plots for logspace paper

%% Plot histograms in log space and linear space 

sig1s=0.4:0.4:1.2;
sig2s=0.4:0.4:1.2;

p.alpha = 0.5;
p.y1 = log(100);
p.y2 = log(1000);

P = @(y,p)( (1-p.alpha)/(p.sig1*sqrt(2*pi))*exp(-((y-p.y1).^2)/(2*p.sig1^2)) + p.alpha/(p.sig2*sqrt(2*pi))*exp(-((y-p.y2).^2)/(2*p.sig2^2)));
Q = @(y,p)(exp(-y).*P(y,p));
figure;
for s1=1:length(sig1s)
    for s2=1:length(sig2s)
        yl = [10^-5, 10^5];
        yrange = log([10^0,10^5]);
        ys = linspace(yrange(1),yrange(2),1001);
        p.sig1 = sig1s(s1);
        p.sig2 = sig2s(s2);
        
        subplot(3,3,s2+(s1-1)*length(sig2s)); 
        valQ = Q(ys,p);
        valQ = valQ / max(valQ(~isinf(valQ)));
        semilogx(exp(ys),valQ,'-b','LineWidth',2);
        hold on
        set(gca,'FontSize',14);   
        valP = P(ys,p);
        valP = valP / max(valP(~isinf(valP)));
        semilogx(exp(ys),valP,'-r','LineWidth',2);
        xlim(exp(yrange));
        title(['$\sigma_1=' num2str(sig1s(s1)) ', \sigma_2=' num2str(sig2s(s2)) '$'],'Interpreter','Latex');
        set(gca,'XTick',10.^(1:2:5));
        set(gca,'YTick',[]);
    end
end

figname = 'PlotsLinearLogTogether';
set(gcf,'Name',figname);
% print(gcf,'-dtiff',['Plots' filesep figname],'-r600');
% print(gcf,'-depsc2',['Plots' filesep figname]);
close(gcf);

%% Plot histograms in log space and linear space separately

sig1s=0.4:0.4:1.2;
sig2s=0.4:0.4:1.2;

p=struct;
p.alpha = 0.5;
p.y1 = log(100);
p.y2 = log(1000);

P = @(y,p)( (1-p.alpha)/(p.sig1*sqrt(2*pi))*exp(-((y-p.y1).^2)/(2*p.sig1^2)) + p.alpha/(p.sig2*sqrt(2*pi))*exp(-((y-p.y2).^2)/(2*p.sig2^2)));
Q = @(y,p)(exp(-y).*P(y,p));

% Log first
figure;
for ss=1:length(sig1s)
    s1 = ss;
    s2 = ss;
    yl = [10^-5, 10^5];
    yrange = log([10^0,10^5]);
    ys = linspace(yrange(1),yrange(2),1001);
    p.sig1 = sig1s(s1);
    p.sig2 = sig2s(s2);
        
    subplot(3,4,ss+1);       
    valP = P(ys,p);
%     valP = valP / trapz(ys,valP);
    valP = valP / max(valP(~isinf(valP)));
    semilogx(exp(ys),valP,'-r','LineWidth',1);
    hold on
    
    set(gca,'FontSize',10);
    xlim(exp(yrange));
    title(['$\sigma_1=\sigma_2=' num2str(sig2s(s2)) '$'],'Interpreter','Latex');
    
    % Draw 10000 "cells" from the distribution and run Hartigan's dip test
%     draw1 = randn(5000,1)*p.sig1+p.y1;
%     draw2 = randn(5000,1)*p.sig2+p.y2;
%     [dip, pval] = HartigansDipSignifTest([draw1;draw2], 500);
%     t = text(1500,0.9,['$p_U=' num2str(pval,3) '$'],'Interpreter','Latex','FontSize',10);
    
    set(gca,'XTick',10.^(1:2:5));
    set(gca,'YTick',[]);
    semilogx([100 100],[0 1],'--k');
    semilogx([1000 1000],[0 1],'--k');
    set(gca,'LineWidth',1);
    set(gca,'GridAlpha',1);
end
subplot(3,4,1);
semilogx(exp(yrange),NaN*ones(2,1));
set(gca,'XTick',[]);
set(gca,'YTick',[]);
xlim(exp(yrange));
set(gca,'FontSize',10);
set(gca,'Box','off');
xlabel('$\log I$','Interpreter','Latex','Color','r');
ylabel('$P(\log I)$','Interpreter','Latex','Color','r');

% Linear next
% figure;
for ss=1:length(sig1s)
    s1 = ss;
    s2 = ss;
    yl = [10^-5, 10^5];
    yrange = log([10^0,10^5]);
    ys = linspace(yrange(1),yrange(2),1001);
    p.sig1 = sig1s(s1);
    p.sig2 = sig2s(s2);
    
    subplot(3,4,4+ss+1);       
    valP = P(ys,p);
    valQ = Q(ys,p);
    valQ = valQ / max(valQ(~isinf(valQ)));
    plot(exp(ys),valQ,'-m','LineWidth',1);
    hold on
    set(gca,'FontSize',10);
%     xlim(exp(yrange));
    xlim([0 2000]);
    title(['$\sigma_1=\sigma_2=' num2str(sig2s(s2)) '$'],'Interpreter','Latex');
    set(gca,'XTick',[0 1000 2000]);
    set(gca,'YTick',[]);
    plot([100 100],[0 1],'--k');
    plot([1000 1000],[0 1],'--k');
        
%     draw1 = lognrnd(p.y1,p.sig1,5000,1);
%     draw2 = lognrnd(p.y2,p.sig2,5000,1);
%     [dip, pval] = HartigansDipSignifTest([draw1;draw2], 500);
%     t = text(1100,0.9,['$p_U=' num2str(pval,3) '$'],'Interpreter','Latex','FontSize',10);
    
    
%     grid('on');
    set(gca,'LineWidth',1);
    set(gca,'GridAlpha',1);
end
subplot(3,4,1+4);
semilogx(exp(yrange),NaN*ones(2,1));
set(gca,'XTick',[]);
set(gca,'YTick',[]);
xlim(exp(yrange));
set(gca,'FontSize',10);
set(gca,'Box','off');
xlabel('$I$','Interpreter','Latex','Color','m');
ylabel('$Q(I)$','Interpreter','Latex','Color','m');


% Linear on semilogx plot
% Linear next
% figure;
for ss=1:length(sig1s)
    s1 = ss;
    s2 = ss;
    yl = [10^-5, 10^5];
    yrange = log([10^0,10^5]);
    ys = linspace(yrange(1),yrange(2),1001);
    p.sig1 = sig1s(s1);
    p.sig2 = sig2s(s2);
    
    subplot(3,4,8+ss+1);       
    valP = P(ys,p);
    valQ = Q(ys,p);
    valQ = valQ / max(valQ(~isinf(valQ)));
    semilogx(exp(ys),valQ,'-b','LineWidth',1);
    hold on
    set(gca,'FontSize',10);
%     xlim(exp(yrange));
    xlim(exp(yrange));
    title(['$\sigma_1=\sigma_2=' num2str(sig2s(s2)) '$'],'Interpreter','Latex');

  % Draw 10000 "cells" from the distribution and run Hartigan's dip test
%     draw1 = lognrnd(p.y1,p.sig1,5000,1);
%     draw2 = lognrnd(p.y2,p.sig2,5000,1);
%     [dip, pval] = HartigansDipSignifTest([draw1;draw2], 500);
%     t = text(1500,0.9,['$p_U=' num2str(pval,3) '$'],'Interpreter','Latex','FontSize',10);
    
%     title(['p=' num2str(pval,3)]);

    set(gca,'XTick',10.^(1:2:5));
    set(gca,'YTick',[]);
    semilogx([100 100],[0 1],'--k');
    semilogx([1000 1000],[0 1],'--k');
    set(gca,'LineWidth',1);
    set(gca,'GridAlpha',1);
end
subplot(3,4,1+8);
semilogx(exp(yrange),NaN*ones(2,1));
set(gca,'XTick',[]);
set(gca,'YTick',[]);
xlim(exp(yrange));
set(gca,'FontSize',10);
set(gca,'Box','off');
xlabel('$\log I$','Interpreter','Latex','Color','b');
ylabel('$Q(I)$','Interpreter','Latex','Color','b');

set(gcf,'Color','w');
set(gcf,'InvertHardcopy','off');
set(gcf,'Units','Inches');
set(gcf,'Position',[3,4,5.2,3]);
set(gcf,'PaperUnits','Inches');
set(gcf,'PaperPosition',get(gcf,'Position'));

figname = 'PlotsLinearLogTogether';
% figname = 'Fig1';
set(gcf,'Name',figname);
print(gcf,'-dtiff',['Plots' filesep figname],'-r600');
print(gcf,'-dpng',['Plots' filesep figname],'-r600');




%% Derivative
% F = @(y,p)(-p.alpha/(p.alpha-1)*(p.sig1^3)/(p.sig2^3)*exp(((y-p.y1).^2)/(2*p.sig1^2)-((y-p.y2).^2)/(2*p.sig2^2)));
% S1 = @(y,p)(-(p.sig1^2+y-p.y1)./(p.sig2^2+y-p.y2));
% S3 = @(y,p)(-(y-p.y1)./(y-p.y2));

F = @(y,p)(p.alpha/(1-p.alpha)*(p.sig1)/(p.sig2)*exp(((y-p.y1).^2)/(2*p.sig1^2)-((y-p.y2).^2)/(2*p.sig2^2)));
S1 = @(y,p)(((y-p.y1)/p.sig1^2+1)./((p.y2-y)/p.sig2^2-1));
S3 = @(y,p)(((y-p.y1)/p.sig1^2)./((p.y2-y)/p.sig2^2));


figure;
for s1=1:length(sig1s)
    for s2=1:length(sig2s)
        p.sig1 = sig1s(s1);
        p.sig2 = sig2s(s2);
        
        subplot(3,3,s2+(s1-1)*length(sig2s));           
        yrange = log([10^1,10^3.1]);
        ys = linspace(yrange(1),yrange(2),1001);
        loglog(exp(ys),F(ys,p),'--k','LineWidth',1,'DisplayName','A(y)');
        hold on
        set(gca,'FontSize',10);   
        % F is strictly positive therefore we are looking for only positive S1 and S3
        vals = S1(ys,p);       
        loglog(exp(ys(vals>0)),vals(vals>0),'-b','LineWidth',1,'DisplayName','B_1(y)');
        vals = S3(ys,p);       
        loglog(exp(ys(vals>0)),vals(vals>0),'-r','LineWidth',1,'DisplayName','B_3(y)');
        
        % Legend
        if(s1==1 && s2 == length(sig2s))
           l = legend('show');
%            l.Location = 'SouthEast';
           l.Units = 'Inches';
           l.Orientation = 'horizontal';
           l.Position = [0.35 0.06 4.5 0.06];
           l.EdgeColor = [1 1 1];
%            l.Position = [l.Position(1)+l.Position(3)*1.3,l.Position(2),l.Position(3)/1.5,l.Position(4)];
        end
        
        % Asymptotes
        yl = [10^-5, 10^5];
        if(yl(2)<F(p.y2-p.sig2^2,p)), yl(2) = F(p.y2-p.sig2^2,p); end
        if(yl(1)>F(p.y1-p.sig1^2,p)), yl(1) = F(p.y1-p.sig1^2,p); end
        pl =loglog([exp(p.y2-p.sig2^2),exp(p.y2-p.sig2^2)],[S1(p.y2-p.sig2^2-0.01,p),yl(2)],':b','LineWidth',1,'DisplayName','B_1');
        set( get( get( pl, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
        pl=loglog([exp(p.y1-p.sig1^2),exp(p.y1-p.sig1^2)],[yl(1),S1(p.y1-p.sig1^2+0.01,p)],':b','LineWidth',1,'DisplayName','B_1');
        set( get( get( pl, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
        pl=loglog([exp(p.y2),exp(p.y2)],[S3(p.y2-0.01,p),yl(2)],':r','LineWidth',1,'DisplayName','B_3');
        set( get( get( pl, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
        pl=loglog([exp(p.y1),exp(p.y1)],[yl(1),S3(p.y1+0.01,p)],':r','LineWidth',1,'DisplayName','B_3');
        set( get( get( pl, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
       
        % Hartigan's test
        % Draw 10000 "cells" from the distribution and run Hartigan's dip test
%         draw1 = lognrnd(p.y1,p.sig1,5000,1);
%         draw2 = lognrnd(p.y2,p.sig2,5000,1);
%         [dip, pval] = HartigansDipSignifTest([draw1;draw2], 500);
%         
%         yshift = log10(yl(2)/yl(1))/5;
%         t = text(15,10.^ceil(log10(yl(2))-yshift),['$p_U=' num2str(pval,3) '$'],'Interpreter','Latex','FontSize',18,'Color','b');
%         
%         draw1 = randn(5000,1)*p.sig1+p.y1;
%         draw2 = randn(5000,1)*p.sig2+p.y2;
%         [dip, pval] = HartigansDipSignifTest([draw1;draw2], 500);
%         t = text(150,10^floor(log10(yl(1))+yshift),['$p_U=' num2str(pval,3) '$'],'Interpreter','Latex','FontSize',18,'Color','r');
        
%         lnUpwardDivergence = p.y2-p.sig2^2;
%         s1minus = S1(lnUpwardDivergence-0.01,p); % If this is < F(same point) then we must have 3 solutions
%         Fminus = F(lnUpwardDivergence-0.01,p);
%         if(s1minus<Fminus)
% %             disp(exp(lnUpwardDivergence));
%             loglog(exp(lnUpwardDivergence),Fminus,'*b','MarkerSize',22);
%         end
    
%         if(s1==3)
%             xlabel('y');
%         end
        
        xlim(exp(yrange));
%         yl = [F(p.y1,p),3*max([F(p.y2-p.sig2^2,p),S1(p.y2-p.sig2^2,p),S3(p.y2-p.sig2^2,p)])];    
        ylim(yl);
        set(gca,'XTick',[10^2,10^3]);
        set(gca,'YTick',10.^[ceil(log10(yl(1))) floor(log10(yl(2)))] );
        title(['$\sigma_1=' num2str(sig1s(s1)) ', \sigma_2=' num2str(sig2s(s2)) '$'],'Interpreter','Latex');        
%         set(gca,'YTick',yl);
    end
end


figname = 'ExtremaSolutions';
% figname = 'Fig2';
set(gcf,'Name',figname);
set(gcf,'Units','Inches');
set(gcf,'Position',[5 4 5.2 3.5]);
set(gcf,'InvertHardcopy','off');
set(gcf,'Color','w');
set(gcf,'PaperPosition',get(gcf,'Position'));

% Add 'y' arrow
ta = annotation('textarrow');
ta.FontSize = 12;
ta.String='y';
ta.Position = [ 0.7    0.1056    0.08    0.0016]; 
ta = annotation('textarrow');
ta.FontSize = 12;
ta.String='y';
ta.Position = [0.424    0.1056    0.08    0.0016];
ta = annotation('textarrow');
ta.FontSize = 12;
ta.String='y';
ta.Position = [0.154    0.1056    0.08    0.0016];

print(gcf,'-dtiff',['Plots' filesep figname],'-r600');
print(gcf,'-dpng',['Plots' filesep figname],'-r600');

%print(gcf,'-depsc2',['Plots' filesep figname]);


%% Read Robert's quench data and show discrepancy

tab = readtable('../Data/20140920-OT1-dynamics/20140920-OT1-dynamics_Specimen_005_E5_E0.txt','Delimiter','\t');
logbn = [0:0.2:9]';
pos = find(tab.x_PE_Gr_A__PpERK>0);
val = tab.x_PE_Gr_A__PpERK(pos);
PlogI = hist(log(val),logbn)';
PlogI = PlogI / trapz(logbn,PlogI);

figure;
% semilogx(exp(logbn),PlogI,'.-');
plot(exp(logbn),PlogI/max(PlogI),'ok','MarkerFaceColor','red','DisplayName','$P(\log I)$');
set(gca,'XScale','log');
set(gca,'YScale','log');
set(gca,'FontSize',10);
hold on
ft = fit(logbn,PlogI,'gauss2');
p2 = struct;
p2.alpha = ft.a2*ft.c2*sqrt(pi);
p2.y1 = ft.b1;
p2.y2 = ft.b2;
p2.sig1 = ft.c1/sqrt(2);
p2.sig2 = ft.c2/sqrt(2);

LogGaussMix = @(y,p)((1-p.alpha)/(p.sig1*sqrt(2*pi))*exp(-(y-p.y1).^2/(2*p.sig1^2)) + p.alpha/(p.sig2*sqrt(2*pi))*exp(-(y-p.y2).^2/(2*p.sig2^2)));
gaussmixfit = LogGaussMix(logbn,p2);
pl = plot(exp(logbn),gaussmixfit/max(gaussmixfit),'.-r','LineWidth',2);
set( get( get( pl, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );

QI = exp(-logbn).*PlogI;
QI = QI / max(QI(find(logbn>0)));
plot(exp(logbn),QI,'ok','MarkerFaceColor','b','DisplayName','$Q(I)$');
vallinearfit = exp(-logbn).*LogGaussMix(logbn,p2);
norm = max(vallinearfit);
vallinearfit = vallinearfit / norm;
pl = plot(exp(logbn),vallinearfit,'-b','LineWidth',2);
set( get( get( pl, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
ylim([0 1]);
xlim([1 10^4]);
yrange = log([10^0,10^4.5]);
ys = linspace(yrange(1),yrange(2),1001);
% [mn,mind] = min(abs(F(ys,p2)-S1(ys,p2)));
% yval = ys(mind);
% plot(exp(yval),exp(-yval)*LogGaussMix(yval,p2)/norm,'*r');
if(p2.y1<p2.y2)
    plot(exp(p2.y1-p2.sig1^2),1,'*b','MarkerSize',20,'DisplayName','$Q(I_*)$');
else
    plot(exp(p2.y2-p2.sig2^2),1,'*b','MarkerSize',20,'DisplayName','$Q(I_*)$');
end

% Hartigan's dip test
[dip, pval] = HartigansDipSignifTest(log10(val), 500);
t = text(10^(1.1),10^(-2.1),['$p_U=' num2str(pval,3) '$'],'Interpreter',...
    'Latex','FontSize',12,'Color','r');
[dip, pval] = HartigansDipSignifTest(val, 500);
t = text(10^1.1,10^(-1.5),['$p_U=' num2str(pval,3) '$'],'Interpreter',...
    'Latex','FontSize',12,'Color','b');
   


l = legend('show');
l.Interpreter = 'Latex';
l.FontSize = 10;
set(l,'Location','SouthWest');
l.Position = [l.Position(1)+0.1 l.Position(2:4)];
text(10^(2.9),10^(-3.25),'(B)');

figname = 'ExperimentHistograms';
set(gcf,'Name',figname);
set(gcf,'Units','Inches');
set(gcf,'Position',[5 4 1.7 1.7]);
set(gcf,'InvertHardcopy','off');
set(gcf,'Color','w');
set(gcf,'PaperPosition',get(gcf,'Position'));
print(gcf,'-dtiff',['Plots' filesep figname],'-r600');
print(gcf,'-dpng',['Plots' filesep figname],'-r600');
% print(gcf,'-depsc2',['Plots' filesep figname]);
%%
ys = linspace(2,yrange(2),1001);
% ys = [log([19:0.00001:20]), ys];
figure;
loglog(exp(ys),F(ys,p2),'--k','LineWidth',1,'DisplayName','$A$');
set(gca,'FontSize',10);
hold on
loglog(exp(ys),S1(ys,p2),'-b','LineWidth',1,'DisplayName','$B_1$');
loglog(exp(ys),S3(ys,p2),'-r','LineWidth',1,'DisplayName','$B_3$');
ylim([10^-2,10^3]); 
xlim([10,10^4]);
l = legend('show');
l.Interpreter = 'Latex';
l.FontSize = 10;
l.Position = [0.19 0.15 0.4 0.3];
text(10^(3.3),10^(-1.6),'(C)');

figname = 'ExperimentExtremaAnalysis';
set(gcf,'Name',figname);
set(gcf,'Units','Inches');
set(gcf,'Position',[5 4 1.7 1.7]);
set(gcf,'InvertHardcopy','off');
set(gcf,'Color','w');
set(gcf,'PaperPosition',get(gcf,'Position'));
print(gcf,'-dtiff',['Plots' filesep figname],'-r600');
print(gcf,'-dpng',['Plots' filesep figname],'-r600');
% print(gcf,'-depsc2',['Plots' filesep figname],'-r300');


%% Take sig1=0.4, sig2 = 0.8 linear case and exhastively check # events and # bootstraps

% nBootstraps = [100:100:1600];
% nEvents = [2500:2500:2*10^5];
% 
% pVals = zeros(length(nBootstraps),length(nEvents));
% 
% p.sig1 = 0.4;
% p.sig2 = 0.8;
% 
% for nB=1:length(nBootstraps)
%     multiWaitbar('Bootstraps',(nB-1)/length(nBootstraps));
%     for nE = 1:length(nEvents)
%         draw1 = lognrnd(p.y1,p.sig1,nEvents(nE),1);
%         draw2 = lognrnd(p.y2,p.sig2,nEvents(nE),1);
%         [dip, pval] = HartigansDipSignifTest([draw1;draw2], nBootstraps(nB));        
%         pVals(nB,nE) = pval;        
%     end
% end
% 
% multiWaitbar('Bootstraps','Close');
% save('pVals','pVals','p','nBootstraps','nEvents')
%%
load('pVals');
figure;
imagesc(2*nEvents,nBootstraps,pVals);
set(gca,'YDir','Normal');
set(gca,'FontSize',12);
colorbar();
xlabel('Number of Events from Q(I)','Interpreter','Latex');
ylabel('Number of Bootstraps');
figname = 'CheckHartigan';
% figname = 'Fig3';
set(gcf,'Name',figname);
set(gcf,'Units','Inches');
set(gcf,'Position',[5 4 5.2 2.5]);
set(gcf,'InvertHardcopy','off');
set(gcf,'Color','w');
set(gcf,'PaperPosition',get(gcf,'Position'));
print(gcf,'-dtiff',['Plots' filesep figname],'-r600');
print(gcf,'-dpng',['Plots' filesep figname],'-r600');

% print(gcf,'-depsc2',['Plots' filesep figname]);

%% ppERK data from Gregoire's 2008 Science paper

ReadERKData_Amir
