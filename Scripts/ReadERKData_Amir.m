% Reads ERK data and plots

A=dir;
nFiles = 8;
concentrations=logspace(-8,-8-nFiles+1,nFiles);


filenames = cell(nFiles,1);
Data = cell(nFiles,1);
for ff=1:nFiles
    filenames{ff} = ['../Data/FeinermanScience2008/ERK1[072407]/ERK1(OVA)_Tube_00' num2str(ff) '.txt'];
    try
        tabData = readtable(filenames{ff},'Delimiter','\t','ReadVariableNames',false);
    catch
        error(['Failed: file ' filenames{ff}]);
    end
    
    tabData.Properties.VariableNames{1} = 'FCSA';
    tabData.Properties.VariableNames{2} = 'SSCA';
    tabData.Properties.VariableNames{3} = 'CD8a';
    tabData.Properties.VariableNames{4} = 'ERK1';            
    tabData.Properties.VariableNames{5} = 'DAPI';
    tabData.Properties.VariableNames{6} = 'ppERK';
    Data{ff} = tabData;    
end

%% Plot histogram of ppERK
log10bn = [0:0.1:6];
figure;
hold on
cmap = colormap(parula(nFiles));
for ff=1:nFiles
   hst = hist(log10(Data{ff}.ppERK),log10bn);
%    hst = hist(log10(Data{ff}.ERK1),log10bn);
   hst = hst / trapz(log10bn,hst);
   plot(log10bn,hst,'.-','Color',cmap(ff,:),'DisplayName',num2str(concentrations(ff)));
end

%% Kernel density estimate of joint ppERK/ERK1 distribution

% for ff=1:nFiles

xcrit = 120; % This is where to slice the histogram
chosen_dose = 4;

% for ff=1:length(Data)
for ff=chosen_dose:chosen_dose
    data = log10(Data{ff}{:,{'ppERK','ERK1'}});
    [bandwidth,density,X,Y]=kde2d(data);
    % plot the data and the density estimate
    figure;
    set(gcf,'Units','normalized');
    set(gcf,'Color','w');
    set(gcf,'InvertHardcopy','off');
    set(gcf,'Name',['Concentration 10^{' num2str(log10(concentrations(ff))) '}']);
    set(gca,'Visible','off');
    ax = axes('Position',[0.12 0.15 0.6 0.6]);
    
%     plot joined 
    plt = density;
    plt(density<0) = NaN;
    hold on
    imagesc(X(1,:),Y(:,1),plt);
    set(gca,'YDir','Normal');
    xcoord = X(1,find(X(1,:)>log10(xcrit),1)-1);
    plot([xcoord,xcoord],[Y(1,1) Y(end,1)],'--r','LineWidth',2);
    set(gca,'FontSize',12);
%     colormap;
    
    set(gca,'XTick',[1 2 3]);
    set(gca,'XTickLabel',{'10^1','10^2','10^3'});
    xlabel('ppERK','Interpreter','Latex');
    l=ylabel('ERK1','Interpreter','Latex');
    set(gca,'YTick',[2 3]);
    set(gca,'YTickLabel',{'10^2','10^3'});
     xlim([0 3.5]);
     ylim([1.5,3]);
     
     %     plot top marginal
    axtop =  axes('Position', [ax.Position(1) ax.Position(2)+ax.Position(4) ax.Position(3) 0.3*ax.Position(4)]);    
    rbins = [0:0.1:6];
    plt = hist(data(:,1),rbins);
    plt = plt / trapz(rbins,plt);    
    semilogx(10.^(rbins),plt,'LineWidth',1,'Color','r'); 
    hold on
    h = area(10.^(rbins),plt,'LineStyle','None');
    h(1).FaceColor = 'r';
    h(1).FaceAlpha = 0.2;
    h(1).Annotation.LegendInformation.IconDisplayStyle = 'off';

    semilogx([xcrit xcrit],[0 plt(find(rbins>log10(xcrit),1)-1)],'--r','LineWidth',1);
    xlim(10.^[0 3.5]);
    set(gca,'XTick',10.^[2 3]);
    set(gca,'Box','off');
    set(gca,'Visible','off');       
    
    %     plot right marginals
    axright =  axes('Position', [ax.Position(1)+ax.Position(3) ax.Position(2) 0.5*ax.Position(3) ax.Position(4)]);    
    rbins = [0:0.1:6];
    plt = hist(data(:,2),rbins);
    plt = plt / trapz(rbins,plt);
    semilogx(10.^(rbins),plt,'LineWidth',1,'Color','k');
    hold on
    h = area(10.^(rbins),plt,'LineStyle','None');
    h(1).FaceColor = 'k';
    h(1).FaceAlpha = 0.2;
    h(1).Annotation.LegendInformation.IconDisplayStyle = 'off';

   
    
    xlim(10.^[1.5 3]);
%     set(gca,'XTick',10.^[2 3]);
    set(gca,'Box','off');
    set(gca,'Visible','off');
    set(gca,'view',[90 90]);
    set(gca,'XDir','Reverse')

    axes(ax);
    c=colorbar;
    c.Location = 'east';
    c.Ticks = [];

    annotation(gcf,'textbox',...
    [0.12 0.88 0.2 0.03],...
    'String','$P(\log I_{ppERK})$',...
    'Interpreter','Latex',...
    'LineStyle','none',...
    'FontSize',12,...
    'FitBoxToText','off');
    
     annotation(gcf,'textbox',...
    [0.72 0.7 0.2 0.03],...
    'String','$P(\log I_{ERK1})$',...
    'Interpreter','Latex',...
    'LineStyle','none',...
    'FontSize',12,...
    'FitBoxToText','off');

    if(ff==chosen_dose)
%        c = colorbar;
       title('');
       grey = [0.6 0.6 0.6];
       % Slanted gate y = mx + c (in log space)
%        m = 1; % Set slope = 1 for the two to be proportional to each other
       AngularGate = @(x) x + 0.25;
       xgate = [1.3 2.7];
%        ppERK = data(:,1);
%        ERK1 = data(:,2);
%        fnActive = find(log10(Data{ff}.ERK1)<AngularGate(log10(Data{ff}.ppERK)));
%        loglog(log10(Data{ff}.ppERK(fnActive)),log10(Data{ff}.ERK1(fnActive)),'.r');
  
       
       plot(xgate,AngularGate(xgate),'--','LineWidth',2,'Color',grey);
       
%        text(1,1,'(A)','FontSize',32,'Color','r');
       annotation(gcf,'textbox',...
    [0.16 0.16 0.07 0.08],...
    'Color',[1 0.1 0.1],...
    'String','(A)',...
    'LineStyle','none',...
    'FontSize',12,...
    'FontName','Arial',...
    'FitBoxToText','off');


       figname = 'FeinermanHeatmap';
%        figname = 'Fig5a';
       set(gcf,'Name',figname);
       set(gcf,'Units','Inches');
       set(gcf,'Position',[5 4 5.2 3.5]);
       set(gcf,'InvertHardcopy','off');
       set(gcf,'Color','w');
       set(gcf,'PaperPosition',get(gcf,'Position'));
       print(gcf,'-dtiff',['Plots' filesep figname],'-r600');
       print(gcf,'-dpng',['Plots' filesep figname],'-r600');              
%        print(gcf,'-depsc2',['Plots' filesep figname]);
    else 
%         close(gcf);
    end
    
%     contour3(X,Y,density,100), hold on
%     plot(data(:,1),data(:,2),'r.','MarkerSize',5)
end

%% Plot P(ERK1), P(ERK1 | ppERK > xcrit), P(ERK1|Activation diagonal gate)

% Calculate p(Activation) in the process and calculate the mutual information

MutualInformation = zeros(length(Data),2); 
% (1) vs  ppERK <> xcrit
% (2) vs  right/left of diagonal gate

alphas = zeros(length(Data),2);
% percent active
% (1) vs  ppERK <> xcrit
% (2) vs  right/left of diagonal gate

DKLs = zeros(length(Data),2,2);
% (1) vs  ppERK xcrit (-,+)
% (2) vs  right/left of diagonal gate (Active,Inactive);


for ff=1:length(Data)
    figure;  
    
    rbins = [1.2:0.1:3.2];
    % P(ERK1)
    plt = hist(log10(Data{ff}.ERK1),rbins);
    plt = plt / trapz(rbins,plt);
    semilogx(10.^(rbins),plt,'LineWidth',1,'Color','k','DisplayName','$P(\log ERK1)$');
    hold on
    set(gca,'FontSize',10);
    PERK = plt;
    
    % P(ERK1 | ppERK > xcrit)
    fnLargeppERK = find(log10(Data{ff}.ppERK) > log10(xcrit));
    plt = hist(log10(Data{ff}.ERK1(fnLargeppERK)),rbins);
    plt = plt / trapz(rbins,plt);
    pl = semilogx(10.^(rbins),plt,'LineWidth',1,'Color','r','DisplayName','$P(\log ERK1 | ppERK > ppERK_*)$');
%     set( get( get( pl, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
    % Calculate DKL and mutual information
    PERK_givenActive = plt;
    val = PERK_givenActive./PERK;
    useinds = find((~isnan(val)).*(val>0).*(~isinf(val)));
    DKL(ff,1,1) = trapz(rbins(useinds),PERK_givenActive(useinds).*log2(val(useinds)));
%     DKL(ff,1,1) = trapz(10.^rbins(useinds),10.^rbins(useinds).*PERK_givenActive(useinds).*log2(val(useinds)));
    
    plt = hist(log10(Data{ff}.ERK1(setdiff(1:height(Data{ff}),fnLargeppERK))),rbins);
    plt = plt / trapz(rbins,plt);
    semilogx(10.^(rbins),plt,'--','LineWidth',1,'Color','r','DisplayName','$P(\log ERK1 | ppERK < ppERK_*)$');
    hold on
    % Get P(Active) = alpha
    alphas(ff,1) = length(fnLargeppERK)/length(Data{ff}.ppERK);    
    % Calculate DKL and mutual information
    PERK_givenInactive = plt;
    val = PERK_givenInactive./PERK;
    useinds = find((~isnan(val)).*(val>0).*(~isinf(val)));
    DKL(ff,1,2) = trapz(rbins(useinds),PERK_givenInactive(useinds).*log2(val(useinds))); 
%     DKL(ff,1,2) = trapz(10.^rbins(useinds),10.^rbins(useinds).*PERK_givenInactive(useinds).*log2(val(useinds)));
    
    % P(ERK1 | Activation)
    grey = [0.6 0.6 0.6];
    fnActive = find(log10(Data{ff}.ERK1)<AngularGate(log10(Data{ff}.ppERK)));
    alphas(ff,2) = length(fnActive)/length(Data{ff}.ppERK);    
    plt = hist(log10(Data{ff}.ERK1(fnActive)),rbins);
    plt = plt / trapz(rbins,plt);
    pl = semilogx(10.^(rbins),plt,'LineWidth',1,'Color',grey,'DisplayName','$P(\log ERK1 | Active)$');
%     set( get( get( pl, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
    % Calculate DKL and mutual information
    PERK_givenActive = plt;
    val = PERK_givenActive./PERK;
    useinds = find((~isnan(val)).*(val>0).*(~isinf(val)));
    DKL(ff,2,1) = trapz(rbins(useinds),PERK_givenActive(useinds).*log2(val(useinds)));
%     DKL(ff,2,1) = trapz(10.^rbins(useinds),10.^rbins(useinds).*PERK_givenActive(useinds).*log2(val(useinds))); 

    plt = hist(log10(Data{ff}.ERK1(setdiff(1:length(Data{ff}.ERK1),fnActive))),rbins);
    plt = plt / trapz(rbins,plt);
    semilogx(10.^(rbins),plt,'--','LineWidth',1,'Color',grey,'DisplayName','$P(\log ERK1 | Inactive)$');
    hold on
    % Calculate DKL and mutual information
    PERK_givenInactive = plt;
    val = PERK_givenInactive./PERK;
    useinds = find((~isnan(val)).*(val>0).*(~isinf(val)));
    DKL(ff,2,2) = trapz(rbins(useinds),PERK_givenInactive(useinds).*log2(val(useinds)));
%     DKL(ff,2,2) = trapz(10.^rbins(useinds),10.^rbins(useinds).*PERK_givenInactive(useinds).*log2(val(useinds))); 


    
    xlim([10 1000]);
    l = legend('show');
    l.Interpreter='Latex';
    l.FontSize = 8;
    % l.Location='NorthWest';
    l.Position = [0.076    0.666    0.39    0.27];
    figname = 'PERK1Comparison';
    set(gcf,'Name',figname);
    set(gca,'YTick',[]);
    xlabel('$ERK1$','Interpreter','latex');
    ylabel('Frequency','Interpreter','latex');
 
    MutualInformation(ff,1) = alphas(ff,1)*DKL(ff,1,1)+(1-alphas(ff,1))*DKL(ff,1,2);
    MutualInformation(ff,2) = alphas(ff,2)*DKL(ff,2,1)+(1-alphas(ff,2))*DKL(ff,2,2);
    if(ff==chosen_dose)
%         set(gcf,'Units','Inches');
%         set(gcf,'Position',[22 4.8 9.7 8.2]);     
%         set(gcf,'PaperOrientation','Landscape');
%         set(gcf,'PaperPosition',[22 4.8 9.7 8.2]);     
        savefig = gcf;    
    else
        close(gcf);
    end
end
figure(savefig); 
% set(gcf,'Color','w');
% set(gcf,'InvertHardcopy','off');
ax = axes('Position',[0.22 0.4 0.22 0.22]);
loglog(concentrations,MutualInformation(:,1),'.-r');
hold on
loglog(concentrations(chosen_dose),MutualInformation(chosen_dose,1),'or','MarkerSize',8);
set(gca,'FontSize',8);
loglog(concentrations,MutualInformation(:,2),'.-','Color',grey);
loglog(concentrations(chosen_dose),MutualInformation(chosen_dose,2),'o','Color',grey,'MarkerSize',8);
xlim([concentrations(end), concentrations(1)]);
ylabel('Mutual information','Interpreter','latex');
xlabel('Concentration','Interpreter','Latex');

figname = 'PERK1Comparison';
% figname = 'FigS1';
set(gcf,'Name',figname);
set(gcf,'Units','Inches');
set(gcf,'Position',[5 4 5.2 3.5]);
set(gcf,'InvertHardcopy','off');
set(gcf,'Color','w');
set(gcf,'PaperPosition',get(gcf,'Position'));
print(gcf,'-dtiff',['Plots' filesep figname],'-r600');
print(gcf,'-dpng',['Plots' filesep figname],'-r600');
% print(gcf,'-depsc2',['Plots' filesep figname]);


%% Plot log10(ppERK/ERK1) histograms
log10foldbn = [-1.5:0.05:1.5];
figure;
hold on
cmap = colormap(parula(nFiles));
for ff=1:nFiles
   hst = hist(log10(Data{ff}.ppERK./Data{ff}.ERK1),log10foldbn);
   hst = hst / trapz(log10foldbn,hst);
   plot(log10foldbn,hst,'.-','Color',cmap(ff,:),'DisplayName',num2str(concentrations(ff)));
end
close(gcf)
%% Plot log2(ppERK/ERK1) histograms
bn = [0:0.05:5];
figure;
hold on
cmap = colormap(parula(nFiles));
for ff=1:nFiles
   hst = hist(Data{ff}.ppERK./Data{ff}.ERK1,bn);
   hst = hst / trapz(bn,hst);
   plot(bn,hst,'.-','Color',cmap(ff,:),'DisplayName',num2str(concentrations(ff)));
end
close(gcf)
%% Pick difficult dose and check

Dosenum = chosen_dose;
logbns = cell(2,1);
logbns{1} = log(10.^[0:0.1:4]');
logbns{2} = log(10.^[-2:0.05:1]');

xlims = cell(2,1);
xlims{1} = [10,300];
xlims{2} = [0.1,1.5];
ylims = cell(2,1);
ylims{1} = [0.00001,500];
ylims{2} = [0.00001,500];

vals = zeros(height(Data{Dosenum}),2);
vals(:,1) = log(Data{Dosenum}.ppERK);
vals(:,2) = log(Data{Dosenum}.ppERK./Data{Dosenum}.ERK1);

varstr = cell(2,1);
varstr{1} = 'I_{ppERK}';
% varstr{2} = '(I_{ppERK}/I_{ERK1})';
varstr{2} = '\tilde{I}';
yl = [10^-5, 10^5]; % Used for plotting the : at the divergence


Hartiganxpos(1) = 10^(2.8);
Hartiganxpos(2) = 10^(0.8);


figure;
for vv=1:2
    logbn = logbns{vv};
    PlogI = hist(vals(:,vv),logbn)';
    PlogI = PlogI / trapz(logbn,PlogI);
    
    % Plot log-binned data
    subplot(2,2,1+(vv-1)*2);
    
    plot(exp(logbn),PlogI/max(PlogI),'ok','MarkerFaceColor','red','DisplayName',['$P(\log ' varstr{vv} ')$']);
    xlim(10.^[min(logbn),max(logbn)]);
    set(gca,'XScale','log');
%     set(gca,'YScale','log');
    set(gca,'FontSize',12);
%     title('Log-binned ppERK');
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
    pl = plot(exp(logbn),gaussmixfit/max(gaussmixfit),'.-r','LineWidth',1);
    set( get( get( pl, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
    
    QI = exp(-logbn).*PlogI;
    QI = QI / max(QI(2:end));
    plot(exp(logbn),QI,'ok','MarkerFaceColor','b','DisplayName',['$Q(' varstr{vv} ')$']);
    vallinearfit = exp(-logbn).*LogGaussMix(logbn,p2);
    norm = max(vallinearfit(2:end));
    vallinearfit = vallinearfit / norm;
    pl = plot(exp(logbn),vallinearfit,'-b','LineWidth',1);
    set( get( get( pl, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
    xlim([1 3000]/100^(vv-1));
%     xlim(xlims{vv});
    ylim([0 1]);
    l = legend('show');
    l.Interpreter = 'Latex';
    l.Position = [l.Position(1)+0.09,l.Position(2:end)];
    l.FontSize = 9;
    
    [dip, pval] = HartigansDipSignifTest(vals(:,vv), 500);
    t = text(Hartiganxpos(vv),0.58,['$p_U=' num2str(pval,3) '$'],'Interpreter',...
        'Latex','FontSize',12,'Color','r');
    [dip, pval] = HartigansDipSignifTest(exp(vals(:,vv)), 500);
    t = text(Hartiganxpos(vv),0.43,['$p_U=' num2str(pval,3) '$'],'Interpreter',...
        'Latex','FontSize',12,'Color','b');
   
    if(vv==1)
        text(10^(3),0.15,'(A)','FontSize',12);
        text(3.33,-0.12,'$I_{ppERK}$','Interpreter',...
            'Latex','FontSize',12,'Color','k');        
    else
        text(10^(1),0.15,'(C)','FontSize',12);
        text(0.05,-0.15,'$\tilde{I}$','Interpreter',...
            'Latex','FontSize',12,'Color','k');               
    end
    % Plot maxima analysis
    F = @(y,p)(-p.alpha/(p.alpha-1)*(p.sig1^3)/(p.sig2^3)*exp(((y-p.y1).^2)/(2*p.sig1^2)-((y-p.y2).^2)/(2*p.sig2^2)));
    S1 = @(y,p)(-(p.sig1^2+y-p.y1)./(p.sig2^2+y-p.y2));
    S3 = @(y,p)(-(y-p.y1)./(y-p.y2));
    subplot(2,2,2+(vv-1)*2);
    title('Extrema analysis');
    yrange = ([min(logbn),max(logbn)]);
    ys = linspace(yrange(1),yrange(2),1001);
    loglog(exp(ys),F(ys,p2),'--k','LineWidth',1,'DisplayName','$A$');
    set(gca,'FontSize',12);
    hold on
    loglog(exp(ys),S1(ys,p2),'-b','LineWidth',1,'DisplayName','$B_1$');
    loglog(exp(ys),S3(ys,p2),'-r','LineWidth',1,'DisplayName','$B_3$');
       
    l = legend('show');
    l.Interpreter = 'Latex';
    l.Location = 'NorthWest';
    l.FontSize = 9;
    
    pl=loglog([exp(p2.y2-p2.sig2^2),exp(p2.y2-p2.sig2^2)],[S1(p2.y2-p2.sig2^2-0.01,p2),yl(2)],':b','LineWidth',1,'DisplayName','B_1');
    set( get( get( pl, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
    pl=loglog([exp(p2.y1-p2.sig1^2),exp(p2.y1-p2.sig1^2)],[yl(1),S1(p2.y1-p2.sig1^2+0.01,p2)],':b','LineWidth',1,'DisplayName','B_1');
    set( get( get( pl, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
    pl=loglog([exp(p2.y2),exp(p2.y2)],[S3(p2.y2-0.01,p2),yl(2)],':r','LineWidth',1,'DisplayName','B_3');
    set( get( get( pl, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
    pl=loglog([exp(p2.y1),exp(p2.y1)],[yl(1),S3(p2.y1+0.01,p2)],':r','LineWidth',1,'DisplayName','B_3');
    set( get( get( pl, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
    
  
    if(vv==1)
        text(150,10^(-4.1),'(B)','FontSize',12);
        text(25,10^(-6),['$I_{ppERK}$'],'Interpreter',...
            'Latex','FontSize',12,'Color','k');
      
    else
        text(0.9,10^(-4.25),'(D)','FontSize',12);
        t = text(0.2,10^(-6.4),['$\tilde{I}$'],'Interpreter',...
            'Latex','FontSize',12,'Color','k');
%         t.Position = [0.05    -0.18         0];
    end
    
%     ylim([10^-3,10^3]);

    
    xlim(xlims{vv});
    yl = ylims{vv};
    if(yl(2)<F(p2.y2-p2.sig2^2,p2)), yl(2) = F(p2.y2-p2.sig2^2,p2); end
    if(yl(1)>F(p2.y1-p2.sig1^2,p2)), yl(1) = F(p2.y1-p2.sig1^2,p2); end
    ylim(yl);
    
    set(gca,'YTick',[]);
    % l = legend('show');
%     l.Interpreter = 'Latex';
%     l.FontSize = 18;
end

figname = 'FixMismatch';
% figname = 'Fig6';
set(gcf,'Name',figname);
set(gcf,'Units','Inches');
set(gcf,'Position',[5 4 5.2 4]);
set(gcf,'InvertHardcopy','off');
set(gcf,'Color','w');
set(gcf,'PaperPosition',get(gcf,'Position'));
print(gcf,'-dtiff',['Plots' filesep figname],'-r600');
print(gcf,'-dpng',['Plots' filesep figname],'-r600');

% print(gcf,'-depsc2',['Plots' filesep figname]);

%% Get P(z=x2/x1) = \int dx1, |x1| Q_x2(x1*z) with x1=ERK1, x2 = ppERK

%  data = log10(Data{ff}{:,{'ppERK','ERK1'}});
logbn = logbns{1};
X1 = Data{chosen_dose}{:,'ppERK'};
X2 = Data{chosen_dose}{:,'ERK1'};
Z = X2./X1;

PlogI = hist(log(X1),logbn)';
PlogI = PlogI / trapz(logbn,PlogI);
ft = fit(logbn,PlogI,'gauss2');
p2 = struct;
p2.alpha = ft.a2*ft.c2*sqrt(pi);
p2.y1 = ft.b1;
p2.y2 = ft.b2;
p2.sig1 = ft.c1/sqrt(2);
p2.sig2 = ft.c2/sqrt(2);
QFromGaussMix = @(y,p)(exp(-y).*((1-p.alpha)/(p.sig1*sqrt(2*pi))*exp(-(y-p.y1).^2/(2*p.sig1^2)) + p.alpha/(p.sig2*sqrt(2*pi))*exp(-(y-p.y2).^2/(2*p.sig2^2))));
  
Zs = linspace(1,max(Z)*2,1001)';
QZ = zeros(length(Zs),1);
for zz=1:length(Zs)
    QZ(zz) = sum(X2.*QFromGaussMix(log(X2.*Zs(zz)),p2));
end

PlogZ = QZ./Zs;

figure; 
semilogx(Zs,PlogZ,'.-');
hold on
semilogx(Zs,QZ,'.-');

close(gcf);