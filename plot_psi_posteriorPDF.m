%%
% Read and plot statistics from the ensemble of models (posterior PDF)
%  (Parametric Slip Inversion, PSI)
%
% Hallo, M., Gallovic, F. (2020): Bayesian self-adapting fault slip inversion 
%  with Green's functions uncertainty and application on the 2016 Mw7.1 Kumamoto
%  earthquake, Journal of Geophysical Research: Solid Earth, 125, e2019JB018703.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Code author: Miroslav Hallo
% Charles University in Prague, Faculty of Mathematics and Physics
% E-mail: hallo@karel.troja.mff.cuni.cz
% Revision 1/2018: The first version of the function
% Revision 2019: Improved graphical representation
%
% Copyright (C) 2018  Miroslav Hallo
%
% This program is published under the GNU General Public License (GNU GPL).
%
% This program is free software: you can modify it and/or redistribute it
% or any derivative version under the terms of the GNU General Public
% License as published by the Free Software Foundation, either version 3
% of the License, or (at your option) any later version.
%
% This code is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY. We would like to kindly ask you to acknowledge the authors
% and don't remove their names from the code.
%
% You should have received copy of the GNU General Public License along
% with this program. If not, see <http://www.gnu.org/licenses/>.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;
clear all;
addpath([pwd,'/lib'])
set(0,'defaulttextinterpreter','tex');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INPUT

% Prefix of the inversion result to plot
fileprefix = 'inv/pop00';

% Prefix for final (saved) figures
prefix = 'PSI_res01';

% Order of fault segments in final plot ( 1 = left->right, -1 = right->left )
order = -1;

% Time diference of two slip-rate time slices (as set by afterprocess.dat) [s]
dsnap = 2.0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% -------------------------------------------------------------------
% Read headers and joint constants
fid = fopen([fileprefix,'_headers.txt'],'r');
tline = fgets(fid);
tline = fgets(fid);
NSeg = str2num(tline);
tline = fgets(fid);
nBins = str2num(tline);
tline = fgets(fid);
sqUvec = str2num(tline);
tline = fgets(fid);
mcount = str2num(tline);
tline = fgets(fid);
corr_thr = str2num(tline);
for kk=1:NSeg
    tline = fgets(fid);
    tline = fgets(fid);
    NL(kk) = str2num(tline);
    tline = fgets(fid);
    NW(kk) = str2num(tline);
    for i=1:13
        tline = fgets(fid);
        allBins(1:nBins,i,kk) = str2num(tline);
    end
    tline = fgets(fid);
    NLc(kk) = str2num(tline);
    tline = fgets(fid);
    NWc(kk) = str2num(tline);
    tline = fgets(fid);
    cellLm(1:NLc(kk),kk) = str2num(tline);
    tline = fgets(fid);
    cellWm(1:NWc(kk),kk) = str2num(tline);
    tline = fgets(fid);
    NL1(kk) = str2num(tline);
    tline = fgets(fid);
    NW1(kk) = str2num(tline);
    tline = fgets(fid);
    posL(1:NL1(kk),kk) = str2num(tline);
    tline = fgets(fid);
    posW(1:NW1(kk),kk) = str2num(tline);
    tline = fgets(fid);
    NPL(kk) = str2num(tline);
    tline = fgets(fid);
    NPW(kk) = str2num(tline);
    tline = fgets(fid);
    cpL(1:NPL(kk),kk) = str2num(tline);
    tline = fgets(fid);
    cpW(1:NPW(kk),kk) = str2num(tline);
end
fclose(fid);
posL = posL/1000;
posW = posW/1000;
cpL = cpL/1000;
cpW = cpW/1000;
cellLm = cellLm/1000;
cellWm = cellWm/1000;

for kk=1:NSeg
    posLcent(1:NL(kk),kk) = posL(1:NL(kk),kk) + (posL(2,kk)-posL(1,kk))/2;
    posWcent(1:NW(kk),kk) = posW(1:NW(kk),kk) + (posW(2,kk)-posW(1,kk))/2;
end

for kk=1:NSeg
    NLint(kk) = NL(kk)*20;
    NWint(kk) = NW(kk)*20;
    cellLint(1:NLint(kk),kk) = posLcent(1,kk) + [0:NLint(kk)-1] * (posLcent(NL(kk),kk)-posLcent(1,kk))/(NLint(kk));
    cellWint(1:NWint(kk),kk) = posWcent(1,kk) + [0:NWint(kk)-1] * (posWcent(NW(kk),kk)-posWcent(1,kk))/(NWint(kk));
    centLint(1:NLint(kk),kk) = cellLint(1:NLint(kk),kk) + (cellLint(2,kk)-cellLint(1,kk))/2;
    centWint(1:NWint(kk),kk) = cellWint(1:NWint(kk),kk) + (cellWint(2,kk)-cellWint(1,kk))/2;
end

display(['Number of models: ',num2str(mcount)])

% Read pop joint
fid = fopen([fileprefix,'_joint.txt'],'r');
for kk=1:NSeg
    for i=1:9
        tline = fgets(fid);
        pop_joint(i,1:nBins,kk) = str2num(tline);
    end
    tline = fgets(fid);
    tline = fgets(fid);
end
fclose(fid);

% Read hypocentre count
fid = fopen([fileprefix,'_hypo2D.txt'],'r');
for kk=1:NSeg
    for i=1:NWc(kk)
        tline = fgets(fid);
        cou_hypo(i,1:NLc(kk),kk) = str2num(tline);
    end
    tline = fgets(fid);
    tline = fgets(fid);
end
fclose(fid);

% Read Neo-Nspline
fid = fopen([fileprefix,'_Nsplines.txt'],'r');
tline = fgets(fid);
tline = fgets(fid);
maxNsplin = str2num(tline);
tline = fgets(fid);
mcount2 = str2num(tline);
for kk=1:NSeg
    tline = fgets(fid);
    tline = fgets(fid);
    neoBins(1:maxNsplin,kk) = str2num(tline);
    tline = fgets(fid);
    pop_neo(1:maxNsplin,kk) = str2num(tline);
end
fclose(fid);


%% -------------------------------------------------------------------
% Plot joint 1
fi = figure('units','normalized','color',[1 1 1],'OuterPosition',[0.2,0.1,0.6,0.8]);
set(fi, 'renderer', 'Painters');
set(fi, 'PaperOrientation', 'portrait');
set(fi, 'PaperPositionMode', 'manual');
set(fi, 'PaperUnits', 'centimeters');
Px = 1 + NSeg*7;
Py = 18;
set(fi,'PaperSize', [Px Py]);
set(fi,'PaperPosition', [0 0 Px Py]);

id = 0;
for i=1:3
    if order==1
        k1 = 1; k2 = 1; k3 = NSeg;
    else
        k1 = NSeg; k2 = -1; k3 = 1;
    end
    for kk=k1:k2:k3
        id = id + 1;
        ax1 = subplot(3,NSeg,id);
        
        if i==1 % Rupture initiation point
            set(gca, 'color', [0.9 0.9 0.9])
            colormap(flipud(hot))
            surf(cellLm(1:NLc(kk),kk),cellWm(1:NWc(kk),kk),cou_hypo(1:NWc(kk),1:NLc(kk),kk)./mcount,'EdgeColor','none');
            axis equal;
            xlim([posL(1,kk) posL(NL(kk),kk)])
            ylim([posW(1,kk) posW(NW(kk),kk)])
            box on;
            grid on;
            view(0,90);
            set(gca, 'Layer', 'top');
            if kk==k3
                h = colorbar('location','eastOutside');
                ylabel(h, 'Probability')
            end
            title(['Rupture initiation #',num2str(kk)])
            xlabel('x - along strike (km)')
            ylabel('y - up dip (km)')
        elseif i==2 % Rupture initiation time
            bins = squeeze(allBins(:,5+1,kk));
            if(sum(bins(:)==bins(1)) == length(bins))
                bins = linspace(1,length(bins),length(bins)) + bins(1);
                bins = bins - floor((length(bins)/2)+1);
                bins = bins - (bins(2)-bins(1))/2;
            end
            binsC = bins + (bins(2)-bins(1))/2;
            bar(binsC,squeeze(pop_joint(1,:,kk))./mcount,'k')
            set(gca,'xlim',[min(bins(:)) max(bins(:))])
            box on; set(gca, 'Layer', 'top');
            title(['Rupture initiation #',num2str(kk)])
            xlabel('Delay to origin (s)')
            ylabel('Probability')
        else % N of splines
            bins = squeeze(neoBins(:,kk));
            bar(bins,squeeze(pop_neo(:,kk))./mcount,'FaceColor',[1 0.4 0.6])
            set(gca,'xlim',[0.5 max(bins(:))+0.5])
            box on; set(gca, 'Layer', 'top');
            title(['Spline points #',num2str(kk)])
            ylabel('Probability')
            xlabel('Number N_\Phi')
        end
    end
end

drawnow
saveas(fi,[prefix,'_postPDF_RUPT.pdf'])


%% -------------------------------------------------------------------
% Plot joint 2
fi = figure('units','normalized','color',[1 1 1],'OuterPosition',[0.2,0.1,0.6,0.8]);
set(fi, 'renderer', 'Painters');
set(fi, 'PaperOrientation', 'portrait');
set(fi, 'PaperPositionMode', 'manual');
set(fi, 'PaperUnits', 'centimeters');
Px = 15;
Py = 5;
set(fi,'PaperSize', [Px Py]);
set(fi,'PaperPosition', [0 0 Px Py]);

ax1 = subplot(1,3,1);
bins = squeeze(allBins(:,5+3,1));
binsC = bins + (bins(2)-bins(1))/2;
binsC = (1-(binsC(:)/sqUvec))*100;
bar(binsC,squeeze(pop_joint(3,:,1))./mcount,'k')
set(gca,'xlim',[min(binsC(:)) max(binsC(:))])
box on; set(gca, 'Layer', 'top');
ylabel('Probability')
xlabel('Data VR (%)')

ax1 = subplot(1,3,2);
bins = squeeze(allBins(:,5+4,1));
binsC = bins + (bins(2)-bins(1))/2;
bar(binsC,squeeze(pop_joint(4,:,1))./mcount,'k')
set(gca,'xlim',[min(bins(:)) max(bins(:))])
box on; set(gca, 'Layer', 'top');
ylabel('Probability')
xlabel('M_0 (Nm)')

ax1 = subplot(1,3,3);
bins = squeeze(allBins(:,5+5,1));
binsC = bins + (bins(2)-bins(1))/2;
bar(binsC,squeeze(pop_joint(5,:,1))./mcount,'k')
set(gca,'xlim',[min(bins(:)) max(bins(:))])
box on; set(gca, 'Layer', 'top');
ylabel('Probability')
xlabel('M_w')

drawnow
saveas(fi,[prefix,'_postPDF_M0.pdf'])


%% -------------------------------------------------------------------
% Plot joint (weighted arithmetic mean)
fi = figure('units','normalized','color',[1 1 1],'OuterPosition',[0.2,0.1,0.6,0.8]);
set(fi, 'renderer', 'Painters');
set(fi, 'PaperOrientation', 'portrait');
set(fi, 'PaperPositionMode', 'manual');
set(fi, 'PaperUnits', 'centimeters');
Px = 12;
Py = 17;
set(fi,'PaperSize', [Px Py]);
set(fi,'PaperPosition', [0 0 Px Py]);

colUse = [0.8594    0.0781    0.2344];

id = 0;
for i=1:4
    if order==1
        k1 = 1; k2 = 1; k3 = NSeg;
    else
        k1 = NSeg; k2 = -1; k3 = 1;
    end
    for kk=k1:k2:k3
        id = id + 1;
        ax1 = subplot(4,NSeg,id);
        
        if(i==1) % Risetime
            bins = squeeze(allBins(:,2,kk));
            binsC = bins + (bins(2)-bins(1))/2;
            bar(binsC,squeeze(pop_joint(6,:,kk))./mcount,'FaceColor',colUse)
            set(gca,'xlim',[min(bins(:)) max(bins(:))])
            box on; set(gca, 'Layer', 'top');
            title(['Average Rise time #',num2str(kk)])
            ylabel('Probability')
            xlabel('s')
        elseif(i==2) % Peaktime
            bins = squeeze(allBins(:,4,kk));
            binsC = bins + (bins(2)-bins(1))/2;
            bar(binsC,squeeze(pop_joint(7,:,kk))./mcount,'FaceColor',colUse)
            set(gca,'xlim',[min(bins(:)) max(bins(:))])
            box on; set(gca, 'Layer', 'top');
            title(['Average Peak time #',num2str(kk)])
            ylabel('Probability')
            xlabel('s')
        elseif(i==3) % Velocity
            bins = squeeze(allBins(:,13,kk));
            binsC = bins + (bins(2)-bins(1))/2;
            bar(binsC/1000,squeeze(pop_joint(9,:,kk))./mcount,'FaceColor',colUse)
            set(gca,'xlim',[min(bins(:)/1000) max(bins(:)/1000)])
            box on; set(gca, 'Layer', 'top');
            title(['Average Rupture velocity #',num2str(kk)])
            ylabel('Probability')
            xlabel('km/s')
        elseif(i==4) % Rake
            bins = squeeze(allBins(:,5,kk));
            edges = [bins/(180/pi); (bins(end)+bins(end)-bins(end-1))/(180/pi)]';
            tmp = squeeze(pop_joint(8,:,kk));
            tmp = tmp / max(tmp);
            azPol = polarhistogram('BinEdges',edges,'BinCounts',tmp);
            set(azPol, 'EdgeColor', 'none');
            set(azPol, 'FaceColor',colUse);
            rlim([0 max(tmp)])
            ax = gca;
            ax.ThetaTick = 0:30:360;
            ax.ThetaTickLabels = {'0','30','60','90','120','150','180','-150','-120','-90','-60','-30'};
            ax.RTickLabel = [];
            title(['  '])
            text(0.5*pi,-2.0,['Average Rake angle #',num2str(kk)],'HorizontalAlign','center')
        end
    end
end

drawnow
saveas(fi,[prefix,'_postPDF_MEAN.pdf'])


%% -------------------------------------------------------------------
% Read Slip pop and plot legend
fid = fopen([fileprefix,'_slip2D.txt'],'r');
bins = allBins(:,1,1);
pop_2D=0;
for kk=1:NSeg
    for p=1:nBins
        tline = fgets(fid);
        for i=1:NW(kk)
            tline = fgets(fid);
            pop_2D(i,1:NL(kk),kk,p) = str2num(tline);
        end
    end
    tline = fgets(fid);
    tline = fgets(fid);
end
fclose(fid);

% Prepare pop plot
nBins2 = nBins; % number of first bins to plot
mycmap =[
0.8594    0.0781    0.2344
0.8594    0.0781    0.2344
];
thBins = 64; % number of edges
width = 2*max(bins(1:nBins2)); % size of circle
thi = linspace(-pi,pi,thBins);
ri = linspace(0,width/2,nBins2);
[TH,R] = meshgrid(thi,ri);
[X,Y] = pol2cart(TH,R);
Z = repmat(zeros(nBins2,1),1,thBins);

% Plot Slip pop legend
figure('units','normalized','color',[1 1 1],'OuterPosition',[0.2,0.1,0.6,0.8]);

for i = 1:3
    subplot(3,1,i);
    colormap(mycmap)
    
    mu = bins(max([1 floor(nBins2*((4-i)/4))]));
    sigm = (nBins2/16)*(1+i/2)*(bins(2)-bins(1));
    tmp = normpdf(bins(1:nBins2),mu,sigm);
    tmp = tmp / max(tmp);
    
    % verze barvenych obruci
    colo = ones(length(tmp),1,3);
    colo(tmp>0.1,1,1) = 0.8594;
    colo(tmp>0.1,1,2) = 0.0781;
    colo(tmp>0.1,1,3) = 0.2344;
    % Maximum posteriory marginal
    colo(tmp>0.9999,1,1) = 0.0;
    colo(tmp>0.9999,1,2) = 0.0;
    colo(tmp>0.9999,1,3) = 0.0;
    
    C = repmat(colo,1,thBins,1);
    
    azSurf = surf(X,Y,Z,C);
    
    set(azSurf, 'EdgeColor', 'none');
    
    set(gca, 'XLim', [-max(bins(1:nBins2)),max(bins(1:nBins2))]);
    set(gca, 'YLim', [-max(bins(1:nBins2)),max(bins(1:nBins2))]);
    xlabel('Slip (m)')
    ylabel('Slip (m)')
    view(2); % flatten to 2-D
    daspect([1 1 1]); % set aspect ratio
    
    hold on
    binsC = bins(1:nBins2) + (bins(2)-bins(1))/2;
    tmp = tmp * (max(bins(1:nBins2))/2);
    bar(binsC,tmp,'facecolor',[0.5 0.5 0.5]);
    hold off
    grid off
    box on; set(gca, 'Layer', 'top');
end
drawnow


%% -------------------------------------------------------------------
% Plot Slip pop plot
for kk=1:NSeg
    fi = figure('units','normalized','color',[1 1 1],'OuterPosition',[0.1,0.05,0.8,0.9]);
    set(fi, 'renderer', 'Painters');
    set(fi, 'PaperOrientation', 'portrait');
    set(fi, 'PaperPositionMode', 'manual');
    set(fi, 'PaperUnits', 'centimeters');
    Px = 19;
    Py = 10;
    set(fi,'PaperSize', [Px Py]);
    set(fi,'PaperPosition', [0 0 Px Py]);
    
    id=0;
    colormap(mycmap)
    for i=NW(kk):-1:1
        for j=1:NL(kk)
            id=id+1;
            subaxis(NW(kk),NL(kk),id,'Spacing',0.005,'Padding',0.0);
            
            tmp = squeeze(pop_2D(i,j,kk,1:nBins2));
            tmp = tmp / max(tmp);
            
            % verze barvenych obruci
            colo = ones(length(tmp),1,3);
            colo(tmp>0.1,1,1) = 0.8594;
            colo(tmp>0.1,1,2) = 0.0781;
            colo(tmp>0.1,1,3) = 0.2344;
            % Maximum posteriory marginal
            colo(tmp>0.9999,1,1) = 0.0;
            colo(tmp>0.9999,1,2) = 0.0;
            colo(tmp>0.9999,1,3) = 0.0;
            
            C = repmat(colo,1,thBins,1);
    
            azSurf = surf(X,Y,Z,C);
            
            set(azSurf, 'EdgeColor', 'none');
            set(gca, 'XLim', [-max(bins(1:nBins2)),max(bins(1:nBins2))]);
            set(gca, 'YLim', [-max(bins(1:nBins2)),max(bins(1:nBins2))]);
            view(2); % flatten to 2-D
            daspect([1 1 1]); % set aspect ratio
            box on; set(gca, 'Layer', 'top');
            grid off
            
            set(gca,'xtickLabel',{})
            set(gca,'ytickLabel',{})
            if (j==1)
                centerW = posW(i,kk) + (posW(2,kk)-posW(1,kk))/2;
                if (i==NW(kk))
                    ylabel('y_i (km)','Rotation',0,'VerticalAlignment','middle','HorizontalAlignment','right')
                else
                    ylabel([num2str(centerW)],'Rotation',0,'VerticalAlignment','middle','HorizontalAlignment','right')
                end
            end
            if (i==1)
                centerL = posL(j,kk) + (posL(2,kk)-posL(1,kk))/2;
                xlabel(num2str(centerL))
                if (j==NL(kk))
                    xlabel('x_i (km)','Rotation',0,'VerticalAlignment','top','HorizontalAlignment','center')
                else
                    xlabel([num2str(centerL)],'Rotation',0,'VerticalAlignment','top','HorizontalAlignment','center')
                end
            end
        end
    end
    
    drawnow
    saveas(fi,[prefix,'_postPDF_SLIP_SEG',num2str(kk),'.pdf'])
end


%% -------------------------------------------------------------------
% Read Peak Slip-Rate pop and plot legend
fid = fopen([fileprefix,'_sliprate2D.txt'],'r');
bins = allBins(:,11,1);
pop_2D=0;
for kk=1:NSeg
    for p=1:nBins
        tline = fgets(fid);
        for i=1:NW(kk)
            tline = fgets(fid);
            pop_2D(i,1:NL(kk),kk,p) = str2num(tline);
        end
    end
    tline = fgets(fid);
    tline = fgets(fid);
end
fclose(fid);

% Prepare pop plot
nBins2 = nBins; % number of first bins to plot
mycmap =[
0.0 0.5 0.7
0.0 0.5 0.7
];
thBins = 64; % number of edges
width = 2*max(bins(1:nBins2)); % size of circle
thi = linspace(-pi,pi,thBins);
ri = linspace(0,width/2,nBins2);
[TH,R] = meshgrid(thi,ri);
[X,Y] = pol2cart(TH,R);
Z = repmat(zeros(nBins2,1),1,thBins);

% Plot Slip-Rate pop legend
figure('units','normalized','color',[1 1 1],'OuterPosition',[0.2,0.1,0.6,0.8]);

for i = 1:3
    subplot(3,1,i);
    colormap(mycmap)
    
    mu = bins(max([1 floor(nBins2*((4-i)/4))]));
    sigm = (nBins2/16)*(1+i/2)*(bins(2)-bins(1));
    tmp = normpdf(bins(1:nBins2),mu,sigm);
    tmp = tmp / max(tmp);
    
    % verze barvenych obruci
    colo = ones(length(tmp),1,3);
    colo(tmp>0.1,1,1) = 0.0;
    colo(tmp>0.1,1,2) = 0.5;
    colo(tmp>0.1,1,3) = 0.7;
    % Maximum posteriory marginal
    colo(tmp>0.9999,1,1) = 0.0;
    colo(tmp>0.9999,1,2) = 0.0;
    colo(tmp>0.9999,1,3) = 0.0;
    
    C = repmat(colo,1,thBins,1);
    
    azSurf = surf(X,Y,Z,C);
    
    set(azSurf, 'EdgeColor', 'none');
    
    set(gca, 'XLim', [-max(bins(1:nBins2)),max(bins(1:nBins2))]);
    set(gca, 'YLim', [-max(bins(1:nBins2)),max(bins(1:nBins2))]);
    xlabel('Slip-rate (m/s)')
    ylabel('Slip-rate (m/s)')
    view(2); % flatten to 2-D
    daspect([1 1 1]); % set aspect ratio
    
    hold on
    binsC = bins(1:nBins2) + (bins(2)-bins(1))/2;
    tmp = tmp * (max(bins(1:nBins2))/2);
    bar(binsC,tmp,'facecolor',[0.5 0.5 0.5]);
    hold off
    grid off
    box on; set(gca, 'Layer', 'top');
end
drawnow


%% -------------------------------------------------------------------
% Plot Slip-Rate pop plot
for kk=1:NSeg
    fi = figure('units','normalized','color',[1 1 1],'OuterPosition',[0.1,0.05,0.8,0.9]);
    set(fi, 'renderer', 'Painters');
    set(fi, 'PaperOrientation', 'portrait');
    set(fi, 'PaperPositionMode', 'manual');
    set(fi, 'PaperUnits', 'centimeters');
    Px = 19;
    Py = 10;
    set(fi,'PaperSize', [Px Py]);
    set(fi,'PaperPosition', [0 0 Px Py]);
    
    id=0;
    colormap(mycmap)
    for i=NW(kk):-1:1
        for j=1:NL(kk)
            id=id+1;
            subaxis(NW(kk),NL(kk),id,'Spacing',0.005,'Padding',0.0);
            
            tmp = squeeze(pop_2D(i,j,kk,1:nBins2));
            tmp = tmp / max(tmp);
            
            % verze barvenych obruci
            colo = ones(length(tmp),1,3);
            colo(tmp>0.1,1,1) = 0.0;
            colo(tmp>0.1,1,2) = 0.5;
            colo(tmp>0.1,1,3) = 0.7;
            % Maximum posteriory marginal
            colo(tmp>0.9999,1,1) = 0.0;
            colo(tmp>0.9999,1,2) = 0.0;
            colo(tmp>0.9999,1,3) = 0.0;
            
            C = repmat(colo,1,thBins,1);
    
            azSurf = surf(X,Y,Z,C);
            
            set(azSurf, 'EdgeColor', 'none');
            set(gca, 'XLim', [-max(bins(1:nBins2)),max(bins(1:nBins2))]);
            set(gca, 'YLim', [-max(bins(1:nBins2)),max(bins(1:nBins2))]);
            view(2); % flatten to 2-D
            daspect([1 1 1]); % set aspect ratio
            box on; set(gca, 'Layer', 'top');
            grid off
            set(gca,'xtickLabel',{})
            set(gca,'ytickLabel',{})
            if (j==1)
                centerW = posW(i,kk) + (posW(2,kk)-posW(1,kk))/2;
                if (i==NW(kk))
                    ylabel('y_i (km)','Rotation',0,'VerticalAlignment','middle','HorizontalAlignment','right')
                else
                    ylabel([num2str(centerW)],'Rotation',0,'VerticalAlignment','middle','HorizontalAlignment','right')
                end
            end
            if (i==1)
                centerL = posL(j,kk) + (posL(2,kk)-posL(1,kk))/2;
                xlabel(num2str(centerL))
                if (j==NL(kk))
                    xlabel('x_i (km)','Rotation',0,'VerticalAlignment','top','HorizontalAlignment','center')
                else
                    xlabel([num2str(centerL)],'Rotation',0,'VerticalAlignment','top','HorizontalAlignment','center')
                end
            end
        end
    end
    
    drawnow
    saveas(fi,[prefix,'_postPDF_SLIP-RATE_SEG',num2str(kk),'.pdf'])
end


%% -------------------------------------------------------------------
% Slip-rate snapshots

bins = allBins(:,11,1);
NTS = 9;
tthr = [1:NTS]*dsnap;

fi = figure('units','normalized','color',[1 1 1],'OuterPosition',[0.1,0.05,0.8,0.9]);
set(fi, 'renderer', 'Painters');
set(fi, 'PaperOrientation', 'portrait');
set(fi, 'PaperPositionMode', 'manual');
set(fi, 'PaperUnits', 'centimeters');
Px = 9;
Py = 16;
set(fi,'PaperSize', [Px Py]);
set(fi,'PaperPosition', [0 0 Px Py]);

% Prepare pop plot
nBins2 = nBins; % number of first bins to plot
thBins = 64; % number of edges
width = 2*max(bins(1:nBins2)); % size of circle
thi = linspace(-pi,pi,thBins);
ri = linspace(0,width/2,nBins2);
[TH,R] = meshgrid(thi,ri);
[X,Y] = pol2cart(TH,R);
Z = repmat(zeros(nBins2,1),1,thBins);

popHbx = max(bins(1:nBins2));  
cnt = 0;
for ts=1:NTS
    % Read Slip-rate snapshot
    fid = fopen([fileprefix,'_sliprate2D_',num2str(ts,'%02i'),'.txt'],'r');
    pop_2D=0;
    for kk=1:NSeg
        for p=1:nBins
            tline = fgets(fid);
            for i=1:NW(kk)
                tline = fgets(fid);
                pop_2D(i,1:NL(kk),kk,p) = str2num(tline);
            end
        end
        tline = fgets(fid);
        tline = fgets(fid);
    end
    fclose(fid);
    
    if order==1
        k1 = 1; k2 = 1; k3 = NSeg;
    else
        k1 = NSeg; k2 = -1; k3 = 1;
    end
    for kk=k1:k2:k3    
        cnt = cnt+1;
        subaxis(NTS,NSeg,cnt,'SpacingVertical',0.005,'SpacingHorizontal',0.0,'MarginTop',0.03,'MarginBottom',0.05)
        set(gca, 'color', [1 1 1])
        for i=1:NW(kk)
            for j=1:NL(kk)
                
                tmp = squeeze(pop_2D(i,j,kk,1:nBins2));
                tmp = tmp / max(tmp);
                
                % verze barvenych obruci
                colo = ones(length(tmp),1,3);
                colo(tmp>0.1,1,1) = 0.0;
                colo(tmp>0.1,1,2) = 0.5;
                colo(tmp>0.1,1,3) = 0.7;
                % Maximum posteriory marginal
                colo(tmp>0.9999,1,1) = 0.0;
                colo(tmp>0.9999,1,2) = 0.0;
                colo(tmp>0.9999,1,3) = 0.0;
                
                C = repmat(colo,1,thBins,1);
                azSurf = surf(X+posLcent(j,kk)*popHbx,Y+posLcent(i,kk)*popHbx,Z,C);
                set(azSurf, 'EdgeColor', 'none');
                view(2); % flatten to 2-D
                daspect([1 1 1]); % set aspect ratio
                hold on
            end
        end
        hold off
        axis equal;
        grid off
        xlim([posL(1,kk) posL(NL(kk),kk)]*popHbx)
        ylim([posW(1,kk) posW(NW(kk),kk)]*popHbx)
        box on; set(gca, 'Layer', 'top');
        view(0,90);
        
        text(posL(NL(kk),kk)*popHbx,posW(NW(kk),kk)*popHbx,[num2str(tthr(ts),'%4.1f'),'s '],'VerticalAlignment','top','HorizontalAlignment','right','fontSize',4)
        
        if (ts==NTS) && kk==k1
            ylabel('Up dip (km)')
            xlabel('Along strike (km)')
            set(gca,'xticklabels',{})
            set(gca,'yticklabels',{})
        elseif ts==NTS
            xlabel('Along strike (km)')
            set(gca,'xticklabels',{})
            set(gca,'yticklabels',{})
        else
            set(gca,'xticklabels',{})
            set(gca,'yticklabels',{})
        end

    end
    
end

drawnow
saveas(fi,[prefix,'_postPDF_SLIP-RATE_TIME.pdf'])


%% -------------------------------------------------------------------
% Rake pop
perc=0.1;
fid = fopen([fileprefix,'_rakes2D.txt'],'r');
pop_2D=0;
for kk=1:NSeg
    for p=1:nBins
        tline = fgets(fid);
        for i=1:NW(kk)
            tline = fgets(fid);
            pop_2D(i,1:NL(kk),kk,p) = str2num(tline);
        end
    end
    tline = fgets(fid);
    tline = fgets(fid);
end
fclose(fid);

% Prepare most probable and std arrays
for kk=1:NSeg
    bins = allBins(:,5,kk);
    binsC = bins + (bins(2)-bins(1))/2;
    for i=1:NW(kk)
        for j=1:NL(kk)
            maxDev = 0;
            pMin = nBins;
            pMax = 1;
            [mV,mI]=max(pop_2D(i,j,kk,:));
            for p=1:nBins
                if (pop_2D(i,j,kk,p)/mV) >= perc
                    if p<pMin
                        pMin = p;
                    end
                    if p>pMax
                        pMax = p;
                    end
                end
            end
            
            maxDev = (pMax-pMin)*(bins(2)-bins(1))/2;
            
            rak = binsC(mI);
            if(rak<-180)
		      rak=rak+360;
		    elseif(rak>180)
		      rak=rak-360;
            end
            
            m_rakes_2D(i,j,kk) = rak;
            d_rakes_2D(i,j,kk) = maxDev;
        end
        
    end
    d_rakes_2D(:,NL(kk)+1,kk) = d_rakes_2D(:,NL(kk),kk);
    d_rakes_2D(NW(kk)+1,:,kk) = d_rakes_2D(NW(kk),:,kk);
end

% Rakes to arrows
for kk=1:NSeg
    Qsize = 0.6;
    for j=1:NW(kk)
        for i=1:NL(kk)
            centerX = posL(i,kk) + (posL(2,kk)-posL(1,kk))/2;
            centerY = posW(j,kk) + (posW(2,kk)-posW(1,kk))/2;
            rak = m_rakes_2D(j,i,kk);
            qdy = sind(rak)*Qsize;
            qdx = cosd(rak)*Qsize;
            qX(j,i,kk) = centerX - qdx/2;
            qY(j,i,kk) = centerY - qdy/2;
            qU(j,i,kk) = qdx;
            qV(j,i,kk) = qdy;
        end 
    end
end

% Plot
for kk=1:NSeg
    fi = figure('units','normalized','color',[1 1 1],'OuterPosition',[0.2,0.1,0.6,0.8]);
    set(fi, 'renderer', 'Painters');
    set(fi, 'PaperOrientation', 'portrait');
    set(fi, 'PaperPositionMode', 'manual');
    set(fi, 'PaperUnits', 'centimeters');
    Px = 14;
    Py = 10;
    set(fi,'PaperSize', [Px Py]);
    set(fi,'PaperPosition', [0 0 Px Py]);
    
    colormap(gca,pink)
    hold on
    quiver(qX(1:NW(kk),1:NL(kk),kk), qY(1:NW(kk),1:NL(kk),kk), qU(1:NW(kk),1:NL(kk),kk), qV(1:NW(kk),1:NL(kk),kk), 0, 'color', 'r' );
    surf( posL(1:NL1(kk),kk), posW(1:NW1(kk),kk), -ones(NW1(kk),NL1(kk)), d_rakes_2D(1:NW1(kk),1:NL1(kk),kk) );
    plot3(repmat(cpL(1:NPL(kk),kk),NPW(kk),1),reshape(repmat(cpW(1:NPW(kk),kk),1,NPL(kk))',[],1),ones(NPL(kk)*NPW(kk),1),'k.','MarkerSize',5);
    hold off
    axis equal;
    xlim([posL(1,kk) posL(NL(kk),kk)])
    ylim([posW(1,kk) posW(NW(kk),kk)])
    box on; set(gca, 'Layer', 'top');
    view(0,90);
    caxis([0 45]);
    h = colorbar;
    ylabel(h, 'Rake scatter (+/- deg)')
    title(['The most probable rake #',num2str(kk)])
    xlabel('x - along strike (km)')
    ylabel('y - up dip (km)')
    
    drawnow
    saveas(fi,[prefix,'_postPDF_RAKE_SEG',num2str(kk),'.pdf'])
end


%% -------------------------------------------------------------------
% Peak slip-rate time pop

perc=0.1;
fid = fopen([fileprefix,'_ratetime2D.txt'],'r');
pop_2D=0;
for kk=1:NSeg
    for p=1:nBins
        tline = fgets(fid);
        for i=1:NW(kk)
            tline = fgets(fid);
            pop_2D(i,1:NL(kk),kk,p) = str2num(tline);
        end
    end
    tline = fgets(fid);
    tline = fgets(fid);
end
fclose(fid);

% Prepare most probable and std arrays
for kk=1:NSeg
    bins = allBins(:,12,kk);
    binsC = bins + (bins(2)-bins(1))/2;
    for i=1:NW(kk)
        for j=1:NL(kk)
            pMin = nBins;
            pMax = 1;
            [mV,mI]=max(pop_2D(i,j,kk,:));
            for p=1:nBins
                if (pop_2D(i,j,kk,p)/mV) >= perc
                    if p<pMin
                        pMin = p;
                    end
                    if p>pMax
                        pMax = p;
                    end
                end
            end
            m_ruptime_2D(i,j,kk) = binsC(mI);
            dl_ruptime_2D(i,j,kk) = binsC(pMin);
            dh_ruptime_2D(i,j,kk) = binsC(pMax);
        end
    end
    m_ruptime_2D(:,NL(kk)+1,kk) = m_ruptime_2D(:,NL(kk),kk);
    m_ruptime_2D(NW(kk)+1,:,kk) = m_ruptime_2D(NW(kk),:,kk);
end

% Plot
for kk=1:NSeg
    fi = figure('units','normalized','color',[1 1 1],'OuterPosition',[0.2,0.1,0.6,0.8]);
    set(fi, 'renderer', 'Painters');
    set(fi, 'PaperOrientation', 'portrait');
    set(fi, 'PaperPositionMode', 'manual');
    set(fi, 'PaperUnits', 'centimeters');
    Px = 14;
    Py = 10;
    set(fi,'PaperSize', [Px Py]);
    set(fi,'PaperPosition', [0 0 Px Py]);
    
    colormap(gca,hot)
    hold on
    surf( posL(1:NL1(kk),kk), posW(1:NW1(kk),kk), -ones(NW1(kk),NL1(kk)), m_ruptime_2D(1:NW1(kk),1:NL1(kk),kk) );
    plot3(repmat(cpL(1:NPL(kk),kk),NPW(kk),1),reshape(repmat(cpW(1:NPW(kk),kk),1,NPL(kk))',[],1),ones(NPL(kk)*NPW(kk),1),'k.','MarkerSize',5);
    hold off
    axis equal;
    xlim([posL(1,kk) posL(NL(kk),kk)])
    ylim([posW(1,kk) posW(NW(kk),kk)])
    box on; set(gca, 'Layer', 'top');
    view(0,90);
    caxis([0 max(binsC)]);
    h = colorbar;
    ylabel(h, 'Peak slip-rate time (s)')
    title(['The most probable peak slip-rate time #',num2str(kk)])
    xlabel('x - along strike (km)')
    ylabel('y - up dip (km)')
    
    drawnow
    saveas(fi,[prefix,'_postPDF_SR-TIME_SEG',num2str(kk),'.pdf'])
end



