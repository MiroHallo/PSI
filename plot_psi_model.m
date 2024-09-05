%%
% Read and plot the final slip and slip-rate model
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
fileprefix = 'inv/map00';

% Prefix for final (saved) figures
prefix = 'PSI_res01';

MaxSlip = 6.0; % Maximum slip [m]
MaxSlipRate = 2.5; % Maximum slip-rate [m/s]
MaxTime = 16.0; % Maximum rupture time [s]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Constants
CoMap = [1,1,1;0.941830039024353,0.979084968566895,0.994117617607117;0.883660137653351,0.958169937133789,0.988235294818878;0.825490176677704,0.937254905700684,0.982352972030640;...
    0.767320275306702,0.916339874267578,0.976470589637756;0.709150314331055,0.895424842834473,0.970588207244873;0.650980412960053,0.874509811401367,0.964705884456635;...
    0.592810451984406,0.853594779968262,0.958823561668396;0.534640550613403,0.832679748535156,0.952941179275513;0.476470589637756,0.811764717102051,0.947058796882629;...
    0.418300658464432,0.790849685668945,0.941176474094391;0.360130727291107,0.769934654235840,0.935294151306152;0.301960796117783,0.749019622802734,0.929411768913269;...
    0.324509799480438,0.769607841968536,0.715686261653900;0.347058832645416,0.790196061134338,0.501960813999176;0.369607865810394,0.810784339904785,0.288235306739807;...
    0.392156869173050,0.831372559070587,0.074509806931019;0.419786095619202,0.839037418365479,0.071122996509075;0.447415322065353,0.846702337265015,0.067736186087131;...
    0.475044578313828,0.854367196559906,0.064349375665188;0.502673804759979,0.862032115459442,0.060962568968534;0.530303061008453,0.869696974754334,0.057575758546591;...
    0.557932257652283,0.877361834049225,0.054188951849937;0.585561513900757,0.885026752948761,0.050802141427994;0.613190710544586,0.892691612243652,0.047415331006050;...
    0.640819966793060,0.900356531143189,0.044028520584106;0.668449223041534,0.908021390438080,0.040641713887453;0.696078419685364,0.915686249732971,0.037254903465509;...
    0.723707675933838,0.923351168632507,0.033868093043566;0.751336872577667,0.931016027927399,0.030481284484267;0.778966128826141,0.938680946826935,0.027094475924969;...
    0.806595385074616,0.946345806121826,0.023707665503025;0.834224581718445,0.954010725021362,0.020320856943726;0.861853837966919,0.961675584316254,0.016934046521783;...
    0.889483094215393,0.969340443611145,0.013547237962484;0.917112290859222,0.977005362510681,0.010160428471863;0.944741547107697,0.984670221805573,0.006773618981242;...
    0.972370743751526,0.992335140705109,0.003386809490621;1,1,0;1,0.941176474094391,0;1,0.882352948188782,0;1,0.823529422283173,0;1,0.764705896377564,0;1,0.705882370471954,0;...
    1,0.64705884456634,0;1,0.58823531866073,0;1,0.52941179275512,0;1,0.47058823704719,0;1,0.41176471114158,0;1,0.35294118523598,0;1,0.294117659330368,0;1,0.235294118523598,0;...
    1,0.17647059261799,0;1,0.11764705926180,0;1,0.05882352963090,0;1,0,0;0.9375,0,0;0.875,0,0;0.8125,0,0;0.75,0,0;0.6875,0,0;0.625,0,0;0.5625,0,0;0.5,0,0];


%% -------------------------------------------------------------------
% Read headers first
fid = fopen([fileprefix,'_headers.txt'],'r');
tline = fgets(fid);
tline = fgets(fid);
NSeg = str2num(tline);
for kk=1:NSeg
    tline = fgets(fid);
    tline = fgets(fid);
    NL(kk) = str2num(tline);
    tline = fgets(fid);
    posL(1:NL(kk),kk) = str2num(tline);
    tline = fgets(fid);
    NW(kk) = str2num(tline);
    tline = fgets(fid);
    posW(1:NW(kk),kk) = str2num(tline);
end
fclose(fid);
posL = posL/1000;
posW = posW/1000;


%% -------------------------------------------------------------------
% Read ALL

% Read slip
fid = fopen([fileprefix,'_slip2D.txt'],'r');
for kk=1:NSeg
    for i=1:NW(kk)
        tline = fgets(fid);
        slip_2D(i,1:NL(kk),kk) = str2num(tline);
    end
    tline = fgets(fid);
    tline = fgets(fid);
end
fclose(fid);

% Read rakes
fid = fopen([fileprefix,'_rakes2D.txt'],'r');
for kk=1:NSeg
    for i=1:NW(kk)
        tline = fgets(fid);
        rakes_2D(i,1:NL(kk),kk) = str2num(tline);
    end
    tline = fgets(fid);
    tline = fgets(fid);
end
fclose(fid);

% Rakes to arrows
for kk=1:NSeg
    if NSeg==1
        Qsize = 0.7;
    else
        Qsize = 1.2;
    end
    for j=1:NW(kk)-1
        for i=1:NL(kk)-1
            centerX = posL(i,kk) + (posL(2,kk)-posL(1,kk))/2;
            centerY = posW(j,kk) + (posW(2,kk)-posW(1,kk))/2;
            rak = rakes_2D(j,i,kk);
            qdy = sind(rak)*Qsize;
            qdx = cosd(rak)*Qsize;
            qX(j,i,kk) = centerX - qdx/2;
            qY(j,i,kk) = centerY - qdy/2;
            qU(j,i,kk) = qdx;
            qV(j,i,kk) = qdy;
        end 
    end
end

% Read ruptime
fid = fopen([fileprefix,'_ruptime2D.txt'],'r');
for kk=1:NSeg
    for i=1:NW(kk)
        tline = fgets(fid);
        ruptime_2D(i,1:NL(kk),kk) = str2num(tline);
    end
    tline = fgets(fid);
    tline = fgets(fid);
end
fclose(fid);

% Read risetime
fid = fopen([fileprefix,'_risetime2D.txt'],'r');   
for kk=1:NSeg
    for i=1:NW(kk)
        tline = fgets(fid);
        risetime_2D(i,1:NL(kk),kk) = str2num(tline);
    end
    tline = fgets(fid);
    tline = fgets(fid);
end
fclose(fid);

% Read peaktime
fid = fopen([fileprefix,'_peaktime2D.txt'],'r');
for kk=1:NSeg
    for i=1:NW(kk)
        tline = fgets(fid);
        peaktime_2D(i,1:NL(kk),kk) = str2num(tline);
    end
    tline = fgets(fid);
    tline = fgets(fid);
end
fclose(fid);

% Read Peak Slip-rate
fid = fopen([fileprefix,'_sliprate2D.txt'],'r');
for kk=1:NSeg
    for i=1:NW(kk)
        tline = fgets(fid);
        sliprate_2D(i,1:NL(kk),kk) = str2num(tline);
    end
    tline = fgets(fid);
    tline = fgets(fid);
end
fclose(fid);

% Read Peak Slip-rate time
fid = fopen([fileprefix,'_ratetime2D.txt'],'r');
for kk=1:NSeg
    for i=1:NW(kk)
        tline = fgets(fid);
        ratetime_2D(i,1:NL(kk),kk) = str2num(tline);
    end
    tline = fgets(fid);
    tline = fgets(fid);
end
fclose(fid);

% Read joint
fid = fopen([fileprefix,'_joint.txt'],'r');
tline = fgets(fid);
tline = fgets(fid);
hypoXY(1,1:NSeg) = str2num(tline);
tline = fgets(fid);
hypoXY(2,1:NSeg) = str2num(tline);
fclose(fid);
hypoXY = hypoXY/1000;

% Display out
disp('Maximim and average slip [m]:')
display([ num2str(max(slip_2D(:))),' ', num2str(mean(slip_2D(:))) ])

disp('Maximim and average slip-rate [m/s]:')
display([ num2str(max(sliprate_2D(:))),' ', num2str(mean(sliprate_2D(:))) ])

disp('Average (weighted) rake angle [deg]:')
for kk=1:NSeg
    avrVal = 0;
    avrNorm = 0;
    for j=1:NW(kk)
        for i=1:NL(kk)
            rak = rakes_2D(j,i,kk);
            amp = slip_2D(j,i,kk);
            if rak<0
                rak=180+(180+rak);
            end
            avrVal = avrVal + (rak*amp);
            avrNorm = avrNorm + amp;
        end
    end
    avrRes = avrVal / avrNorm;
    if avrRes>180
        avrRes = avrRes-360;
    end
    display(['Rake #',num2str(kk),': ',num2str(avrRes)])
end


%% ------------------------------------------------------------------------
% Plot SLIP and SLIP-RATE figure

for kk=1:NSeg
    fi = figure('units','normalized','color',[1 1 1],'OuterPosition',[0.1,0.1,0.8,0.8]);
    set(fi, 'renderer', 'Painters');
    set(fi, 'PaperOrientation', 'portrait');
    set(fi, 'PaperPositionMode', 'manual');
    set(fi, 'PaperUnits', 'centimeters');
    Px = 25;
    Py = 17;
    set(fi,'PaperSize', [Px Py]);
    set(fi,'PaperPosition', [0 0 Px Py]);
    
    subplot(2,1,1)
    set(gca, 'color', [0.9 0.9 0.9])
    colormap(gca,CoMap)
    Z = zeros(NW(kk),NL(kk));
    hold on
    surf( posL(1:NL(kk),kk), posW(1:NW(kk),kk), Z, slip_2D(1:NW(kk),1:NL(kk),kk) );
    quiver(qX(1:NW(kk)-1,1:NL(kk)-1,kk), qY(1:NW(kk)-1,1:NL(kk)-1,kk), qU(1:NW(kk)-1,1:NL(kk)-1,kk), qV(1:NW(kk)-1,1:NL(kk)-1,kk), 0, 'color', 'b' );
    plot(hypoXY(1,kk),hypoXY(2,kk),'*','MarkerFaceColor','b','MarkerEdgeColor','b','MarkerSize',5)
    hold off
    axis equal;
    xlim([posL(1,kk) posL(NL(kk),kk)])
    ylim([posW(1,kk) posW(NW(kk),kk)])
    box on;
    grid off;
    view(0,90);
    set(gca, 'Layer', 'top');
    h = colorbar('location','eastoutside');
    ylabel(h, 'Slip (m)')
    caxis([0 MaxSlip])
    if NSeg==1
        title('Slip')
    else
        title(['Slip #',num2str(kk)])
    end
    xlabel('x - along strike (km)')
    ylabel('y - up dip (km)')
    
    subplot(2,1,2)
    set(gca, 'color', [0.9 0.9 0.9])
    colormap(gca,CoMap)
    hold on
    surf( posL(1:NL(kk),kk), posW(1:NW(kk),kk), Z, sliprate_2D(1:NW(kk),1:NL(kk),kk));
    plot(hypoXY(1,kk),hypoXY(2,kk),'*','MarkerFaceColor','b','MarkerEdgeColor','b','MarkerSize',5)
    hold off
    axis equal;
    xlim([posL(1,kk) posL(NL(kk),kk)])
    ylim([posW(1,kk) posW(NW(kk),kk)])
    box on;
    grid off;
    view(0,90);
    set(gca, 'Layer', 'top');
    h = colorbar('location','eastoutside');
    ylabel(h, 'Peak slip-rate (m/s)')
    caxis([0 MaxSlipRate])
    if NSeg==1
        title('Slip-rate')
    else
        title(['Slip-rate #',num2str(kk)])
    end
    xlabel('x - along strike (km)')
    ylabel('y - up dip (km)')
    
    drawnow
    saveas(fi,[prefix,'_SLIP_SEG',num2str(kk),'.pdf'])
end


%% -------------------------------------------------------------------
% Plot SLIP-RATE snapshots

% Prepare Slip-Rate functions
dt = 0.2;
N = (MaxTime/dt)+1;
time = [0:N-1]*dt;
SRfunc = zeros(N,NW(kk),NL(kk),kk);
for kk=1:NSeg
    for j=1:NW(kk)
        for i=1:NL(kk)
            SRfunc(1:N,j,i,kk) = RYoffe(N,dt,ruptime_2D(j,i,kk),risetime_2D(j,i,kk),peaktime_2D(j,i,kk),slip_2D(j,i,kk));
        end
    end
end

% Prepare time slices
NTS = 8; % Number of time slices
tthr = [1:NTS]*(MaxTime/NTS);
Nthr = floor([1:NTS]*(N/NTS));
for kk=1:NSeg
    for j=1:NW(kk)
        for i=1:NL(kk)
            for ts=1:NTS
                slip_TS_2D(j,i,kk,ts) = SRfunc(Nthr(ts),j,i,kk);
            end
        end
    end
end

% Plot slip-Rate snapshots
for kk=1:NSeg
    fi = figure('units','normalized','color',[1 1 1],'OuterPosition',[0.1,0.1,0.8,0.8]);
    set(fi, 'renderer', 'Painters');
    set(fi, 'PaperOrientation', 'portrait');
    set(fi, 'PaperPositionMode', 'manual');
    set(fi, 'PaperUnits', 'centimeters');
    Px = 9;
    Py = 16;
    set(fi,'PaperSize', [Px Py]);
    set(fi,'PaperPosition', [0 0 Px Py]);
    
    cnt = 0;
    for ts=1:NTS
        cnt = cnt+1;
        subaxis(NTS,1,cnt,'SpacingVertical',0.005,'SpacingHorizontal',0.0,'MarginTop',0.03,'MarginBottom',0.05)
        
        set(gca, 'color', [0.9 0.9 0.9])
        colormap(gca,CoMap)
        Z = zeros(NW(kk),NL(kk));
        hold on;
        surf( posL(1:NL(kk),kk), posW(1:NW(kk),kk), Z, slip_TS_2D(1:NW(kk),1:NL(kk),kk,ts));
        shading flat;
        plot(hypoXY(1,kk),hypoXY(2,kk),'*','MarkerFaceColor','b','MarkerEdgeColor','b','MarkerSize',3)
        hold off;
        axis equal;
        %xlim([-0.1 34.1])
        %ylim([-0.1 18.1])
        xlim([posL(1,kk) posL(NL(kk),kk)])
        ylim([posW(1,kk) posW(NW(kk),kk)])
        box on;
        view(0,90);
        grid off;
        set(gca, 'Layer', 'top');
        caxis([0 MaxSlipRate])
        
        text(posL(NL(kk),kk),posW(NW(kk),kk),[num2str(tthr(ts),'%4.1f'),'s '],'VerticalAlignment','top','HorizontalAlignment','right','fontSize',4)
        
        if ts==NTS
            xlabel('Along strike (km)')
            ylabel('Up dip (km)')
        else
            set(gca,'xticklabels',{})
            set(gca,'yticklabels',{})
        end
    end
    
    drawnow
    saveas(fi,[prefix,'_SLIP-RATE_SEG',num2str(kk),'.pdf'])
end


%% -------------------------------------------------------------------
% Plot SLIP-RATE additional plots

for kk=1:NSeg
    
    fi = figure('units','normalized','color',[1 1 1],'OuterPosition',[0.1,0.1,0.8,0.8]);
    set(fi, 'renderer', 'Painters');
    set(fi, 'PaperOrientation', 'portrait');
    set(fi, 'PaperPositionMode', 'manual');
    set(fi, 'PaperUnits', 'centimeters');
    Px = 25;
    Py = 17;
    set(fi,'PaperSize', [Px Py]);
    set(fi,'PaperPosition', [0 0 Px Py]);
    
    subplot(2,2,1)
    set(gca, 'color', [0.9 0.9 0.9])
    colormap(gca,CoMap)
    Z = zeros(NW(kk),NL(kk));
    hold on
    surf( posL(1:NL(kk),kk), posW(1:NW(kk),kk), Z, slip_2D(1:NW(kk),1:NL(kk),kk) );
    plot(hypoXY(1,kk),hypoXY(2,kk),'*','MarkerFaceColor','b','MarkerEdgeColor','b','MarkerSize',5)
    for i=NW(kk):-1:1
        for j=1:NL(kk)
            Xtmp = posL(j,kk) + (posL(2,kk)-posL(1,kk)) .* (time-time(1))./(time(end)-time(1));
            Ytmp = posW(i,kk) + (posW(2,kk)-posW(1,kk)) .* (SRfunc(1:N,i,j,kk)./MaxSlipRate);
            plot(Xtmp,Ytmp,'b')
        end
    end
    
    hold off
    axis equal;
    xlim([posL(1,kk) posL(NL(kk),kk)])
    ylim([posW(1,kk) posW(NW(kk),kk)])
    box on;
    grid off;
    view(0,90);
    set(gca, 'Layer', 'top');
    h = colorbar('location','eastoutside');
    ylabel(h, 'Total slip (m)')
    caxis([0 MaxSlip])
    if NSeg==1
        title('Slip')
    else
        title(['Slip #',num2str(kk)])
    end
    xlabel('x - along strike (km)')
    ylabel('y - up dip (km)')
    
    subplot(2,2,2)
    set(gca, 'color', [0.9 0.9 0.9])
    hold on;
    surf( posL(1:NL(kk),kk), posW(1:NW(kk),kk), Z, ruptime_2D(1:NW(kk),1:NL(kk),kk));
    plot(hypoXY(1,kk),hypoXY(2,kk),'*','MarkerFaceColor','b','MarkerEdgeColor','b','MarkerSize',5)
    hold off;
    colormap(gca,hot)
    axis equal;
    xlim([posL(1,kk) posL(NL(kk),kk)])
    ylim([posW(1,kk) posW(NW(kk),kk)])
    box on;
    grid off;
    view(0,90);
    set(gca, 'Layer', 'top');
    h = colorbar('location','eastoutside');
    ylabel(h, 'Rupture time (s)')
    caxis([0 MaxTime])
    title('Rupture time')
    xlabel('x - along strike (km)')
    ylabel('y - up dip (km)')
    
    subplot(2,2,3)
    set(gca, 'color', [0.9 0.9 0.9])
    surf( posL(1:NL(kk),kk), posW(1:NW(kk),kk), risetime_2D(1:NW(kk),1:NL(kk),kk));
    colormap(gca,redblue)
    axis equal;
    xlim([posL(1,kk) posL(NL(kk),kk)])
    ylim([posW(1,kk) posW(NW(kk),kk)])
    box on;
    grid off;
    view(0,90);
    set(gca, 'Layer', 'top');
    h = colorbar('location','eastoutside');
    ylabel(h, 'Rise time (s)')
    title('Rise time')
    xlabel('x - along strike (km)')
    ylabel('y - up dip (km)')
    
    subplot(2,2,4)
    surf( posL(1:NL(kk),kk), posW(1:NW(kk),kk), peaktime_2D(1:NW(kk),1:NL(kk),kk));
    colormap(gca,redblue)
    axis equal;
    xlim([posL(1,kk) posL(NL(kk),kk)])
    ylim([posW(1,kk) posW(NW(kk),kk)])
    box on;
    grid off;
    view(0,90);
    set(gca, 'Layer', 'top');
    h = colorbar('location','eastoutside');
    ylabel(h, 'Peak time (s)')
    title('Peak time')
    xlabel('x - along strike (km)')
    ylabel('y - up dip (km)')
    
    drawnow
    saveas(fi,[prefix,'_SRF_SEG',num2str(kk),'.pdf'])
end



