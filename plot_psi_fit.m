%%
% Shows final fit of real and synthetic data
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

% Plot raw or standartized data? ( 0 = raw, 1 = standartized )
RAWorSTD = 1;

% Order of waveforms in final plot ( 1 = bottom->top, -1 = top->bottom )
order = 1;

% Directory with the inversion result to plot
resdir = 'inv';

% Prefix for final (saved) figures
prefix = 'PSI_res01';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% -------------------------------------------------------------------

% Read data from ASCII files

% Read target (observed) Uvec data
if RAWorSTD==0
    fid = fopen('var_Unez.tmp','r'); % raw data
else
    fid = fopen('var_Unezstd.tmp','r'); % standartized
end
tline = fgets(fid);
ind = 0;
while ischar(tline)
    ind = ind+1;
    Unez(ind,:) = str2num(tline);
    tline = fgets(fid);
end
fclose(fid);
N = length(Unez(1,:));

% Read ML model Uvec data
if RAWorSTD==0
    fid = fopen([resdir,'/best00_Unez.txt'],'r'); % raw data
else
    fid = fopen([resdir,'/best00_Unezstd.txt'],'r'); % standartized
end
tline = fgets(fid);
ind = 0;
while ischar(tline)
    ind = ind+1;
    Unez2(ind,:) = str2num(tline);
    tline = fgets(fid);
end
fclose(fid);

% Read MAP model Uvec data
if RAWorSTD==0
    fid = fopen([resdir,'/map00_Unez.txt'],'r'); % raw data
else
    fid = fopen([resdir,'/map00_Unezstd.txt'],'r'); % standartized
end
tline = fgets(fid);
ind = 0;
while ischar(tline)
    ind = ind+1;
    Unez3(ind,:) = str2num(tline);
    tline = fgets(fid);
end
fclose(fid);

% Read stations
fid = fopen('input/station_list.dat','r');
tline = fgets(fid);
st = 0;
while ischar(tline)
    st = st+1;
    allCell = textscan(tline,'%f %f %s %f %s');
    statLaLoDis(st,1) = allCell{1,1};
    statLaLoDis(st,2) = allCell{1,2};
    statLaLoDis(st,3) = allCell{1,4};
    statName(st,1) = (allCell{1,3});
    tline = fgets(fid);
end
fclose(fid);

% Read station info
fid = fopen('input/station_info.dat','r');
tline = fgets(fid);
st = 0;
while ischar(tline)
    st = st+1;
    allCell = textscan(tline,'%f %f %f %f %f %f %s');
    statUsed(st,1) = allCell{1,1};
    statUsed(st,2) = allCell{1,2};
    statUsed(st,3) = allCell{1,3};
    statUsed(st,4) = allCell{1,4};
    statUsed(st,5) = allCell{1,5};
    statUsed(st,6) = allCell{1,6};
    tline = fgets(fid);
end
fclose(fid);

staMax = length(statLaLoDis(:,1));
t = Unez(:,1)-Unez(1,1);


%% -------------------------------------------------------------------
% Plot waveforms

% figure
fi = figure('units','normalized','color',[1 1 1],'OuterPosition',[0.1,0.1,0.8,0.8]);
set(fi, 'renderer', 'Painters');
set(fi, 'PaperOrientation', 'portrait');
set(fi, 'PaperPositionMode', 'manual');
set(fi, 'PaperUnits', 'centimeters');
Px = 14;
Py = 17;
set(fi,'PaperSize', [Px Py]);
set(fi,'PaperPosition', [0 0 Px Py]);

ax(1) = subaxis(1,3,1,'Spacing',0.02,'Padding',0.0); hold on; title('N-S'); 
ax(2) = subaxis(1,3,2,'Spacing',0.02,'Padding',0.0); hold on; title('E-W');
ax(3) = subaxis(1,3,3,'Spacing',0.02,'Padding',0.0); hold on; title('U-D');

maxAll = max(max(abs(Unez(:,2:end))));

% Plot waveforms
ypos = order;
id = 1;
statNameUsed = {};
for i = 1:staMax
    maxSt = 0;
    cmc = 0;
    for cmp=1:3
        if(statUsed(i,cmp)~=0)
            cmc = cmc + 1;
            maxSt = max([maxSt,max(abs(Unez(:,id+cmc)))]);
            %maxAll = maxSt % Normalize by each station separately
        end
    end
    for cmp=1:3
        if(statUsed(i,cmp)~=0)
            if cmp==1 || (cmp==2 && statUsed(i,1)==0) || (statUsed(i,1)==0 && statUsed(i,2)==0)
                ypos = ypos+order;
                statNameUsed(abs(ypos),1) = statName(i,1);
            end
            axes(ax(cmp));
            id = id + 1;
            % Plot real data
            wave = Unez(:,id);
            wave = wave / (1.0*maxAll);
            plot(t,wave + ypos,'k')
            
            % Plot text
            %text(t(end),ypos,num2str(max(abs(Unez(:,id))),'%7.1e'),'VerticalAlignment','bottom','HorizontalAlignment','right','fontSize',4)
            
            % Plot ML model data
            wave2 = Unez2(:,id);
            wave2 = wave2 / (1.0*maxAll);
            plot(t,wave2 + ypos,'c')
            
            % Plot MAP model data
            wave3 = Unez3(:,id);
            wave3 = wave3 / (1.0*maxAll);
            plot(t,wave3 + ypos,'b')
            
        end
    end
end

% Plot axes
for cmp=1:3
    axes(ax(cmp));
    hold off;
    grid off;
    box on;
    set(gca, 'Layer', 'top');
    ylim([min(order,ypos+order),max(order,ypos+order)])
    if order == 1
        yticks(order+order:1:ypos)
    else
        yticks(ypos:1:order+order)
    end
    
    if cmp==1
        if order == 1
            yticklabels(statNameUsed(2:abs(ypos),1))
        else
            yticklabels(statNameUsed(abs(ypos):-1:2,1))
        end
    else
        yticklabels({})
    end
    
    if cmp==3
        if order == 1
            legend('Data','ML model','MAP model','Location','northeast')
        else
            legend('Data','ML model','MAP model','Location','southeast')
        end
    end
    
    xlabel('Time (s)','FontSize',8)
    set(gca,'FontSize',8)
end
linkaxes(ax)
drawnow;

% Save
if RAWorSTD==0
    stext = 'RAW';
else
    stext = 'STD';
end
saveas(fi,[prefix,'_FIT_',stext,'.pdf'])

