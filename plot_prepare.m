%%
% Read and plot files prepared for GF computation
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Code author: Miroslav Hallo
% Charles University in Prague, Faculty of Mathematics and Physics
% E-mail: hallo@karel.troja.mff.cuni.cz
% Revision 1/2018: The first version of the function
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
set(0,'defaulttextinterpreter','tex');

%% -------------------------------------------------------------------

% Read sources.dat
fid = fopen('sources.dat','r');
tline = fgets(fid);
souN = 0;
while ischar(tline)
    souN = souN+1;
    NEZ_source(souN,1:7) = str2num(tline);
    tline = fgets(fid);
end
fclose(fid);

% Read stations.dat
fid = fopen('stations.dat','r');
tline = fgets(fid);
stN = 0;
while ischar(tline)
    stN = stN+1;
    NEZ_station(stN,1:3) = str2num(tline);
    tline = fgets(fid);
end
fclose(fid);

% Read fault.dat
fid = fopen('fault.dat','r');
tline = fgets(fid);
fauN = 0;
fauP = 0;
while ischar(tline)
    fauP = fauP+1;
    if fauP == 1
        fauN = fauN+1;
    end
    NEZ_fault(1:3,fauP,fauN) = str2num(tline);
    if fauP == 5 % Skip line to next fault
    	tline = fgets(fid);
        fauP = 0;
    end
    tline = fgets(fid);
end
fclose(fid);

% Read dirac.dat
fid = fopen('dirac.dat','r');
tline = fgets(fid);
N = 0;
while ischar(tline)
    N = N+1;
    tmp = str2num(tline);
    dirac(N) = complex(tmp(1),tmp(2));
    tline = fgets(fid);
end
fclose(fid);


%% -------------------------------------------------------------------
% Plot map

figure('units','normalized','color',[1 1 1],'OuterPosition',[0.2,0.1,0.6,0.8]);
hold on
plot3(NEZ_source(:,3),NEZ_source(:,2),NEZ_source(:,4),'bo')
plot3(NEZ_station(:,2),NEZ_station(:,1),NEZ_station(:,3),'r^')

%for i = 1 :length(NEZ_source(:,1))
%    text(NEZ_source(i,3),NEZ_source(i,2),NEZ_source(i,4),num2str(i))
%end

for f = 1 :fauN
    plot3(NEZ_fault(2,:,f),NEZ_fault(1,:,f),NEZ_fault(3,:,f),'b-')
end
hold off
set(gca,'ZDir','reverse')
xlabel('Easting')
ylabel('Northing')
zlabel('Depth')
box on
axis equal


%% -------------------------------------------------------------------
% Plot source function and its FFT

figure('units','normalized','color',[1 1 1],'OuterPosition',[0.2,0.1,0.6,0.8]);
subplot(2,1,1)
plot(abs(dirac))
xlabel('Discrete frequencies No.')
title('Amplitude spectrum of S-T function (in-file)')
box on

subplot(2,1,2)
plot(real(ifft(dirac)))
xlabel('Time sample No.')
title('S-T function')
box on


