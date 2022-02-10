%% Porkchop Plot from Itterative Lambert Arc Calculation
% C: 13SEP19 LM: 10JAN22
clear; close all; clc; format long g; rpOMLstart();

%% Inputs

% David Eagle's Lambert Code Folder Path:
DELC = '/Users/rohan/dev/mathworks_FileExchange/davideagle_lambert/';       % MACOS   Path 
%DELC = 'C:\Users\rohan\dev\mathworks_FileExchange\davideagle_lambert';       % Windows Path

plt = 1;    % 1=Plot Porkchop

% Bodies
depbdy = '3';                   % NAIF Body ID of departure body (3=Earth)
arrbdy = '4';                   % NAIF Body ID of arrival body (4=Mars)
[ctr_bdy] = mice_bodc2n(0);     % <-- 0=Sun

% Days and Bounds
et1 = cspice_str2et( {'Jun 04, 2005', 'Feb 20, 2007'} );
et2 = et1;
num_of_Pts = 1000;
dvmaxd = 10;
dvmaxa = 15;


dispC3 = true;






%% Lambert Calculation
tic
addpath(genpath(DELC));
pcd = constants();
mu = pcd.Sun.mu;

% cspice_etcal(number) <-- conv to UTC date/time
t1 = (0:num_of_Pts-1) * ( et1(2) - et1(1) )/num_of_Pts + et1(1);
t2 = (0:num_of_Pts-1) * ( et2(2) - et2(1) )/num_of_Pts + et2(1);
pb1 = mice_spkezr(depbdy, t1, 'J2000', 'NONE', ctr_bdy.name );
pb2 = mice_spkezr(arrbdy, t1, 'J2000', 'NONE', ctr_bdy.name );

for i=1:length(t1)
    for j=1:length(t2)
        tof(i,j) = t2(j) - t1(i);
        if tof(i,j) <= 0
            vimag(i,j) = NaN;
            vfmag(i,j) = NaN;
        else
            lambcall = {l0(1,pb1(i).state,pb2(j).state,tof(i,j),mu)};
            vi_mag = lambcall{1,1}(1);
            vf_mag = lambcall{1,1}(2);
            if vi_mag > dvmaxd
                vimag(i,j) = NaN;
            else 
                vimag(i,j) = vi_mag;
            end
            if vf_mag > dvmaxa
                vfmag(i,j) = NaN;
            else 
                vfmag(i,j) = vf_mag;
            end            
        end
    end
end


% TOF Post Processing
for i=1:length(tof)
    for j=1:length(tof)
        if tof(i,j) < 0
            tof(i,j) = NaN;
        else
            tof(i,j) = tof(i,j)/86400;
        end
    end
end

%% Plotting 
if plt==1
    
    t1jd = cspice_et2utc( t1, 'J', 3 ); t1jd = str2num(t1jd(:,4:end));
    t2jd = cspice_et2utc( t2, 'J', 3 ); t2jd = str2num(t2jd(:,4:end));

    for i=1:length(t1jd)
        t1day = datetime(t1jd(i),'ConvertFrom','juliandate');
        t1days(i) = datenum(t1day);
        t2day = datetime(t2jd(i),'ConvertFrom','juliandate');
        t2days(i) = datenum(t2day);
    end

    hold on
    
    
    if dispC3
       vimag = vimag.^2; 
    end
    
    
    [c1,h1] = contour(t1days,t2days,vimag',8,'b');
    clabel(c1,h1,'fontname','courier new','color','b');
    %colorbar
    [C,h] = contour(t1days,t2days,vfmag',8,'k');
    clabel(C,h,'fontname','courier new');
    [CC,hh] = contour(t1days,t2days,tof',10,'r','linewidth',2);
    clabel(CC,hh,'fontname','courier new','color','r');
    datetick('x','dd/mm/yyyy'); datetick('y','dd/mm/yyyy');
    hold off
    ax = gca;
    set(gcf,'Color',[0.95 0.95 0.95]);
    set(gca,'Color',[0.95 0.95 0.95],'fontname','courier new');
    ytickangle(45)
    ax.YAxis.TickLabelFormat = '%.5f';
    ax.XAxis.TickLabelFormat = '%.5f';
    grid on
    legend({'Departure V_\infty (km/s)','Arrival V_\infty (km/s)','Time of Flight (Days)'},'fontsize',14,'location','southeast')
    title(['Single Rev. Type 1 and 2 Transfer Between Bodies: ',depbdy,' and ',arrbdy],'fontsize',16,'fontname','courier new')
    xlabel([depbdy,' Departure Date'],'fontsize',14,'fontname','courier new');
    ylabel([arrbdy,' Arrival Date'],'fontsize',14,'fontname','courier new');
end
toc