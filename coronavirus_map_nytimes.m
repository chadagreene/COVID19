% This script was created by Chad A. Greene of NASA Jet Propulsion 
% Laboratory, in March 2020, to generate maps of COVID-19 cases 
% in the United States. 
% 
% The script can automatically download the latest data from the NY Times GitHub. 
% Data are shown real pretty here:
% https://www.nytimes.com/interactive/2020/us/coronavirus-us-cases.html 
% or you can get into the nitty gritty on the NYT Github here: 
% https://github.com/nytimes/covid-19-data
%
% Known issues: 
% * The NYT data lump all NYC data into one, unassociated with a fips, so 
% it doesn't get mapped here. Kansas City also gets left out. 
% 
%% Enter preferences:

DownloadLatestData = true; 
makegif = true; % make a gif using gif function from Climate Data Tools
makemp4 = true; % makes an mp4 using videoWriter 

FixNYC = true;  % evenly distributes New York City values among their five counties (scaled by population). If this is false, no values will appear in the map for NYC 

%% Download new data if requested: 

if DownloadLatestData
   websave('us-counties.csv','https://raw.githubusercontent.com/nytimes/covid-19-data/master/us-counties.csv'); 
end

%% Load data:

T = readtable('us-counties.csv'); % from NY Times

% Distribute NYC values evenly among the five counties with fips codes: 
if FixNYC
   ind = find(strcmp(table2cell(T(:,2)),'New York City'));
   
   for k = 1:length(ind)
      Ttmp = repmat(T(ind(k),:),[5 1]); 
      Ttmp(:,2) = {'New York';'Kings';'Queens';'Bronx';'Richmond'}; 
      Ttmp(:,4) = {36061;36047;36081;36005;36085}; % fips code
      Ttmp(:,5) = {table2array(Ttmp(1,5))*0.1939;table2array(Ttmp(2,5))*0.3075;table2array(Ttmp(3,5))*0.2713;table2array(Ttmp(4,5))*0.1705;table2array(Ttmp(4,5))*0.0567}; % scaled by population
      Ttmp(:,6) = {table2array(Ttmp(1,6))*0.1939;table2array(Ttmp(2,6))*0.3075;table2array(Ttmp(3,6))*0.2713;table2array(Ttmp(4,6))*0.1705;table2array(Ttmp(4,6))*0.0567}; % scaled by population
   end
   
   % Tack it on at the bottom of the table: (Don't worry about being out of chronological order, and don't worry about deleting the New York City lines)  
   T = [T;Ttmp]; 
end

% Put the important columns of the table into their own arrays: 
Tt = table2array(T(:,1));     % array of times in datetime format
Tfips = table2array(T(:,4));  % county fips codes
Tcases = table2array(T(:,5)); % confirmed cases 
Tdeath = table2array(T(:,6)); % deaths

% load county outline and population data:
load US_counties.mat 

%% Reformat data: 

% Create a daily array of times: 
t = min(Tt):max(Tt); 

% Prellocate matrices with a row for each county and a column for each day:
C = NaN(length(x),length(t)); % confirmed cases 
D = C;                        % deaths

% Relate the fips values in each dataset: 
[~,Locb] = ismember(Tfips,fips);

% Loop through each day and populate matrices 
for k = 1:length(t) 
   
   % Get indices of all data corresponding to this day: 
   ind = Tt==t(k) & Locb>0;
   
   C(Locb(ind),k) = Tcases(ind); 
   D(Locb(ind),k) = Tdeath(ind); 
end

% Scale by population of each county: 
C_per10k = 10000*C./pop; % number of confirmed cases per 10,000 residents
D_per10k = 10000*D./pop; % number of deaths per 10,000 residents

%% Make maps: 

% The patchsc color mapping is NOT dynamic, meaning the color axes must be set
% before plotting and can't be changed afterward. You might try something like:
% 
%  cax = median(log10(C_per10k(:,end)),'omitnan') + 2.5*std(log10(C_per10k(:,end)),'omitnan')*[-1 1];
% 
% to set it automatically, or you can manually set the axis limits to values
% of your liking. It might take some tinker to find color axis values that 
% show the minor spatial variations without clipping too much on either end. 

caxc = log10([0.1 40]);  % case color axis limits
caxd = log10([0.003 3]); % death color axis limits

% Open a figure: 
figure('color','k','pos',[50 50 1600 495])

% Start with the most recent data: 
k=size(C_per10k,2);

% PLOT THE CONFIRMED CASE DATA: 
subplot(1,2,1) 

% Indices of finite data: 
ind = C_per10k(:,k)>0; 

% The patchsc function is in Climate Data Toolbox for Matlab (Greene et al, 2019):
hc = patchsc(x(ind),y(ind),log10(C_per10k(ind,k)),...
   'edgecolor','none','colormap',cmocean('thermal'),'caxis',caxc);

% Plot state boundaries for context: 
plot(statex,statey,'color',0.65*[1 1 1])

% Set axis limits: (sorry Alaska) 
axis([ -2294168    2415097    2820366    5725879])
daspect([1 1 1]) 
axis off 
set(gca,'position',[0 0 0.5 1])

% Set tick marks for the colorbar. (Here I'm setting a wider range than necessary, which should hopefully cover everything ever needed) 
caxticks = [0.001 0.003 0.01 0.03 0.1 0.3 1 3 10 30 100 300 1000 3000 10000 30000]; 

cbc = colorbar('south'); 
cbc.Position = [0.03 0.05 0.15 0.017];
cbc.Color = 0.65*[1 1 1]; 
set(cbc,'fontsize',14,'xtick',log10(caxticks),'xticklabel',num2str(caxticks'),'xaxislocation','top')
xlabel(cbc,'known COVID-19 cases per 10k residents') 

ax(1) = gca; % This will come in handy later. 

% PLOT THE DEATH DATA: 
subplot(1,2,2) 

ind = D_per10k(:,k)>0; 

% The patchsc function is in Climate Data Toolbox for Matlab (Greene et al, 2019):
hd = patchsc(x(ind),y(ind),log10(D_per10k(ind,k)),...
   'edgecolor','none','colormap',cmocean('thermal'),'caxis',caxd);

% Plot state boundaries for context: 
plot(statex,statey,'color',0.65*[1 1 1])

% Set axis limits: (sorry Alaska and other colonies) 
axis([ -2294168    2415097    2820366    5725879])
daspect([1 1 1]) 
axis off 
set(gca,'position',[0.5 0 0.5 1])

cbd = colorbar('south'); 
cbd.Position = [0.53 0.05 0.15 0.017];
cbd.Color = 0.65*[1 1 1]; 
set(cbd,'fontsize',14,'xtick',log10(caxticks/100),'xticklabel',num2str(caxticks'/100),'xaxislocation','top')
xlabel(cbd,'COVID-19 deaths per 10k residents') 

ax(2) = gca; % This will come in handy later. 

% The sgtitle function was introduced in 2018b. Replace with "text" if you need to: 
txt = sgtitle(datestr(t(k),'yyyy-mmm-dd'),...
   'fontsize',18,'color',0.65*[1 1 1],'fontweight','bold');

ntitle({'animation by Chad A. Greene';'data from The New York Times'},...
   'location','se','fontsize',10,'color',0.65*[1 1 1])

linkaxes(ax,'xy') % lets you zoom in on one map and automatically adjust the extents of the other. 

%% Animate: 

if makemp4
   % Initialize an mp4 video: 
   v = VideoWriter(['coronavirus_NYT_',datestr(t(end),'yyyymmdd'),'.mp4'],'MPEG-4'); 
   v.FrameRate=8; % 8 frames per second 
   v.Quality=100; 
   open(v)
   writeVideo(v,getframe(gcf)); % writes the first frame
end

if makegif
   % Initialize a gif and write the first frame (using Climate Data Tools):
   gif(['coronavirus_NYT_',datestr(t(end),'yyyymmdd'),'.gif'],'frame',gcf,'delaytime',1/8)   
end

% Loop through each day: 
for k = 1:length(t) 
   
   % Make the case axes current: 
   axes(ax(1))
   
   % Update case data: 
   ind = C_per10k(:,k)>0; 
   delete(hc(:))
   hc = patchsc(x(ind),y(ind),log10(C_per10k(ind,k)),...
      'edgecolor','none','colormap',cmocean('thermal'),'caxis',caxc);
   uistack(hc,'bottom') % places the county data under state boundaries 
   
   % Make the case axes current: 
   axes(ax(2))
   
   % Update case data: 
   ind = D_per10k(:,k)>0; 
   delete(hd(:))
   hd = patchsc(x(ind),y(ind),log10(D_per10k(ind,k)),...
      'edgecolor','none','colormap',cmocean('thermal'),'caxis',caxd);
   uistack(hd,'bottom') % places the county data under state boundaries 
   
   % Update date: 
   txt.String = datestr(t(k),'yyyy-mmm-dd'); 
   
   % Write the video frame: 
   if makegif 
      gif
   end
   if makemp4
      writeVideo(v,getframe(gcf));
   end
end

if makegif
   % Repeat the last frame: (this for the gif, to give viewers time to see the final pattern) 
   for k = 1:12
      gif
   end
end
if makemp4
   close(v)
end

%% Analyze local trends 
% 

% Enter preferences for analysis: 
dur = 10;     % length of window of investigation (days). Longer durations give more stable numbers, but provide a less current view of what's happening, and a shorter window also allows more counties to be included.  
thresh = 20; % minimum # of cases on first day of the analysis window. The numbers tend to be less stable when starting with fewer than 10 or 20 cases.   

% Convert to datenum for polynomial fit relative to days: 
td = datenum(t); 

% Get indices of the window length of interest: 
indw = (length(td)-dur+1):length(td); 

% Sum of cases in all counties (these numbers are lower than totals for all of US bc some cases haven't been assigned to counties yet.)  
sumC = 10000*sum(C,'omitnan')/sum(pop,'omitnan'); 

% Fit polynomials to log data: (Rates of doubling are 1/p(1) days)
p_national = polyfit(td(indw),log(sumC(indw)),1)
p_seattle = polyfit(td(indw),log(C_per10k(82,indw)),1)
p_miami = polyfit(td(indw),log(C_per10k(2923,indw)),1)

figure
plot(t,sumC,'k','linewidth',2) 
axis tight
box off
ylabel('known cases of COVID-19 per 10k residents') 
set(gca,'yscale','log') 
hold on
plot(t,C_per10k(82,:),'color',[0.22 0.49 0.72]) % King County, WA (Seattle) 
plot(t,C_per10k(2923,:),'color',[0.89 0.1 0.11]) % Miami-Dade, FL

% Plot the linear fits: 
plot(t(indw),exp(polyval(p_national,td(indw))),'k-')
plot(t(indw),exp(polyval(p_seattle,td(indw))),'-','color',[0.22 0.49 0.72])
plot(t(indw),exp(polyval(p_miami,td(indw))),'-','color',[0.89 0.1 0.11])

legend('United States','King, WA (Seattle)','Miami-Dade FL','location','southeast')
legend boxoff 

%% Map the recent doubling rates for all counties

% Preallocate an array for doubling rates: 
r = NaN(size(pop)); 

% Loop through each county: 
for k = 1:length(x)

   % Solve for rates only if day 1 of the analysis window has at least thresh # of cases.  
   if C(k,indw(1))>=thresh 
      pv = polyfit(td(indw),log(C(k,indw)),1); 
      r(k) = 1/pv(1); 
   end
end

ind = isfinite(r); 

figure('color',0.5*[1 1 1],'pos',[50 50 840 495])
hp = patchsc(x(ind),y(ind),r(ind),...
   'edgecolor','none','colormap',cmocean('-balance'),...
   'caxis',1/p_national(1) + 3*[-1 1]);

plot(statex,statey,'color',0.65*[1 1 1])

axis([-2294168    2415097    2820366    5725879])
daspect([1 1 1]) 

axis off 

set(gca,'position',[0 0 1 1])

cb = colorbar('south'); 
cb.Position = [0.05 0.05 0.3 0.02];
cb.Color = 0.15*[1 1 1]; 
set(cb,'xdir','reverse','xaxislocation','top') 
xlabel(cb,['doubling rate (days) of known cases, ',datestr(t(indw(1)),'ddmmm'),'-',datestr(t(indw(end)),'ddmmm')])

