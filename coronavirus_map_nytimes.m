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
% * As of Jan 2021, this code has gotten a somewhat messy as I've made a bunch of little
% tweaks over the course of the past year. Sorry about that. 
% 
%% Enter preferences:

DownloadLatestData = true; 
makegif = false; % make a gif using gif function from Climate Data Tools
makemp4 = true; % makes an mp4 using videoWriter 
use_export_fig = true; % export_fig makes higher resolution animations than the default getframe, but takes longer to make the animations.  
newcases = true; % Analyze and plot "new cases" meaning the number of new cases each day rather than cumulative sums. 
despike = true; % removes single-day spikes via a  3 day moving median after accounting for the typical weekly cycle. (This effectively throws out all the cases that get reported when a state dumps a giant load in a single day, and that usually happens when a state finds or starts counting a type of case they hadnt previously reported. Anyway, this option discards those major data dumps in favor of creating a smooth time series without distracting flashes of data in the visualization.   
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
if newcases
   isnc = isnan(C); 
   isnd = isnan(D); 
   C(isnd) = 0; 
   D(isnc) = 0; 
   C_per10k = 10000*diff(movmean([zeros(length(x),1) C],[6 0],2),1,2)./pop; 
   D_per10k = 10000*diff(movmean([zeros(length(x),1) D],[6 0],2),1,2)./pop;
   %C_per10k(C_per10k==0) = NaN; 
   %D_per10k(D_per10k==0) = NaN; 

   if despike 
      day = repmat(1:7,1,100); 
      day = day(1:length(t)); 
      Cd = diff([zeros(length(x),1) C],1,2);
      Cdm = movmean(Cd,7,2); 

      % Remove the weekly cycle by assuming that for each county, each day of week has a level of reporting that is a scalar multiple of that week's average. For example, Sundays might typically have half the cases of the overall week average, but Monday might have double: 
      Cd2 = NaN(size(C));
      for cnty = 1:size(C,1)
         for k = 1:7
            ind = day==k; 
            tmp = median(Cd(cnty,ind)./Cdm(cnty,ind),'omitnan'); 

            Cd2(cnty,ind) = Cd(cnty,ind) - Cd(cnty,ind).*tmp; 
         end
      end
      
      % Despike: 
      Cd2 = movmedian(Cd,3,2); 
      
      % Smooth ith a gaussian filter. Previously I was using a 7 day unweighted moving average, but such a boxcar window results in jittery data because the 7 day average is always sensitive to each new datapoint the boxcar encounters.  
      Cd2f = filter(gausswin(15)/sum(gausswin(15)),1,Cd2,[],2); 
      
      C_per10k = 10000*Cd2f./pop; 
      %C_per10k = 10000*movmean(Cd2,[6 0],2)./pop; 

   end
   
else
   % Cumulative: 
   C_per10k = 10000*C./pop; % number of confirmed cases per 10,000 residents
   D_per10k = 10000*D./pop; % number of deaths per 10,000 residents
end

%% Make maps: 

% The patchsc color mapping is NOT dynamic, meaning the color axes must be set
% before plotting and can't be changed afterward. You might try something like:
% 
%  cax = median(log10(C_per10k(:,end)),'omitnan') + 2.5*std(log10(C_per10k(:,end)),'omitnan')*[-1 1];
% 
% to set it automatically, or you can manually set the axis limits to values
% of your liking. It might take some tinker to find color axis values that 
% show the minor spatial variations without clipping too much on either end. 

linear = false; 

if newcases
   if linear 
      caxc = [0 10]; 
   else
      caxc = log10([3 150]/10);  % case color axis limits
      caxd = log10([0.01 10]/10); % death color axis limits
   end
else
   caxc = log10([0.1 100]);  % case color axis limits
   caxd = log10([0.01 10]); % death color axis limits
end

% Open a figure: 
sc = 1; % size scaling factor. bigger number makes a bigger figure
figure('color','k','pos',[50 50 800*sc 495*sc])

% Start with the most recent data: 
k=size(C_per10k,2);

cmap = cmocean('thermal'); 

% % The patchsc function is in Climate Data Toolbox for Matlab (Greene et al, 2019):
% hc = patchsc(x(ind),y(ind),log10(C_per10k(ind,k)),...
%    'edgecolor','none','colormap',cmocean('thermal'),'caxis',caxc);

% Update case data: 
ind = C_per10k(:,k)>0; 

% The patchsc function is in Climate Data Toolbox for Matlab (Greene et al, 2019):
if linear 
   hc = patchsc(x(ind),y(ind),(C_per10k(ind,k)),...
      'edgecolor','none','colormap',cmocean('thermal'),'caxis',caxc);
else 
   hc = patchsc(x(ind),y(ind),log10(C_per10k(ind,k)),...
   'edgecolor','none','colormap',cmocean('thermal'),'caxis',caxc);
end

% Plot state boundaries for context: 
plot(statex,statey,'color',0.65*[1 1 1])

% Set axis limits: (sorry Alaska) 
axis([ -2294168    2415097    2820366    5725879])
daspect([1 1 1]) 
axis off 
set(gca,'position',[0 0 1 1])

% Set tick marks for the colorbar. (Here I'm setting a wider range than necessary, which should hopefully cover everything ever needed) 
caxticks = [0.001 0.003 0.01 0.03 0.1 0.3 1 3 10 30 100 300 1000 3000 10000 30000]; 

cbc = colorbar('south'); 
cbc.Position = [0.05 0.05 0.3 0.017];
cbc.Color = 0.65*[1 1 1]; 
set(cbc,'fontsize',14*sc,'xtick',log10(caxticks),'xticklabel',num2str(caxticks'),'xaxislocation','top')
xlabel(cbc,'daily COVID-19 cases per 10k residents') 

ax(1) = gca; % This will come in handy later. 

% The sgtitle function was introduced in 2018b. Replace with "text" if you need to: 
txt = ntitle(datestr(t(k),'mmm dd'),...
   'fontsize',20*sc,'color',0.65*[1 1 1],'fontweight','bold');

ntitle({'animation by Chad A. Greene';'data from The New York Times'},...
   'location','se','fontsize',10*sc,'color',0.65*[1 1 1])

%% Animate: 

if makemp4
   % Initialize an mp4 video: 
   v = VideoWriter(['coronavirus_NYT_daily_',datestr(t(end),'yyyymmdd'),'despike.mp4'],'MPEG-4'); 

   v.FrameRate=16; %
   v.Quality=100; 
   open(v)
   if use_export_fig
      fr = export_fig('-nocrop','-r500'); % can use getframe(gcf) instead, but export_fig offers more control over resolution, etc.  
   else
      fr = getframe(gcf); 
   end
   writeVideo(v,fr); % writes the first frame
end

if makegif
   % Initialize a gif and write the first frame (using Climate Data Tools):
   gif(['coronavirus_NYT_daily_',datestr(t(end),'yyyymmdd'),'despike2.gif'],'frame',gcf,'delaytime',1/16)   
end

% Loop through each day: 
for k = 1:length(t) 
   
   delete(hc(:))
   
   % Update case data: 
   ind = C_per10k(:,k)>0; 

   % The patchsc function is in Climate Data Toolbox for Matlab (Greene et al, 2019):
   hc = patchsc(x(ind),y(ind),log10(C_per10k(ind,k)),...
      'edgecolor','none','colormap',cmocean('thermal'),'caxis',caxc);

   % Update date: 
   txt.String = datestr(t(k),'mmm dd'); 
   uistack(hc,'bottom')
   
   % Write the video frame: 
   if makegif 
      gif
   end
   if makemp4
      if use_export_fig
         fr = export_fig('-nocrop','-r500'); % can use getframe(gcf) instead, but export_fig offers more control over resolution, etc.  
      else
         fr = getframe(gcf); 
      end
      writeVideo(v,fr); 
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

disp done