% matpybias.m
% MATLAB function to read mat-file with IRU calibration data
% and compute plot bias terms
%
% Input 
%      mat-file with following quantities for each maneuver interval
%         start_datevec(num,6)
%         stop_datevec(num,6)
%         start_J2000sec(num,1)
%         stop_J2000sec(num,1)
%         initquat(4,num)
%         finalquat(4,num)
%         manvrquat(4,num)
%         intratebody(3,num)
%         diffchancnt(4,num)
%         ini2fintim(1,num)
%         ini2finang(1,num)
%         sumprop(num,3,3)
%         sumproprot(num,3,9)
%         pcadbias(3,num)
%         cntratebiasbef(4,num)
%         cntratebiasaft(4,num)
%      plotopt : logical; plot option (true or false)
% Output IRU calibration
% V01 wsdavis add more information to plots


function message = matpybias(irumatfile,plotopt)

% initialization
  version = '01';
  message = {['matpybias version ' version ', run date ' datestr(now,31)]};
  message = [message {['irumatfile: "' irumatfile '"']}];
  if plotopt
      message = [message {'plot option is true'}]; 
  else
      message = [message {'plot option is false'}]; 
  end
% rad = 1;
  deg = pi/180; % deg to rad conversion, 180*deg = pi (radians)
% asec = deg/3600; % arcsec to rad conversion
  hr = 3600;
  minang = 30.0; % minimum maneuver angle for processing
% Uplinked M-matrix calibrations
% nummat = 4;
% CAP 891, 09/19/2003, uplinked 2003:274:13:19:00
% YrMat(1) = 2003 + (274 + 13/24 + 19/1440 - 1)/365.0;
  MmatJ2000sec(1) = (datenum([2003 1 274 13 19 0]) - datenum([2000 1 1 12 0 0]))*86400.0;
  Mmat{1} = [ 3.3203451E-6 -1.3894030E-4 -3.8983014E-5
              9.3199606E-5 -1.3754384E-5  4.5790156E-6
              1.3764573E-5  9.8274797E-6  3.3727364E-6 ];
% CAP 1021, 12/13/2006, uplinked 2006:352:14:50:00
% YrMat(2) = 2006 + (352 + 14/24 + 50/1440 - 1)/365.0;
  MmatJ2000sec(2) = (datenum([2006 1 352 14 50 0]) - datenum([2000 1 1 12 0 0]))*86400.0;
  Mmat{2} = [ 0.392656e-04  0.920426e-04  0.047589e-04
             -1.509408e-04  0.429630e-04  0.155405e-04
             -0.477453e-04  0.168878e-04  0.679124e-04 ];
% Matrix Patch Request #: 283, uplinked 2010:350:22:10:00
% YrMat(3) = 2010 + (350 + 22/24 + 10/1440 - 1)/365.0;
  MmatJ2000sec(3) = (datenum([2010 1 350 22 0 0]) - datenum([2000 1 1 12 0 0]))*86400.0;
  Mmat{3} = [ 8.792448e-05  1.409469e-04  3.321078e-05
             -1.200405e-04  9.856613e-05  2.051482e-05
             -3.014042e-05  2.017329e-05  1.311966e-04 ];
% Matrix Patch Request #: 289, uplinked 2011:105:21:20:00
% YrMat(4) = 2011 + (104 + 21/24 + 20/1440 - 1)/365.0;
  MmatJ2000sec(4) = (datenum([2011 1 105 21 20 0]) - datenum([2000 1 1 12 0 0]))*86400.0;
  Mmat{4} = [  1.651433e-04  1.888956e-04  6.763121e-05
              -5.143320e-05  1.669320e-04  2.689127e-05
              -2.455693e-06  1.191769e-05  2.150693e-04 ];

  MmatJ2000sec(5) = (datenum([2030 1 1 0 0 0]) - datenum([2000 1 1 12 0 0]))*86400.0;
  
% on-board G-matrix:  CAP 0867, 2003:200
  Gmat = [ -0.499539493  0.500015266  0.500455729 -0.500504173
           -0.254059137  0.609733116 -0.253191860  0.610258254
           -0.557983976 -0.053139506 -0.556465488 -0.053843139 ];

% scale factors in radians/sec/minor-cycle, 1 minor-cycle = 1.025/16 sec in duration
% CAP 0867, 2003:200
  fpos = [0.1555267e-05; 0.1564191e-05; 0.1552225e-05; 0.1569751e-05];
  fneg = [0.1555571e-05; 0.1564383e-05; 0.1552371e-05; 0.1570025e-05];
  kpos = fpos*1.025/16; % radians/count, vector
  kneg = fneg*1.025/16; % radians/count, vector
  kave = (kpos + kneg)/2;
  
% read mat-file
  if  exist(irumatfile,'file')
      message = [message {['Loading file "' irumatfile '"' ]}];
      load(irumatfile);
      numman = size(start_datevec,1);
      message = [message {sprintf('Number of maneuvers = %d',numman)}];
      start_datestr = datestr(start_datevec(1,:),31);
      message = [message {['Data start: ' start_datestr]}];
      stop_datestr = datestr(stop_datevec(end,:),31);
      message = [message {['Data stop: ' stop_datestr]}];
  else
      mess = ['file: "' irumatfile '" does not exist'];
      message = [message {mess}];
      return
  end
  
% processing
% eliminate maneuvers less than minang
  message = [message {sprintf('Eliminate maneuvers less than %7.3f deg',minang)}];
  idx = (ini2finang >= minang);
  start_datevec = start_datevec(idx,:);
  stop_datevec = stop_datevec(idx,:);
  start_J2000sec = start_J2000sec(idx,:);
  stop_J2000sec = stop_J2000sec(idx,:);
  initquat = initquat(:,idx);
  finalquat = finalquat(:,idx);
  manvrquat = manvrquat(:,idx);
  intratebody = intratebody(:,idx);
  diffchancnt = diffchancnt(:,idx);
  ini2fintim = ini2fintim(:,idx);
  ini2finang = ini2finang(:,idx);
  sumprop = sumprop(idx,:,:);
  sumproprot = sumproprot(idx,:,:);
  pcadbias = pcadbias(:,idx);
  cntratebiasbef = cntratebiasbef(:,idx);
  cntratebiasaft = cntratebiasaft(:,idx);
  numman = size(start_datevec,1);
  message = [message {sprintf('Number of maneuvers = %d',numman)}];
  start_datestr = datestr(start_datevec(1,:),31);
  message = [message {['Data start: ' start_datestr]}];
  stop_datestr = datestr(stop_datevec(end,:),31);
  message = [message {['Data stop: ' stop_datestr]}];

% compute average 4-channel biases for each meneuver
  chan4bias = (cntratebiasbef + cntratebiasaft)/2;
  chan4bias(1,:) = kave(1)*chan4bias(1,:);
  chan4bias(2,:) = kave(2)*chan4bias(2,:);
  chan4bias(3,:) = kave(3)*chan4bias(3,:);
  chan4bias(4,:) = kave(4)*chan4bias(4,:);
  
% compute bias difference between before and after 
  diff4bias = cntratebiasbef - cntratebiasaft;
  diff4bias(1,:) = kave(1)*diff4bias(1,:);
  diff4bias(2,:) = kave(2)*diff4bias(2,:);
  diff4bias(3,:) = kave(3)*diff4bias(3,:);
  diff4bias(4,:) = kave(4)*diff4bias(4,:);
  
% compute 4-channel 3-biases
  chan3bias = Gmat*chan4bias;
  
% Plots
  years = datevec2year(stop_datevec);
  if plotopt
      figure()
      plot(years,chan4bias(1,:)*hr/deg,'r-',...
           years,chan4bias(2,:)*hr/deg,'g-',...
           years,chan4bias(3,:)*hr/deg,'b-',...
           years,chan4bias(4,:)*hr/deg,'m-');
      title('Bias for Each Channel of IRU-2');
      xlabel('year'); ylabel('bias (deg/hr)'); grid on
      legend('Channel-1','Channel-2','Channel-3','Channel-4','Location','Best');
      saveas(gcf,'4ChannelBias.fig','fig');
      saveas(gcf,'4ChannelBias.png','png');
      
      figure()
      plot(years,diff4bias(1,:)*hr/deg,'r-',...
           years,diff4bias(2,:)*hr/deg,'g-',...
           years,diff4bias(3,:)*hr/deg,'b-',...
           years,diff4bias(4,:)*hr/deg,'m-');
      title('Bias Diff for Each Maneuver for Each Channel of IRU-2');
      xlabel('year'); ylabel('bias diff (deg/hr)'); grid on
      legend('Channel-1','Channel-2','Channel-3','Channel-4','Location','Best');
      saveas(gcf,'4ChannelBiasDiff.fig','fig');
      saveas(gcf,'4ChannelBiasDiff.png','png');
      
      figure()
      plot(years,chan3bias(1,:)*hr/deg,'r-',...
           years,chan3bias(2,:)*hr/deg,'g-',...
           years,chan3bias(3,:)*hr/deg,'b-');
      title('3-Vector Bias from 4 Channels of IRU-2');
      xlabel('year'); ylabel('bias (deg/hr)'); grid on
      legend('X-axis (Roll)','Y-axis (Pitch)','Z-axis (Yaw)','Location','Best');
      saveas(gcf,'3VecBias4Channel.fig','fig');
      saveas(gcf,'3VecBias4Channel.png','png');
      
      figure()
      plot(years,pcadbias(1,:)*hr/deg,'r-',...
           years,pcadbias(2,:)*hr/deg,'g-',...
           years,pcadbias(3,:)*hr/deg,'b-');
      title('3-Vector Bias from PCAD for IRU-2');
      xlabel('year'); ylabel('bias (deg/hr)'); grid on
      legend('X-axis (Roll)','Y-axis (Pitch)','Z-axis (Yaw)','Location','Best');
      saveas(gcf,'3VecBiasPCAD.fig','fig');
      saveas(gcf,'3VecBiasPCAD.png','png');
            
  end
return

% datevec2year.m
% converts date vector to year and fraction of year
% input date vector of [year mon day hr min sec]
% output year
% year = datevec2year(timvec)
% uses built-in MATLAB functions which don't account for leap seconds
function year = datevec2year(timvec)
  year = timvec(:,1);
  doy = datenum(timvec(:,1:3)) - ...
        datenum([timvec(:,1) repmat([1 0],size(timvec,1),1)]) + ...
        (timvec(:,4)*3600 + timvec(:,5)*60 + timvec(:,6))/86400;
  year = year + (doy - 1)./(datenum(year+1,1,1) - datenum(year,1,1));
return