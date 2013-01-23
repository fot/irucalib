% matpycal.m
% MATLAB function to read mat-file with IRU calibration data
% and compute updated 9-parameter calibration with forward/backward
% Kalman Filter.
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
%      matflag : logical; true compute relative to on-board M-matrix
% Output IRU calibration
% V01 wsdavis add more information to plots
% V02 wsdavis add relative to on-board calib
% V03 wsdavis add fourth calibration uplink
% V04 wsdavis enhance filters
% V05 wsdavis add 3-SF consistency check, change channel alignment plots
%             add angle plots
% V06 wsdavis use varargin, prompt for arguments
%             remove x and Px outputs
% V07 wsdavis 
% V08 wsdavis add option for automatic output of fig & png files
% V09 wsdavis add m-matrix uplink
% V10 wsdavis 
% [message] = matpycal10(irumatfile,plotopt,matflag)

function [message] = matpycal10(varargin)

% initialization
  version = '10';
  message = {['matpycal version ' version ', run date ' datestr(now,31)]};
  disp(message)
% input iru file name
  if ((nargin < 1) || not(exist(varargin{1},'file')))
      FilterSpec = {'*.mat','Binary mat file (*.mat)';...
                    '*.*',  'All Files (*.*)' };
      [FileName,PathName] = getfile(FilterSpec,'Select IRU mat File');
      if ischar(FileName)
          irumatfile = [PathName FileName];
      else
          mess = 'No mat file name, matpycal run aborted';
          disp(mess)
          message = [message {mess}];
          return
      end
  else
      irumatfile = varargin{1};
  end
  message = [message {['irumatfile: "' irumatfile '"']}];
% input plot option
  if ((nargin < 2) || (not(islogical(varargin{2}))))
      button = questdlg('Plot results?','Plot Option');
      if strcmp(button,'Yes')
          plotopt = true;
      elseif strcmp(button,'No')
          plotopt = false;
      else
          disp('Run cancelled')
          return
      end
  else
      plotopt = varargin{2};
  end
  if plotopt 
      button = questdlg('Save Plots?','Plot Save','Yes','No','Yes');
      if strcmp(button,'Yes')
          saveopt = true;
      else
          saveopt = false;
      end
  end
% input adjust for M-file option
  if ((nargin < 3) || not(islogical(varargin{3})))
      button = questdlg('Adjust for M-matrix?','M-matrix Adjust Option','No');
      if strcmp(button,'Yes')
          matflag = true;
      elseif strcmp(button,'No')
          matflag = false;
      else
          disp('Run cancelled')
          return
      end
  else
      matflag = varargin{3};
  end
  rad = 1;
  deg = pi/180; % deg to rad conversion, 180*deg = pi (radians)
  asec = deg/3600; % arcsec to rad conversion
  sec = 1;
  hr = 3600*sec; % hour to second conversion 
  day = 24*hr;
  attstd = [ 10 0.1 0.1 ]*asec; % initial attitude standard deviation
  Patt = diag(attstd).^2;
  sig_d = 2.5*asec;          % digital uncertainty per axis 
  sig_v = 4.3E-4*deg/hr^0.5; % gyro white noise
  sig_b = 1.0E-7*deg/sec;    % uncertainty of initial bias
  sig_u = 9.0E-4*deg/hr^1.5; % bias white noise
  Q1p = (0.00003*rad)^2/(100*day);
  Qcoef = [0 Q1p 0 0 ];
  minang = 30.0; % minimum maneuver angle for processing
  P0 = 1.0E-9;
  Px0 = eye(9,9)*P0; % initial state vector covariance
  x0 = [ 3.32*10^-6; -1.39*10^-4; -3.90*10^-5;...
         9.32*10^-5; -1.38*10^-5;  4.58*10^-6;...
         1.38*10^-5;  9.83*10^-6;  3.37*10^-6]; % initial state vector
  mess = sprintf('Q1p = %10.3e, P0 = %10.3e',Q1p,P0);
  message = [message {mess}];
% Uplinked M-matrix calibrations
  nummat = 5;

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

% Matrix Patch Request PR-309, CAP 1227, IRU-2 M-matrix (uplinked 2012-Mar-02 15:26)
% YrMat(5) = 2012 + (062 + 15/24 + 26/1440 - 1)/366.0;
  MmatJ2000sec(5) = (datenum([2012 1 62 15 26 0]) - datenum([2000 1 1 12 0 0]))*86400.0;
  Mmat{5} = [  2.484911e-04  1.933052e-04  5.790450e-05 
               2.882515e-05  2.505313e-04  4.593649e-05 
               4.250966e-05  8.943458e-06  2.930867e-04  ];

  MmatJ2000sec(6) = (datenum([2030 1 1 0 0 0]) - datenum([2000 1 1 12 0 0]))*86400.0;

% alpha filter init values
  X0 = [-1.0; 2.8; -2.2]*asec;
  PX0 = eye(3)* (pi/100)^2;
  
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
  
% Kalman Filtering  
  num = size(start_datevec,1);
  x = zeros(9,num); % preallocate state vector
  Px = zeros(9,9,num); % preallocate state vector covariance
  manerr = zeros(3,num); % preallocate resid manerr vector
% process first maneuver separately without process noise
  x(:,1) = x0; % initialize state vector x
  Px(:,:,1) = Px0; % initialize covariance
  H = reshape(sumproprot(1,:,:),3,9);
  dur = ini2fintim(1,1); % duration of maneuver
  manvrmat = qtomat(manvrquat(:,1)); % maneuver rotation matrix from init to final
  R = eye(3)*Rcov(sig_d,sig_v,sig_b,sig_u,dur); % iru noise during man
  R = R + Patt + manvrmat*Patt*manvrmat'; % add init and final att covar
  Kg = Px(:,:,1)*H'/(H*Px(:,:,1)*H' + R); % Kalman gain
  rotatequat = quatmult(quatconj(initquat(:,1)),finalquat(:,1)); % rotate quat from pcad att
  zee = quat2vect(quatmult(quatconj(manvrquat(:,1)),rotatequat)); % delta vector
  x(:,1) = x(:,1) + Kg*(zee - H*x(:,1)); % update state
  Px(:,:,1) = (eye(9) - Kg*H)*Px(:,:,1); % update covariance
  manerr(:,1) = zee - H*x(:,1); % compute residual error
% Continue processing forward filter
  for n = 2:num % forward filter
      x(:,n) = x(:,n-1);
      Qx = Qprocess(Qcoef,stop_J2000sec(n-1,1),start_J2000sec(n,1),Px0);
      Px(:,:,n) = Px(:,:,n-1) + Qx; % covariance propagation
      H = reshape(sumproprot(n,:,:),3,9);
      dur = ini2fintim(1,n); % duration of maneuver
      manvrmat = qtomat(manvrquat(:,n)); % maneuver rotation matrix from init to final
      R = eye(3)*Rcov(sig_d,sig_v,sig_b,sig_u,dur); % iru noise during man
      R = R + Patt + manvrmat*Patt*manvrmat'; % add init and final att covar
      Kg = Px(:,:,n)*H'/(H*Px(:,:,n)*H' + R); % Kalman gain
      rotatequat = quatmult(quatconj(initquat(:,n)),finalquat(:,n)); % rotate quat from pcad att
      zee = quat2vect(quatmult(quatconj(manvrquat(:,n)),rotatequat)); % delta vector
      x(:,n) = x(:,n) + Kg*(zee - H*x(:,n)); % update state
      Px(:,:,n) = (eye(9) - Kg*H)*Px(:,:,n); % update covariance
      manerr(:,n) = zee - H*x(:,n); % compute residual error
  end

% backward filter
  xb = x(:,num); % initialize backward state
  Pxb = Px(:,:,num); % initialize backward covariance
% Pxb = Px0;
  for n = (num-1):(-1):1
      xf = x(:,n); % forward state at time tn
      Pxf = Px(:,:,n); % forward covariance at time tn
      Qx = Qprocess(Qcoef,stop_J2000sec(n,1),start_J2000sec(n+1,1),Px0);
      Pxb = Pxb + Qx; % propagate backwards cov
%     xb = xb; % propagate backward state
      Wf = Pxb/(Pxf + Pxb); % forward weight
      Wb = Pxf/(Pxf + Pxb); % backward weight
      x(:,n) = Wf*xf + Wb*xb; % combine forward and backward states
      Px(:,:,n) = Pxf/(Pxf + Pxb)*Pxb; % smoothed covariance      
      rotatequat = quatmult(quatconj(initquat(:,n)),finalquat(:,n)); % rotate quat from pcad att
      zee = quat2vect(quatmult(quatconj(manvrquat(:,n)),rotatequat)); % delta vector
      H = reshape(sumproprot(n,:,:),3,9);
      manerr(:,n) = zee - H*x(:,n); % compute residual error      
      dur = ini2fintim(1,n); % duration of maneuver
      R = eye(3)*Rcov(sig_d,sig_v,sig_b,sig_u,dur); % iru noise during man
      R = R + Patt + manvrmat*Patt*manvrmat'; % add init and final att covar
      Kg = Pxb*H'/(H*Pxb*H' + R); % Kalman gain
      xb = xb + Kg*(zee - H*xb); % update backward state
      Pxb = Pxb - Kg*H*Pxb; % update backward covariance
  end
  
% output results of Davenport F/B 9-parameter filter
  mess = 'Output last result of Davenport F/B 9-parameter filter';
  message = [message {mess}];
  mess = 'Mmat = '; message = [message {mess}];
  mess = sprintf('  %16.8e  %16.8e  %16.8e',x(1,end-1),x(2,end-1),x(3,end-1));
  message = [message {mess}];
  mess = sprintf('  %16.8e  %16.8e  %16.8e',x(4,end-1),x(5,end-1),x(6,end-1));
  message = [message {mess}];
  mess = sprintf('  %16.8e  %16.8e  %16.8e',x(7,end-1),x(8,end-1),x(9,end-1));
  message = [message {mess}];
  
% adjust relative to on-board M-matrix
  if matflag
     for n = 1:num
          for m = 1:nummat
              if ((MmatJ2000sec(m) <= stop_J2000sec(n)) && ...
                  (stop_J2000sec(n)) < MmatJ2000sec(m+1))
                  Mx = vectomat(x(:,n),3,3);
                  D = (eye(3) + Mx)/(eye(3) + Mmat{m}) - eye(3);
                  x(:,n) = mattovec(D);
              end
          end
     end
     mess = 'Dmat = '; message = [message {mess}];
     mess = sprintf('  %16.8e  %16.8e  %16.8e',x(1,end-1),x(2,end-1),x(3,end-1));
     message = [message {mess}];
     mess = sprintf('  %16.8e  %16.8e  %16.8e',x(4,end-1),x(5,end-1),x(6,end-1));
     message = [message {mess}];
     mess = sprintf('  %16.8e  %16.8e  %16.8e',x(7,end-1),x(8,end-1),x(9,end-1));
     message = [message {mess}]; 
  end
  
% Plots
%  years = 2003.564 + (stop_J2000sec - stop_J2000sec(1,1))/(86400*365.25);
  years = datevec2year(stop_datevec);
  if plotopt
%     plot scale factor adjustment, entire time span
      figure()
      plot(years,x(1,:)/asec*pi,'r-',years,x(5,:)/asec*pi,'g-',years,x(9,:)/asec*pi,'b-');%hold on;
      title('Scale-Factor Adjustment for X, Y, & Z Axes')
      xlabel('year'); ylabel('180-deg Man Err (arcsec)'); grid on
      legend('X-axis (Roll)','Y-axis (Pitch)','Z-axis (Yaw)','Location','Best')
      if saveopt
          saveas(gcf,'XYZaxesScaleFact.fig','fig');
          saveas(gcf,'XYZaxesScaleFact.png','png');
      end
%     plot residual y & z maneuver errors
      figure()
      yzerr = manerr(2,:).^2 + manerr(3,:).^2;
      rms = sqrt(mean(yzerr));
      yzerr = sqrt(yzerr);
      plot(years,yzerr/asec);
      title(['Resid YZ Maneuver errors, RMS = ' sprintf('%9.3f',rms/asec)]);
      xlabel('year'); ylabel('Man Err (arcsec)'); grid on
      mess = sprintf('RMS YZ resid error = %9.3f asec',rms/asec);
      message = [message {mess}];
      if saveopt
          saveas(gcf,'YZResidManErr.fig','fig');
          saveas(gcf,'YZResidManErr.png','png');
      end
%     plot scale factor adjustment stdev
      figure()
      Pxx = sqrt(reshape(Px(1,1,:),1,num)); 
      Pyy = sqrt(reshape(Px(5,5,:),1,num)); 
      Pzz = sqrt(reshape(Px(9,9,:),1,num)); 
      plot(years,Pxx/asec*pi,'r-',years,Pyy/asec*pi,'g-',years,Pzz/asec*pi,'b-');
      title('Scale-Factor Adjustment for X, Y, & Z Axes StDev')
      xlabel('year'); ylabel('180-deg Man Err (arcsec)'); grid on
      legend('X-axis (Roll)','Y-axis (Pitch)','Z-axis (Yaw)','Location','Best')
      if saveopt
          saveas(gcf,'XYZScaleFactStd.fig','fig');
          saveas(gcf,'XYZScaleFactStd.png','png');
      end
  end

  if matflag
      fprintf('%s\n',message{:})
      return
  end
  
% Compute 12-parameter calibration ------------------------------------
% Renormalize G
  Umat = Gmat'/(Gmat*Gmat');
  Umat = diag(sqrt(dot(Umat,Umat,2)))\Umat; % normalize rows
  Gmat = (Umat'*Umat)\Umat';

% Compute null vector
  Qmat = eye(4) - Umat*Gmat; % Projection into null space
  Hmat = Qmat/diag(sqrt(dot(Qmat,Qmat))); % normalize columns
  Hvec = Hmat(:,1); % null vector
  
% initialize alpha filter
  X  = zeros(3,num); % preallocate state vector, alpha
  PX = zeros(3,3,num); % preallocate state vector covariance
  hiderr = zeros(1,num); % preallocate resid hidden error vector
  Z = zeros(1,num); % preallocate Z

% forward filter
  X(:,1) = X0; % initialize state vector X
  PX(:,:,n) = PX0; % initialize covariance
  for n = 2:num 
      X(:,n) = X(:,n-1);
      QX = Qprocess(Qcoef,stop_J2000sec(n-1,1),start_J2000sec(n,1),PX0);
      PX(:,:,n) = PX(:,:,n-1) + QX;
%     4-vector angle difference
      dur = ini2fintim(1,n); % duration of maneuver
      DeltaCnts = diffchancnt(:,n) - dur*(cntratebiasbef(:,n) + cntratebiasaft(:,n))/2;
      poscnt = (DeltaCnts >= 0);
      negcnt = (DeltaCnts < 0);
      DeltaTheta = DeltaCnts.*(poscnt.*kpos + negcnt.*kneg);
%     observation Z
      Z(n) = dot(Hvec,DeltaTheta);
      R = (sig_d^2 + sig_v^2*dur + sig_b^2*dur^2 + sig_u^2*dur^3); % iru noise during man
      H = ((eye(3) + vectomat(x(:,n),3,3))*intratebody(:,n))';
      Kg = PX(:,:,n)*H'/(H*PX(:,:,n)*H' + R); % Kalman gain
      X(:,n) = X(:,n) + Kg*(Z(n) - H*X(:,n)); % update state
      PX(:,:,n) = PX(:,:,n) - Kg*H*PX(:,:,n); % update covariance
      hiderr(n) = Z(n) - H*X(:,n); % compute residual error
  end

% backward filter
  Xb = X(:,num); % initialize backward state
  PXb = PX(:,:,num); % initialize backward covariance
% PXb = PX0;
  for n = (num-1):(-1):1
      Xf = X(:,n); % forward state
      PXf = PX(:,:,n); % forward covariance
      Xb = X(:,n+1); % backward state not yet updated
      QX = Qprocess(Qcoef,stop_J2000sec(n,1),start_J2000sec(n+1,1),PX0);
      PXb = PXb + QX;
      dur = ini2fintim(1,n); % duration of maneuver
      R = Rcov(sig_d,sig_v,sig_b,sig_u,dur); % iru noise during man
      H = ((eye(3) + vectomat(x(:,n),3,3))*intratebody(:,n))';
      Kg = PXb*H'/(H*PXb*H' + R); % Kalman gain
      Xb = Xb + Kg*(Z(n) - H*Xb); % update backward state
      PXb = PXb - Kg*H*PXb; % update backward covariance
      PX(:,:,n) = eye(3)/(eye(3)/PXf + eye(3)/PXb); % combine forward and backward covariances
      X(:,n) = PX(:,:,n)*(PXf\Xf + PXb\Xb); % combine forward and backward states
      hiderr(n) = Z(n) - H*X(:,n); % compute residual error
  end
  
  if plotopt
%     plot alpha
      figure()
      plot(years,X(1,:)/asec,'r-',years,X(2,:)/asec,'g-',years,X(3,:)/asec,'b-');
      title('Alpha for X, Y, & Z Axes')
      xlabel('year'); ylabel('Alpha (arcsec)'); grid on
      legend('X-axis (Roll)','Y-axis (Pitch)','Z-axis (Yaw)','Location','Best')
      if saveopt
          saveas(gcf,'AlphaXYZ.fig','fig');
          saveas(gcf,'AlphaXYZ.png','png');
      end

%     plot residual hidden errors
      figure()
      rms = sqrt(mean(hiderr)^2 + std(hiderr)^2);
      plot(years,Z/asec,'r-',years,hiderr/asec,'b-');
      title(['Resid Hidden Errors, RMS = ' sprintf('%9.3f',rms/asec)]);
      xlabel('year'); ylabel('Hid Err (arcsec)'); grid on
      legend('Before Calib','AfterCalib','Location','Best')
      mess = sprintf('RMS hidden error = %9.3f asec',rms/asec);
      message = [message {mess}];
      if saveopt
          saveas(gcf,'AlphaResid.fig','fig');
          saveas(gcf,'AlphaResid.png','png');
      end
  end

% compute individual axis calibration
  kave = zeros(4,num);
  Alignvec = zeros(3,num);
  uvec1 = zeros(3,num);
  uvec2 = zeros(3,num);
  uvec3 = zeros(3,num);
  uvec4 = zeros(3,num);
  vecang12 = zeros(1,num);
  vecang13 = zeros(1,num);
  vecang14 = zeros(1,num);
  vecang23 = zeros(1,num);
  vecang24 = zeros(1,num);
  vecang34 = zeros(1,num);
  if not(matflag)
      for n = 1:num
          Mmat = vectomat(x(:,n),3,3);
          alpvec = X(:,n);
          SUmat = Umat/(eye(3) + Mmat) + Hvec*alpvec';
          [Smat Umatnew Alignmat] = SUParams(SUmat);
          uvec1(:,n) = Umatnew(1,:)';
          uvec2(:,n) = Umatnew(2,:)';
          uvec3(:,n) = Umatnew(3,:)';
          uvec4(:,n) = Umatnew(4,:)';
          if (n == 1), Alignmat0 = Alignmat; end
          Svec = diag(Smat);
          kave(:,n) = 0.5*(kpos + kneg)./Svec; 
          Alignvec(:,n) = vec(Alignmat*Alignmat0');
          vecang12(n) = vecang(Umatnew(1,:)',Umatnew(2,:)');
          vecang13(n) = vecang(Umatnew(1,:)',Umatnew(3,:)');
          vecang14(n) = vecang(Umatnew(1,:)',Umatnew(4,:)');
          vecang23(n) = vecang(Umatnew(2,:)',Umatnew(3,:)');
          vecang24(n) = vecang(Umatnew(2,:)',Umatnew(4,:)');
          vecang34(n) = vecang(Umatnew(3,:)',Umatnew(4,:)');
      end
%     plot residual scale factors for each axis
      figure()
      plot(years,(kave(1,:)-kave(1,1))/asec,'r-',...
           years,(kave(2,:)-kave(2,1))/asec,'g-',...
           years,(kave(3,:)-kave(3,1))/asec,'b-',...
           years,(kave(4,:)-kave(4,1))/asec,'m-');
      title('Scale Factor Change for Each Channel');
      xlabel('year'); ylabel('SF Change (arcsec/cnt)'); grid on
      legend('Channel-1','Channel-2','Channel-3','Channel-4','Location','Best');
      if saveopt
          saveas(gcf,'ChanSFchange.fig','fig');
          saveas(gcf,'ChanSFchange.png','png');
      end
  
%     plot overall alignment change
      figure()
      plot(years,Alignvec(1,:)/asec,'r-',...
           years,Alignvec(2,:)/asec,'g-',...
           years,Alignvec(3,:)/asec,'b-');
      title('IRU/ACA Relative Alignment Change');
      xlabel('year'); ylabel('Align Change (arcsec)'); grid on
      legend('Body-X','Body-Y','Body-Z','Location','Best');
      if saveopt
          saveas(gcf,'AlignChange.fig','fig');
          saveas(gcf,'AlignChange.png','png');
      end
    
%     plot XY-Interal Axis Change for Each Gyro
      figure()
      plot(years,(vecang12-vecang12(1))/asec,'r-',...
           years,(vecang34-vecang34(1))/asec,'g-');
      title('XY Interal Axis Change for Each Gyro');
      xlabel('year'); ylabel('XY Align Change (arcsec'); grid on
      legend('Gyro-1, Channels-12','Gyro-2, Channels-34','Location','Best');
      if saveopt
          saveas(gcf,'XYalignChange.fig','fig');
          saveas(gcf,'XYalignChange.png','png');
      end
  
%     plot relative axis change
      figure()
      plot(years,(vecang13-vecang13(1))/asec,'r-',...
           years,(vecang14-vecang14(1))/asec,'g-',...
           years,(vecang23-vecang23(1))/asec,'b-',...
           years,(vecang24-vecang24(1))/asec,'m-');
      title('Relative Axis Change between Gyros');
      xlabel('year'); ylabel('Relative Change (arcsec)'); grid on
      legend('Channels-13','Channels-14','Channels-23','Channels-24','Location','Best');
      if saveopt
          saveas(gcf,'GyroAlignChange.fig','fig');
          saveas(gcf,'GyroAlignChange.png','png');
      end

%     plot maneuver quaternion manvrquat(:,idx)
      figure()
      plot(years,manvrquat(1,:),'r-',...
           years,manvrquat(2,:),'g-',...
           years,manvrquat(3,:),'b-');
      title('q1-3 maneuver quat');
      xlabel('year'); ylabel('quat'); grid on
      legend('quat-1','quat-2','quat-3','Location','Best');
  end  
      if saveopt
          saveas(gcf,'Q123ManQuat.fig','fig');
          saveas(gcf,'Q123ManQuat.png','png');
      end
  
% Compute M = (G*S*U')^-1 - I for consistency check
  xn = zeros(3,num);
  for n = 1:num
      Svec = 0.5*(kpos + kneg)./kave(:,n);
      Umatnew = [uvec1(:,n),uvec2(:,n),uvec3(:,n),uvec4(:,n)]';
      M = inv(Gmat*diag(Svec)*Umatnew) - eye(3);
      xn(:,n) = diag(M);
  end

  if plotopt
%     plot 3-scale factor adjustment consistency check
      figure()
      plot(years,(xn(1,:)-x(1,:))/asec*pi,'r-',years,(xn(2,:)-x(5,:))/asec*pi,'g-',years,(xn(3,:)-x(9,:))/asec*pi,'b-');
      title([{'Scale-Factor Adjustment for X, Y, & Z Axes'} {'Consistency Check'}])
      xlabel('year'); ylabel('180-deg Man Err (arcsec)'); grid on
      legend('X-axis (Roll)','Y-axis (Pitch)','Z-axis (Yaw)','Location','Best')
      if saveopt
          saveas(gcf,'SFConsisChk.fig','fig');
          saveas(gcf,'SFConsisChk.png','png');
      end
  
  end
  
% finalization
% write messages to workspace
  fprintf('%s\n',message{:})
% write messages to file
  [~, name] = fileparts(irumatfile);
  sumfilename = ['matpycal_' version '_' name '.sum'];
  fsum = fopen(sumfilename,'w');
  fprintf(fsum,'%s\n',message{:});
  fclose(fsum);
   fprintf('Summary file written to %s\n', [pwd, filesep, sumfilename])

return

% function to compute noise stdev per IRU axis
% 
function Rn = Rcov(sig_d,sig_v,sig_b,sig_u,dur)
  Rn = sig_d^2 + sig_v^2*dur + sig_b^2*dur^2 + sig_u^2*dur^3;
return

% Qprocess
% function to compute process noise
% input:    Qcoef(1:dim) : coefficients of cubic process noise in time
%         stop_J2000sec: stop time of previous maneuver interval in sec from 2000:001:12:00:00 UTC
%        start_J2000sec: start time of next maneuver interval in sec from 2000:001:12:00:00 UTC
%          P0(1:num,1:num): covariance after reset
% output: Q(1:dim,1:dim) : numxnum process noise matrix
% 
function Q = Qprocess(Qcoef,stop_J2000sec,start_J2000sec,P0)
% Times of reset (e.g. IRU turned off) seconds from J2000 UTC (2000:001:12:00:00)
  dim = size(P0,1);
% IRU-2 made operational
  reset_J2000sec(1) = 585910800.0 - 473342400.0;
% rough time of mid-safe mode starting 2011:187:12:28:47 UTC, ending 2011:192:03:54:30
  reset_J2000sec(2) = 836697600.0 - 473342400.0; 
% rough time of mid-safe mode starting 2012:150:03:33:29 UTC, ending 2012:152:00:00:00
  reset_J2000sec(3) = 864950400.0 - 473342400.0;
  if any(and((stop_J2000sec < reset_J2000sec),(reset_J2000sec < start_J2000sec)))
      Q = P0; % reset to large uncertainty
  else
      dur = abs(start_J2000sec - stop_J2000sec);
      Q = (Qcoef(1) + Qcoef(2)*dur + Qcoef(3)*dur^2 + Qcoef(4)*dur^3);
      Q = Q*eye(dim);
  end
return

% function to convert nxm matrix to nmx1 vector, row ordered
function vec = mattovec(mat)
  [nrows ncols] = size(mat);
  vec = reshape(mat',nrows*ncols,1);
return

% function to convert nmx1 vector to nxm matrix column ordered
function mat = vectomat(vec,nrows,ncols)
  num = max(size(vec));
  if (num == nrows*ncols)
      mat = reshape(vec,ncols,nrows)';
  else
      mat = [];
  end
return

% SUParams
% compute derived parameters from CXO IRU SU calibration
% input SU(1:4,1:3)
function [S U A] = SUParams(SU)
  S = diag(sqrt(dot(SU,SU,2)));
  U = S\SU;
  u = U';
  x = sum(u,2);
  x = x/sqrt(dot(x,x));
  y =(eye(3) - x*x')*u;
  y = y/diag(sqrt(dot(y,y)));
  R = x*x'-[0,-x(3),x(2);x(3),0,-x(1);-x(2),x(1),0];
  y = [y(:,1), R^2*y(:,2), R^3*y(:,3), R*y(:,4)];
  y = sum(y,2);
  y = y/sqrt(dot(y,y));
  z = cross(x,y);
  A = [x,y,z];
return

% qtomat.m
% function to convert quaternion to rotation matrix
% mat = qtomat(q)
function mat = qtomat(q)

  q1 = q(1); q2 = q(2); q3 = q(3); q4 = q(4);

  mat = [ q1^2-q2^2-q3^2+q4^2  2*(q1*q2+q3*q4)   2*(q1*q3-q2*q4)
          2*(q1*q2-q3*q4)  -q1^2+q2^2-q3^2+q4^2  2*(q2*q3+q1*q4)
          2*(q1*q3+q2*q4)      2*(q2*q3-q1*q4)   -q1^2-q2^2+q3^2+q4^2 ];

return

% quat2vect
% returns the rotation vector of a quaternion
% vect = quat2vect(quat)
function vect = quat2vect(quat)
  quat = quatnorm(quat);
  vect = quat(1:3,:);
  sinang = sqrt(dot(vect,vect)); % sine of rotation angle
  vect = 2.0*vect;
  idx = sinang > 1.0E-7; % true if cannot make small angle approximation
  angadj = atan2(sinang(1,:),quat(4,:))./sinang;
  vect(:,idx) = vect(:,idx).*repmat(angadj,3,1);
return

% quatconj
% returns conjugate of quaternion
% 
function q = quatconj(p)
  q = zeros(size(p));
  q(1:3,:) = -p(1:3,:);
  q(4,:) =  p(4,:);
return

% quatmult
% returns product of two quaterions                                       
% q = quatmult(r,p)
function q = quatmult(r,p)
  q = zeros(4,max(size(r,2),size(p,2)));
  q(1,:) = + r(1,:).*p(4,:) + r(2,:).*p(3,:) - r(3,:).*p(2,:) + r(4,:).*p(1,:);
  q(2,:) = - r(1,:).*p(3,:) + r(2,:).*p(4,:) + r(3,:).*p(1,:) + r(4,:).*p(2,:);
  q(3,:) = + r(1,:).*p(2,:) - r(2,:).*p(1,:) + r(3,:).*p(4,:) + r(4,:).*p(3,:);
  q(4,:) = - r(1,:).*p(1,:) - r(2,:).*p(2,:) - r(3,:).*p(3,:) + r(4,:).*p(4,:);
return

% quatnorm.m
% normalize to unit length and adjust for sign (q4 >= 0)
function q = quatnorm(p)
  mag = sqrt(dot(p,p));
  mag = repmat(mag,4,1);
  q = p./mag;
  idx = q(4,:) < 0;
  q(:,idx) = -q(:,idx);
return

function ang = vecang(u,v)
  x = u'*v;
  y = cross(u,v);
  ang = atan2(sqrt(y'*y),x);
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

% VEC
% function to compute rotation vector from rotation matrix
% function V = vec(A)
% input A is 3x3 orthonormal matrix
% output V is 3x1 rotation vector in radians
function V = vec(A)
  sv = 0.5*[ A(2,3)-A(3,2)
             A(3,1)-A(1,3)
           A(1,2)-A(2,1)];
  s = sqrt(sv'*sv);
  c = 0.5*(A(1,1) + A(2,2) + A(3,3) -1);
  a = atan2(s,c);
  if a < 0.00001
     V = sv;
  else
     V = sv/s*a;
  end
return

% getfile.m
% Expands uigetfile GUI to additional prompts and information
% title = 'Select File';
% FilterSpec = {'*.bin;*.BIN','Binary Files (*.bin;*.BIN)';...
%               '*.txt;*.TXT','Text Files (*.txt;*.TXT)';...
%               '*.*',  'All Files (*.*)' };
% FilterSpec = '*.*';
% function [FileName,PathName] = getfile(FilterSpec,DialogTitle,DefaultName)

function [FileName,PathName] = getfile(varargin)
  if ((nargin >= 1) && (ischar(varargin{1}) || iscell(varargin{1})))
      FilterSpec = varargin{1};
  else
      FilterSpec = '*.*';
  end
  if ((nargin >= 2) && (ischar(varargin{2})))
      DialogTitle = varargin{2};
  else
      DialogTitle = '';
  end
  if ((nargin >= 3) && (ischar(varargin{3})))
      DefaultName = varargin{3};
  else
      DefaultName = '';
  end
  done = false;
  while (not(done))
      [FileName,PathName] = uigetfile(FilterSpec,DialogTitle,DefaultName);
      if (not(ischar(FileName)))
          qstring = [{'No file was selected.'}...
                     {'Try Again?'}...
                     {'  Press "Yes" to select file'}...
                     {'  Press "No" to exit'}];
          button = questdlg(qstring,'No file Selected','Yes','No','Yes');
          if strcmp(button,'No')
              return
          end
      elseif (not(exist([PathName FileName],'file')))
          qstring = [{'File does not exist.'}...
                     {'Continue?'}...
                     {'  Press "Yes" to reselect'}...
                     {'  Press "No" to exit'}];
          button = questdlg(qstring,'title','Yes','No','Yes');
          if strcmp(button,'No')
              return
          end
      else
          qstring = [{['File "' FileName '" was selected.']}...
                     {'Use this file?'}...
                     {'  Press "Yes" to use file'}...
                     {'  Press "No" to reselect'}...
                     {'  Press "Cancel" to Exit'}];          
          button = questdlg(qstring,'Confirm File','Yes','No','Cancel','Yes');
          if strcmp(button,'Cancel')
              PathName = 0;
              FileName = 0;
              return
          elseif strcmp(button,'Yes')
              done = true;
          end
      end
  end
return
