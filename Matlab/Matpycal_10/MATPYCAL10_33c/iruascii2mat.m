% iruascii2mat.m
% convert iru ascii data (getirudata33b) to mat file
% add option to turn off aberration adjustment
% [message] = iruascii2mat(filename,aber_adjust)

function [message] = iruascii2mat(filename,aber_adjust)

%initialization
  sol = 299792.458; % speed of light (km/sec)
  message = {'iruascii2mat'};
  
% read data
  data = load(filename,'-ascii');
  numrows = size(data,1);
  numcols = size(data,2);
  mess = sprintf('read file "%s", %d rows, %d columns',filename,numrows,numcols);
  disp(mess)
  message = [message {mess}];
  if aber_adjust
      mess = 'Adjust for aberration';
  else
      mess = 'No aberration adjustment';
  end
  disp(mess)
  message = [message {mess}];
  
% put data in arrays
  num = data(:,1);
% put interval start time into date vector
  start_datevec = zeros(numrows,6); % preallocate MATLAB date vector
  start_datevec(:,1) = data(:,2); % get year => [year 0 0 0 0 0]
  start_datevec(:,2) = 1; % set mon to 1 => [year 1 0 0 0 0]
  start_datevec(:,3) = data(:,3); % get doy => [year 1 doy 0 0 0]
  start_datenum = datenum(start_datevec); % [year 1 doy 0 0 0] => datenum
  start_datevec = datevec(start_datenum); % datenum => [year mon day 0 0 0]
  start_datevec(:,4:6) = data(:,4:6); % get hrs min sec => [year mon day hrs min sec]
  start_datenum = datenum(start_datevec); % [year 1 doy hrs min sec] => datenum
  start_J2000sec = (start_datenum - datenum([2000 1 1 12 0 0]))*86400.0;
  J2000sec_20060101 = (datenum([2006 1 1 0 0 0]) - datenum([2000 1 1 12 0 0]))*86400.0; 
  J2000sec_20090101 = (datenum([2009 1 1 0 0 0]) - datenum([2000 1 1 12 0 0]))*86400.0; 
  J2000sec_20120701 = (datenum([2012 7 1 0 0 0]) - datenum([2000 1 1 12 0 0]))*86400.0; 
  idx = (start_J2000sec >= J2000sec_20060101); % 2005-12-31 leap second
  start_J2000sec(idx) = start_J2000sec(idx) + 1.0;
  idx = (start_J2000sec >= J2000sec_20090101); % 2009-12-31 leap second
  start_J2000sec(idx) = start_J2000sec(idx) + 1.0;
  idx = (start_J2000sec >= J2000sec_20120701); % 2012-06-30 leap second
  start_J2000sec(idx) = start_J2000sec(idx) + 1.0;
% put interval start time into date vector
  stop_datevec = ones(numrows,6); % preallocate MATLAB date vector
  stop_datevec(:,1) = data(:,7); % get year => [ year 0 0 0 0 0 ]
  stop_datevec(:,2) = 1; % set mon to 1 => [ year 1 0 0 0 0 ]
  stop_datevec(:,3) = data(:,8); % get doy => [year 1 doy 0 0 0 ]
  stop_datenum = datenum(stop_datevec); % [year 1 doy 0 0 0 ] => datenum
  stop_datevec = datevec(stop_datenum); % datenum => [year mon day 0 0 0 ]
  stop_datevec(:,4:6) = data(:,9:11); % get hrs min sec => [year mon day hrs min sec ]
  stop_datenum = datenum(stop_datevec); % [year 1 doy hrs min sec] => datenum
  stop_J2000sec = (stop_datenum - datenum([2000 1 1 12 0 0]))*86400.0;
  idx = (stop_J2000sec >= J2000sec_20060101); % end of 2005 leap second
  stop_J2000sec(idx) = stop_J2000sec(idx) + 1.0;
  idx = (stop_J2000sec >= J2000sec_20090101); % end of 2009 leap second
  stop_J2000sec(idx) = stop_J2000sec(idx) + 1.0;
% get other data
  initquat = data(:,12:15)'; % initial pcad quaternion
  finalquat = data(:,16:19)'; % final pcad quaternion
  manvrquat = data(:,20:23)'; % maneuver quat from propagated iru rates
  intratebody = data(:,24:26)'; % integrated body rates (radians)
  diffchancnt = data(:,27:30)'; % difference in channel counts
  ini2finquat = data(:,31:34)'; % rotation quat from pcad initial to final quat
  ini2finvect = data(:,35:37)'; % rotation vector from ini2finquat
  finpropquat = data(:,38:41)'; % final rate propagated quaternion
  deltaquat   = data(:,42:45)'; % delta quaternion at end of interval
  ini2fintim  = data(:,46)'; % duration of interval (sec)
  ini2finang  = data(:,47)'; % maneuver angle from pcad quat (deg)
  sumprop = zeros(numrows,3,3); % sum of propagation matrices
  sumprop(:,1,1:3) = data(:,48:50);
  sumprop(:,2,1:3) = data(:,51:53);
  sumprop(:,3,1:3) = data(:,54:56);
  sumproprot = zeros(numrows,3,9); % sum of propagated rotations
  sumproprot(:,1,1:9) = data(:,57:65);
  sumproprot(:,2,1:9) = data(:,66:74);
  sumproprot(:,3,1:9) = data(:,75:83);
  pcadbias = data(:,84:86)'; % pcad bias at start of maneuver (rad/sec)
  cntratebiasbef = data(:,87:90)'; % count rate bias before man (cnt/sec)
  cntratebiasaft = data(:,91:94)'; % count rate bias after man (cnt/sec)
  clear data
  
% normalize all quaternions
  initquat = quatnorm(initquat);
  finalquat = quatnorm(finalquat);
  manvrquat = quatnorm(manvrquat);
  ini2finquat = quatnorm(ini2finquat);
  finpropquat = quatnorm(finpropquat);
  deltaquat = quatnorm(deltaquat);
  
% check data consistency
% initquat, finalquat, and ini2finquat
  testquat = quatmult(quatconj(initquat),finalquat);
  diffquat = quatmult(quatconj(testquat),ini2finquat);
  diffmag = sqrt(dot(diffquat(1:3,:),diffquat(1:3,:)));
  mess = sprintf('max testquat to ini2finquat diff = %6e',max(diffmag));
  disp(mess);
  message = [message {mess}];

% ini2finquat and ini2finvect
  testvect = quat2vect(ini2finquat);
  diffvect = testvect - ini2finvect;
  diffmag = sqrt(dot(diffvect,diffvect));
  mess = sprintf('max testvect to ini2finvect diff = %6e',max(diffmag));
  disp(mess);
  message = [message {mess}];

% initquat,manvrquat,& finpropquat
  testquat = quatmult(initquat,manvrquat);
  diffquat = quatmult(quatconj(testquat),finpropquat);
  diffmag = sqrt(dot(diffquat(1:3,:),diffquat(1:3,:)));
  mess = sprintf('max testquat to finpropquat diff = %6e',max(diffmag));
  disp(mess);
  message = [message {mess}];

% finalquat, finpropquat,& deltaquat
  diffquat = quatmult(quatconj(finpropquat),finalquat);
  diffquat = quatmult(quatconj(diffquat),deltaquat);
  diffmag = sqrt(dot(diffquat(1:3,:),diffquat(1:3,:)));
  mess = sprintf('max testquat to deltaquat diff = %6e',max(diffmag));
  disp(mess);
  message = [message {mess}];

% maneuver duration
  difftime = abs((stop_J2000sec - start_J2000sec) - ini2fintim');
  mess = sprintf('max mantime diff = %6e',max(difftime));
  disp(mess);
  message = [message {mess}];
  
% Adjust initial and final PCAD attitudes for stellar aberration
  if aber_adjust
      for n = 1:numrows
          [ ~, Vinit ] = EarthPosVel(start_J2000sec(n)-43200.0);
          Xinit = QuatXAxis(initquat(:,n));
          [ ~, Vfinal ] = EarthPosVel(stop_J2000sec(n)-43200.0);
          Xfinal = QuatXAxis(finalquat(:,n));
          initquat(:,n) =  quatmult(vec2quat(cross(Xinit, Vinit/sol )), initquat(:,n));
          finalquat(:,n) = quatmult(vec2quat(cross(Xfinal,Vfinal/sol)),finalquat(:,n));
      end
  end
    
% save mat-file with subset of information
  [~, name] = fileparts(filename);
  matfilename = [name '.mat'];
  save(matfilename,'start_datevec','stop_datevec','start_J2000sec',...
      'stop_J2000sec','initquat','finalquat','manvrquat','intratebody',...
      'diffchancnt','ini2fintim','ini2finang','sumprop','sumproprot',...
      'pcadbias','cntratebiasbef','cntratebiasaft');
  mess = ['Output saved to ' matfilename];
  disp(mess);
  message = [message {mess}];
  sumfilename = [name '_mat.sum'];
  fid = fopen(sumfilename,'w');
  fprintf(fid,'%s\n',message{:});
  fclose(fid);
  
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

% returns the inertial X-Axis unit 3-vector of the attitude represented by q
% input  q : double(4,num); array of quaternions
% output X : double(3,num); array of 3-vectors
% XAxis(1) = q(1)^2 - q(2)^2 - q(3)^2 + q(4)^2;
% XAxis(2) = 2*(q(1)*q(2) + q(3)*q(4));  
% XAxis(3) = 2*(q(1)*q(3) - q(2)*q(4));
function X = QuatXAxis(q)
  X(1,:) = q(1,:).^2 - q(2,:).^2 - q(3,:).^2 + q(4,:).^2;
  X(2,:) = 2*(q(1,:).*q(2,:) + q(3,:).*q(4,:));
  X(3,:) = 2*(q(1,:).*q(3,:) - q(2,:).*q(4,:));
return

% vec2quat
% returns quaternion of a rotation vector 
function q = vec2quat(v)
angle = sqrt(v(1)*v(1) + v(2)*v(2) + v(3)*v(3));
  if angle == 0.0
     q = [0; 0; 0; 1];
  else
     sina2 = sin(angle/2.0);
     q = [ v(1)/angle*sina2
           v(2)/angle*sina2
           v(3)/angle*sina2
           cos(angle/2.0)  ];
     if q(4) < 0.0
       q = -q;
     end
  end
  
% EarthPosVel
% Earth position and velocity
% input time as seconds from 2000:001:00:00:00
function [Pos, Vel] = EarthPosVel(time)
deg = pi/180;
T = (time - 43200)/(86400*36525); % time from J2000 in cneturies
aop = 102.9373*deg + 0.324284*deg*T;
ecc = 0.01670862 - 0.00004204*T;
lan = 100.4658*deg + 35999.372851*deg*T;
sma = 1.00001161*149597870;
inc = 23.439291*deg - 0.0130042*deg*T;
man = mod(lan - aop,2*pi);
ean = man + ecc*sin(man) + 0.5*ecc^2*sin(2*man);
GM = 1.32712438E+11;
sininc = sin(inc);
cosinc = cos(inc);
sinaop = sin(aop);
cosaop = cos(aop);
S = [ 1     0      0
      0  cosinc -sininc
      0  sininc  cosinc];
S = S*[ cosaop  -sinaop   0
        sinaop   cosaop   0
           0        0     1 ];
sinean = sin(ean);
cosean = cos(ean);
x = sqrt(1-ecc*ecc);
Pos = S*[ cosean - ecc
          x*sinean
          0  ]*sma;
Vel = S*[ -sinean
          x*cosean  
          0  ]*sqrt(GM/sma)/(1-ecc*cosean);       

  

