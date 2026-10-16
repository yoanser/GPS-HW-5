clear;
clc;
close all;

c = 299792458; %m/s , speed of light
NIST_ECEF = [-1288398.567 -4721696.932	4078625.350]; %Nist ecef location

%% Question 1 Analayzing Multiple PRN Satellites

PRN = [3; 4; 7; 8; 9; 16; 26; 27; 31];

%Read in rinex file
fileinfo = rinexinfo('NIST00USA_R_20262310000_01D_30S_MO.rnx');
GPS_Data = rinexread('NIST00USA_R_20262310000_01D_30S_MO.rnx');
%Read ephemeris values

[gps_ephem,ionoparams] = read_clean_GPSbroadcast('brdc2310.26n');

%% Question 3
ECEFGuess = lla2ecef([90 0 0])
straightLineECEF = NIST_ECEF - ECEFGuess
guessStraightLineError = norm(ECEFGuess - NIST_ECEF) %Straght line error in meters

i=0; %Init counter

while true

i = i+1



%Create seperate data files for each satellite and do all the calculations!
for k = 1:length(PRN)
    PRNUM = sprintf('PRN_%d', PRN(k));
    dataPRN.(PRNUM) = GPS_Data.GPS(GPS_Data.GPS.SatelliteID == PRN(k), :); %Filter for ith PRN data
    [wn.(PRNUM),tow.(PRNUM) ,~] = cal2gps(dataPRN.(PRNUM).Time);
    timestep = seconds(30);

    %Filtering, fill in spots that satellites isnt visible with NANS to make
    %plots more accurate
    %dataPRN.(dynamicName) = retime(dataPRN.(dynamicName),'regular', 'fillwithmissing', 'TimeStep', timestep);

    %Find iono range at first epoch
    PIF3b(k,:) = ionocorr(dataPRN.(PRNUM).C1C(1), 1575.42e6, dataPRN.(PRNUM).C2W(1), 1227.60e6);

    %Find expected range
    expectedRange3b(k,:) = findExpectedRange(ECEFGuess(i,:),tow.(PRNUM)(1),wn.(PRNUM)(1),gps_ephem,PRN(k));

    %Find expected elevation and azmiuth
    %read ephemeris data for each satellite and get clock correction!
    [~,satPos3b.(PRNUM),~,bsv3b(k,:),relCorr3b(k,:),~,~] = eph2pvt2025(gps_ephem,[wn.(PRNUM)(1) tow.(PRNUM)(1)],PRN(k));

    %Find Azmiuth and Elevation
    [AZ3b(k,:), EL3b(k,:), ~] = compute_azelrange(ECEFGuess(i,:),satPos3b.(PRNUM));


    %Find tropospheric correction (Assuming Zd = 2m again!)
    tropo3b(k,:) = tropomodel(2,EL3b(k,:));


end

resultsTable3b = table(PRN, PIF3b, expectedRange3b, AZ3b, EL3b, bsv3b, relCorr3b, tropo3b, ...
    'VariableNames', {'PRN','PIF (m)','ExpectedRange (m)','AZ (deg)','Elevation (deg)','BSV (m)','RELAT (m)','TROP (m)'});
disp(resultsTable3b);

prefitResid3b = PIF3b - (expectedRange3b - bsv3b - relCorr3b + tropo3b)

for j = 1:length(PRN)
    PRNUM = sprintf('PRN_%d', PRN(j));
    %Find G
    G3b(j,:) = [-(satPos3b.(PRNUM)(1) - ECEFGuess(i,1))./ expectedRange3b(j), -(satPos3b.(PRNUM)(2) - ECEFGuess(i,2))./ expectedRange3b(j),   ...
        -(satPos3b.(PRNUM)(3) - ECEFGuess(i,3))./ expectedRange3b(j), 1];

end


%Error 
solutionQuestion3c(:,i) = inv(G3b'*G3b)*G3b'*prefitResid3b


postfitResid = prefitResid3b - G3b*solutionQuestion3c(:,i)

ECEFGuess(i+1,:) = ECEFGuess(i,:) + (solutionQuestion3c(1:3,i))'

newError(i) = norm(ECEFGuess(i+1,:)-NIST_ECEF)
normofSolution(i) = norm(solutionQuestion3c(1:3,i))

if norm(solutionQuestion3c(1:3,i))<0.01 %Breaks while loop if error state vector is less than 1m
   
    break

end



end

finalErrorVec = NIST_ECEF - ECEFGuess(end,:) 
finalClockbias = solutionQuestion3c(end,end)
finalPostfitresid = postfitResid



NISTLLA = ecef2lla(NIST_ECEF);
R_ECEFtoENU = ECEF2ENU(NISTLLA(1),NISTLLA(2));
finalErrorVecENU = R_ECEFtoENU * finalErrorVec'

H_ECEF = inv(G3b'*G3b);
 
H_ENUSpatial = R_ECEFtoENU*H_ECEF(1:3,1:3)*R_ECEFtoENU';

HDOP = sqrt(H_ENUSpatial(1,1) + H_ENUSpatial(2,2))
VDOP = sqrt(H_ENUSpatial(3,3))