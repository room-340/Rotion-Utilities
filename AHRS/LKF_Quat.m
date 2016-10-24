function [quat, angles, dw, P, Q, R, wPrime] = LKF_Quat(acc, mag, anr,...
    accl_coefs, magn_coefs, gyro_coefs,...
    dw, dT, quat, accBasic, magBasic,...
    P, Q, R)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%% Data correction
%     g = 9.81; % m/s^2
    accelerationThreshold = 0.3; % m/s^2
    fieldThreshold = 0;
    % Sensors data
    m = mag;
    a = acc;
    w = [anr(1) -anr(2) anr(3)];

    % Calibrate magnetometer
    B = [magn_coefs(1) magn_coefs(4) magn_coefs(5)
        magn_coefs(6) magn_coefs(2) magn_coefs(7)
        magn_coefs(8) magn_coefs(9) magn_coefs(3)]; % 3x3
    B0 = [magn_coefs(10); magn_coefs(11); magn_coefs(12)]; % 3x1
    mCal = ((eye(3,'single')-B)*(m'-B0))';

    % Calibrate accelerometer
    B = [accl_coefs(1) accl_coefs(4) accl_coefs(5)
        accl_coefs(6) accl_coefs(2) accl_coefs(7)
        accl_coefs(8) accl_coefs(9) accl_coefs(3)]; % 3x3
    B0 = [accl_coefs(10);accl_coefs(11);accl_coefs(12)]; % 3x1
    aCal = ((eye(3,'single')-B)*(a'-B0))';

    % Calibrate gyroscopes
    B = [gyro_coefs(1) gyro_coefs(4) gyro_coefs(5)
        gyro_coefs(6) gyro_coefs(2) gyro_coefs(7)
        gyro_coefs(8) gyro_coefs(9) gyro_coefs(3)]; % 3x3
    B0 = [gyro_coefs(10);gyro_coefs(11);gyro_coefs(12)]; % 3x1
    wCal = ((eye(3,'single')-B)*(w'-B0))';

    % Correct angular rate
    wPrime = wCal - dw;
    
%     absWPrime = sqrt(wPrime(1)^2 + wPrime(2)^2 + wPrime(3)^2);
    
    for i=1:3
        if abs(wPrime(i)) < 0.2
            wPrime(i) = 0;
        end;
    end;
%     if(absWPrime < 0.2)
%         wPrime = [0 0 0];
%     end;
    %% Calculate quaternion
    qPrime  = mrotate_eml(quat, wPrime, dT);
%     qPrime  = quatrotate_eml(quat, wPrime, dT);
%     qPrime = qPrime./norm(qPrime);
%     Fx = wPrime(1)*dT;
%     Fy = wPrime(2)*dT;
%     Fz = wPrime(3)*dT;
%     Fm = sqrt(Fx*Fx+Fy*Fy+Fz*Fz);
%     sinFm2 = sin(Fm/2);
%     cosFm2 = cos(Fm/2); 
%     if Fm~=single(0)
%         Nd=[cosFm2, Fx/Fm*sinFm2, Fy/Fm*sinFm2, Fz/Fm*sinFm2];
%     else
%         Nd=single([ 1, 0, 0, 0]);
%     end
%     qPrime = [Nd(1)*quat(1)-Nd(2:4)*quat(2:4)'...
%         Nd(1)*quat(2:4)+quat(1)*Nd(2:4)+cross(quat(2:4),Nd(2:4))];
    %% Kalman filter
    % Initial matrix values
    I = eye(6,'single');
    A = zeros(6,6,'single');
    Z = zeros(6,1,'single');
    H = zeros(6,6,'single');
    
    %System matrix
    F = [   1            wPrime(3)*dT  -wPrime(2)*dT -dT  0   0;
        -wPrime(3)*dT       1           wPrime(1)*dT  0  -dT  0;
         wPrime(2)*dT   -wPrime(1)*dT      1          0   0  -dT;
            0               0              0          1   0   0;
            0               0              0          0   1   0;
            0               0              0          0   0   1];


    % Estimated vector values 
%     accBasicPrime = quatrotate_eml( qPrime, accBasic);
%     magBasicPrime = quatrotate_eml( qPrime, magBasic);
    accBasicPrime = rotateVQ( accBasic, qPrime);
    magBasicPrime = rotateVQ( magBasic, qPrime);
%     qPrimeConj = [qPrime(1) -qPrime(2) qPrime(3) qPrime(4)];
%     qTemp = [-accBasic*qPrimeConj(2:4)' ...
%         +qPrimeConj(1)*accBasic+cross(qPrimeConj(2:4),accBasic)];
%     rQ = [qPrime(1)*qTemp(1)-qPrime(2:4)*qTemp(2:4)' ...
%         qPrime(1)*qTemp(2:4)+qTemp(1)*qPrime(2:4)+cross(qTemp(2:4),qPrime(2:4))];
%     accBasicPrime = rQ(2:4);
%     
%     qTemp = [-magBasic*qPrimeConj(2:4)' ...
%         +qPrimeConj(1)*magBasic+cross(qPrimeConj(2:4),magBasic)];
%     rQ = [qPrime(1)*qTemp(1)-qPrime(2:4)*qTemp(2:4)' ...
%         qPrime(1)*qTemp(2:4)+qTemp(1)*qPrime(2:4)+cross(qTemp(2:4),qPrime(2:4))];
%     magBasicPrime = rQ(2:4);
    
    % Chek if there is non-zero data to prevent divergence
    if (norm(aCal) > single(0)) && (norm(mCal) > single(0))
        % Measurement vector
        dA = aCal./norm(aCal) - accBasicPrime;
        dM = mCal - magBasicPrime;
        Z(1:3,1) = dA;
        Z(4:6,1) = dM;

        % Measurement matrix with gating
        threshold_a =  ...
            abs(sqrt(aCal(1)^2 + aCal(2)^2 + aCal(3)^2) - 9.81);
        if (threshold_a < accelerationThreshold)
            H(1,2) = -accBasicPrime(3); H(1,3) =  accBasicPrime(2);
            H(2,1) =  accBasicPrime(3); H(2,3) = -accBasicPrime(1);
            H(3,1) = -accBasicPrime(2); H(3,2) =  accBasicPrime(1);
        end

        threshold_m = abs(sqrt(mCal(1)^2 + mCal(2)^2 + mCal(3)^2) - single(1.0));
        if (threshold_m < fieldThreshold)
            H(4,2) = -magBasicPrime(3); H(4,3) =  magBasicPrime(2);
            H(5,1) =  magBasicPrime(3); H(5,3) = -magBasicPrime(1);
            H(6,1) = -magBasicPrime(2); H(6,2) =  magBasicPrime(1);
        end

        % Kalman Filter equations
        P = F*P*F' + Q;
        K = (P*H')/(H*P*H' + R);
        P = (I - K*H)*P*(I - K*H)' + K*R*K';
        X = K*Z;

        % Update Gyroscopes Bias estimation
        dw = dw + X(4:6)';

        % Correct quaternion
        qe  = single([X(1)/2 X(2)/2 X(3)/2]);
        qqe = [sqrt(single(1) - (qe(1)^2 + qe(2)^2 + qe(3)^2)),...
            qe(1), qe(2), qe(3)];
        quat = multQQ( qPrime, qqe );
%         quat = quat./norm(quat);
%         quat = [qqe(1)*qPrime(1)-qqe(2:4)*qPrime(2:4)' ...
%             qqe(1)*qPrime(2:4)+qPrime(1)*qqe(2:4)+cross(qPrime(2:4),qqe(2:4))];
        
        %% Angles
        % from Wiki
%         angles = [0 0 0];
        [angles(1), angles(2), angles(3)] = quat2angle_eml(quat);
%         %Normalize quaternion 
%         qin = quat./(sqrt(sum(quat.^2,2))* ones(1,4));
% 
%         %DCM elements
%         c21 = 2*(qin(2)*qin(3) + qin(1)*qin(4));
%         c11 = qin(1)^2 + qin(2)^2 - qin(3)^2 - qin(4)^2;
%         c31 = 2*(qin(2)*qin(4) - qin(1)*qin(3));
%         c32 = 2*(qin(3)*qin(4) + qin(1)*qin(2));
%         c33 = qin(1)^2 - qin(2)^2 - qin(3)^2 + qin(4)^2;
%         
%         angles = [0 0 0];
%         %Angles from DCM
%         angles(1) = atan2(c21,c11);
%         angles(2) = asin(-c31);
%         angles(3) = atan2(c32,c33);
    end
end

