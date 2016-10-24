function [X, P, Q, R, debug] = LKF_Euler(acc, mag, anr,...
    accl_coefs, magn_coefs, gyro_coefs,...
    dT,... 
    X, P, Q, R, accThreshold, magThreshold, magDeclination)
%   Основной алгоритм
%   На вход идут текущие измерения и матрицы с предыдущей итерации
%% Data correction
    g = 9.81; % m/s^2
    debug = [0 0 0];
%   X = [Psi Gamma Thetta Wx Wy Wz dWx dWy dWz]
    % Sensors data
    m = mag;
    a = acc;
%     w = anr;
    % ------------------------------- CRUTCH
    w = [anr(1) -anr(2) anr(3)];
    gyro_coefs(11) = -gyro_coefs(11);
    % ---------------------------------------

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

    %% Prediction update
      F = [1 0 0 dT  0  0  0  0  0;
           0 1 0  0 dT  0  0  0  0;
           0 0 1  0  0 dT  0  0  0;
           0 0 0  1  0  0 -1  0  0;
           0 0 0  0  1  0  0 -1  0;
           0 0 0  0  0  1  0  0 -1;
           0 0 0  0  0  0  1  0  0;
           0 0 0  0  0  0  0  1  0;
           0 0 0  0  0  0  0  0  1];
      % system noise
%       gyroDrift = (30/3600)*pi/180; % deg/h -> deg/sec -> rad/sec
      Qv = [Q(1,1) Q(2,2) Q(3,3) Q(4,4) Q(5,5) Q(6,6) Q(7,7) Q(8,8) Q(9,9)]';
      X = F*X' + Qv;
      P = F*P*F' + Q;
      
    %% Measurement update
    % Initial matrix values
    I = eye(9,'single');
    Z = zeros(6,1,'single');
    H = zeros(6,9,'single');
    
%     H(1:3,1:3) = eye(3, 'single');
    H(4:6,4:6) = eye(3, 'single');
    
    
    % DCM basic axis -> object axis
%     Mbo = angles2dcm(X(1), X(2), X(3));
    
    
    % Estimated vector values 
%     accBasicPrime = Mbo*accBasic;
%     magBasicPrime = Mbo*magBasic;
    
    % Chek if there is non-zero data to prevent divergence
    if (norm(aCal) > single(0)) && (norm(mCal) > single(0))
        % Measurement vector
        Z(4:6) = wCal;
%         absANR = sqrt(wCal(1)^2 + wCal(2)^2 + wCal(3)^2);
%         if absANR < 0.1
%             Z(4:6) = 0;
%             debug(3) = 3;
%         end;
        for i=1:3
            if abs(wCal(i)) < 0.3
                Z(3+i) = 0;
                debug(3) = 3;
            end;
        end;
        H(1,1) = 1; H(2, 2) = 1; H(3, 3) = 1;
%         R(1,1) = 2e0; R(2,2) = 2e0; R(3,3) = 2e0;
        
        rollS = -atan2(aCal(2),sqrt(aCal(1)^2 + aCal(3)^2)); % крен
        pitchS = atan2(aCal(1),sqrt(aCal(2)^2 + aCal(3)^2)); % тангаж
        Moh = angle2dcm_eml(0, pitchS, rollS); % object -> horizon
        
        % Fill up Z vector for no-correction situations
        currentOrienthor = Moh*X(1:3);
        roll = currentOrienthor(1);
        pitch = currentOrienthor(2);
        yaw = currentOrienthor(3);
        
        % If object isnt moving - correct pitch & roll
        threshold_a = abs(sqrt(aCal(1)^2 + aCal(2)^2 + aCal(3)^2) - g);
        if (threshold_a < accThreshold)
%             R(1,1) = 2e-2; R(2,2) = 2e-2; R(3,3) = 2e-2;
            % Note that corrected roll angle has different sign (opposite to rollS)
            roll = atan2(aCal(2),sqrt(aCal(1)^2 + aCal(3)^2)); % крен
            % ----------------------------------------------------- CRUTCH
            pitch = -atan2(aCal(1),sqrt(aCal(2)^2 + aCal(3)^2)); % тангаж
            debug(1) = 1;
        end
        % If magnetic field is stable - correct yaw
        threshold_m = abs(sqrt(mCal(1)^2 + mCal(2)^2 + mCal(3)^2) - 1.0);
        if (threshold_m < magThreshold)
            magHeading = Moh*mCal';
            yaw = -atan2(magHeading(2),magHeading(1)) + magDeclination;
            debug(2) = 2;
        end
        
        orientHor = [roll pitch yaw];
        % Convert orientation angles from horizon to object axis
        orientObj = Moh'*orientHor'; % horizon -> object X Y Z
        Z(3) = orientObj(3);
        Z(2) = orientObj(2);
        Z(1) = orientObj(1);
        
        
        
%         debug = Z(1:3);
%         Z(1:3) = [-roll pitch yaw];
%         debug = [roll pitch yaw];
        %DCM horizon axis -> object axis
%         Mho = DCMfromAngles(X(1), X(2), 0);
%         magHeading = Mho'*mCal';
%         Z(3) = -atan2(magHeading(2),magHeading(1)) + magDeclination;
%         Z(1) = atan2(aCal(2),sqrt(aCal(1)^2+aCal(3)^2));
%         Z(2) = atan2(aCal(1),sqrt(aCal(2)^2+aCal(3)^2));

        


        
    
        % Kalman Filter equations
        Y = Z - H*X;
        K = (P*H')/(H*P*H' + R);
        P = (I - K*H)*P*(I - K*H)' + K*R*K';
        X = X + K*Y;

        % Update Gyroscopes Bias estimation
%         dw = dw + X(4:6);


    end
end

