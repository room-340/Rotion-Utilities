function [X, P, Q, R, debug, quat] = LKF_QuatExt(acc, mag, anr,...
    accl_coefs, magn_coefs, gyro_coefs,...
    dT, quat, accBasic, magBasic,... 
    X, P, Q, R, accThreshold, magThreshold)
%   �������� ����������
%
%   INPUT:
%   ������� ���������:
%   acc - ������� ��������� [1x3]
%   mag - ������� ������������� ������ ���������� ���� [1x3]
%         ��� ���������� ����� |mag| = 1
%   anr - ������� ������� �������� [1x3]
%   ���������� ������� �������:
%   quat - ���������� ���������� ������� [w i j k]
%   X - ������ ��������� [1x9]=[qi qj qk wx wy wz dwx dwy wdz]
%       [qi qj qk] - ���������� ��������
%       [wx wy wz] - ������ ������� ��������
%       [dwx dwy dwz] - ������ ����� ���� ���������
%   P - ������� ���������� [9x9]
%   Q - ������� ����� ������� [9x9]
%   R - ������� ����� ��������� [9x9]
%   ���������� ��������:
%   accl_coefs - ������������� ������������ ������������� [1x12]
%   magn_coefs - ������������� ������������ ������������ [1x12]
%   gyro_coefs - ������������� ������������ ��������� [1x12]
%   dT - ��� ������������� dT = 1/�������
%   accBasic - ������������� ������ ��������� � ��������� ������ �������
%              [1x3], |accBasic| = 1
%   magBasic - ������������� ������ ���������� ���� � ��������� ������ �������
%              [1x3], |magBasic| = 1
%   accThreshold - ����� ��������� ��������� �� �������������
%   magThreshold - ����� ��������� ��������� �� ������������
%
%   OUTPUT:
%   X - ������ ��������� [1x9]
%   P - ������� ���������� [9x9]
%   Q - ������� ����� ������� [9x9]
%   R - ������� ����� ��������� [9x9]
%   debug - ��������� ������� ������ [1x3]
%           debug(1) = 1 - ��������� �� �������������
%           debug(2) = 2 - ��������� �� ������������
%           debug(3) = 3 - zero-update �� ������� ��������
%   quat - ���������� ���������� ������� [w i j k]
%% Data correction
    debug = [0 0 0];
%   X = [qe(1) qe(2) qe(3) Wx Wy Wz dWx dWy dWz]
    
    % Sensors data
    m = mag;
    a = acc;
%     w = anr;
%   ---------------------------------- CRUTCH
    w = [anr(1) -anr(2) anr(3)];
    gyro_coefs(11) = -gyro_coefs(11);

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
    
    wPrime = X(4:6);
    F = [1              wPrime(3)*dT   -wPrime(2)*dT 0  0  0  -dT  0  0;
         -wPrime(3)*dT  1              wPrime(1)*dT  0  0  0  0  -dT  0;
         wPrime(2)*dT   -wPrime(1)*dT  1             0  0  0  0  0  -dT;
         0              0              0            1  0  0 -1  0  0;
         0              0              0            0  1  0  0 -1  0;
         0              0              0            0  0  1  0  0 -1;
         0              0              0            0  0  0  1  0  0;
         0              0              0            0  0  0  0  1  0;
         0              0              0            0  0  0  0  0  1];
    
%       gyroDrift = (30/3600)*pi/180; % deg/h -> deg/sec -> rad/sec
    % system noise
    Qv = [Q(1,1) Q(2,2) Q(3,3) Q(4,4) Q(5,5) Q(6,6) Q(7,7) Q(8,8) Q(9,9)]';
    X = F*X' + Qv;
    P = F*P*F' + Q;
      
    wPrime = X(4:6); % yup, its already changed!
    % Rotate quat by qe
    quatPrime  = mrotate_eml(quat, wPrime, dT);
      
    %% Measurement update
    % Initial matrix values
    I = eye(9,'single');
    Z = zeros(9,1,'single');
    H = zeros(9,9,'single');
    
    % Estimated vector values 
    accBasicPrime = quatrotate( quatPrime, accBasic);
    magBasicPrime = quatrotate( quatPrime, magBasic);
    
    H(7:9,4:6) = eye(3, 'single');
    % Might consider deep correction mode
    % H(7:9,7:9) = eye(3, 'single');
    
    % Chek if there is non-zero data to prevent divergence
    if (norm(aCal) > single(0)) && (norm(mCal) > single(0))
        % Measurement vector
        dA = aCal./norm(aCal) - accBasicPrime;
        dM = mCal - magBasicPrime;
        Z(1:3) = dA;
        Z(4:6) = dM;
        Z(7:9) = wCal;
        % Angular rate zero-update
        for i=1:3
            if abs(wCal(i)) < 0.25
                Z(6+i) = 0;
                debug(3) = 3;
            end;
        end;
        
        % Measurement matrix with gating
        threshold_a =  ...
            abs(sqrt(aCal(1)^2 + aCal(2)^2 + aCal(3)^2) - 9.81);
        if (threshold_a < accThreshold)
            H(1,2) = -accBasicPrime(3); H(1,3) =  accBasicPrime(2);
            H(2,1) =  accBasicPrime(3); H(2,3) = -accBasicPrime(1);
            H(3,1) = -accBasicPrime(2); H(3,2) =  accBasicPrime(1);
            debug(1) = 1;
        end

        threshold_m = abs(sqrt(mCal(1)^2 + mCal(2)^2 + mCal(3)^2) - single(1.0));
        if (threshold_m < magThreshold)
            H(4,2) = -magBasicPrime(3); H(4,3) =  magBasicPrime(2);
            H(5,1) =  magBasicPrime(3); H(5,3) = -magBasicPrime(1);
            H(6,1) = -magBasicPrime(2); H(6,2) =  magBasicPrime(1);
            debug(2) = 2;
        end
        % qe is a minor quat rotation. It's calculated from scratch on
        % every step:
        X(1:3) = [0 0 0];
        
        % Kalman Filter equations
        Y = Z - H*X;
        K = (P*H')/(H*P*H' + R);
        P = (I - K*H)*P*(I - K*H)' + K*R*K';
        X = X + K*Y;
        
        % Correct quaternion
        qe  = single([X(1)/2 X(2)/2 X(3)/2]); % This is qe
        qqe = [sqrt(single(1) - (qe(1)^2 + qe(2)^2 + qe(3)^2)),...
            qe(1), qe(2), qe(3)]; % and this is full quaternion built out of it
        quat = multQQ( quatPrime, qqe );

    end
end

