function [ang_, anr_, bias_, quaternion] = ahrsAlg_Quat(acc, magNorm, anr, accCoefs, magCoefs, gyrCoefs)
    % �������� �������
    % �� ���� ������� ��� ������ (������ �������)
    initP = [1e-1 1e-1 1e-1 1e-3 1e-3 1e-3]; % euler
    initQ = [1e-3 1e-3 1e-3 1e-7 1e-7 1e-7]; % euler
    initR = [1e0 1e0 1e0 1e0 1e0 1e0]; % euler
    initQuat = [1 0 0 0];
    frequency = 95;
    % Filter step
    dT = 1/single(frequency);
%     accBasic = [-0.15 -0.11 9.73]; % m/s^2
%     accBasic = [0 0 1]; % m/s^2
%     magBasic = [0.367 -0.051 -0.928];
    accelerationThreshold = 0.2; % m/s^2
    fieldThreshold = 0.2; % 0.06
    magDeclination = 0.7854;
    
    % Automatic magDeclination calculation attempt
    magDeclination = getMagDeclination(acc(1,:), magNorm(1,:), accCoefs, magCoefs);
%     g = 9.81; % m/s^2

    % Initial Covariance matrix
    P = single(diag(initP));
    % System noises
    Q = single(diag(initQ));
    % Measurement noises
    R = single(diag(initR));
    % Initial Gyroscope bias
%     dw = single(0.01);

    len = length(acc);
    start = 1;
    quaternion = zeros(len, 4, 'single');
    angles = zeros(len, 3, 'single');
    debug = zeros(len, 3, 'single');
    bias_ = zeros(len, 3, 'single');
    anr_ = zeros(len, 3, 'single');
    X = zeros(len, 9, 'single');
    T = zeros(len, 1, 'single');
    accBasic = acc(2,:);
    accBasic = accBasic./norm(accBasic);
    len = length(acc);
    magBasic = magNorm(2,:);
    dw = [0 0 0];
    quaternion(1,:) = initQuat;

    for k=start:len-1
%         if k > 3400
%             sd = 0;
%         end;
        T(k+1) = T(k) + dT;

        [quaternion(k+1,:), angles(k+1, :), dw, P, Q, R, anr_(k+1,:)] = LKF_Quat(acc(k,:), magNorm(k,:), anr(k,:),...
        accCoefs, magCoefs, gyrCoefs,...
        dw, dT, quaternion(k,:), accBasic, magBasic,...
        P, Q, R);

%         bias_(k,:) = dw;
%         triad(k,:) = triadOrientation(acc(k,:), magNorm(k,:), accBasic, magBasic);
%         [angles(k+1,1), angles(k+1,2), angles(k+1,3)] = quat2angle(quaternion(k+1,:));
    end;
    
%     anr_ = 0;
    ang_ = [angles(:,3) angles(:,2) angles(:,1)];
    anr_ = X(:,4:6);
    bias_ = X(:,7:9);
    %euler
    figure
    plot(T, ang_*180/pi);
    title('angles');
%     figure
%     plot(T, anr_);
%     title('angular rate');
%     figure
%     plot(T, bias_);
%     title('gyro bias');
    figure
    plot(T, quaternion);
%     figure
%     plot(T, anr_);
% %     hold on
% %     plot(T, anr(:,1),'k');
% %     title('anr(x) comparison')
%     figure
%     plot(T, debug);
%     title('debug');
%     legend accCorrection magCorrection gyrCorrection
%     figure
%     plot(triad*180/pi);
%     title('triad');

    function magDec = getMagDeclination(accl, magn, accl_coefs, magn_coefs)
        % Calibrate magnetometer
        B = [magn_coefs(1) magn_coefs(4) magn_coefs(5)
            magn_coefs(6) magn_coefs(2) magn_coefs(7)
            magn_coefs(8) magn_coefs(9) magn_coefs(3)]; % 3x3
        B0 = [magn_coefs(10); magn_coefs(11); magn_coefs(12)]; % 3x1
        magCal = ((eye(3,'single')-B)*(magn'-B0))';

        % Calibrate accelerometer
        B = [accl_coefs(1) accl_coefs(4) accl_coefs(5)
            accl_coefs(6) accl_coefs(2) accl_coefs(7)
            accl_coefs(8) accl_coefs(9) accl_coefs(3)]; % 3x3
        B0 = [accl_coefs(10);accl_coefs(11);accl_coefs(12)]; % 3x1
        accCal = ((eye(3,'single')-B)*(accl'-B0))';
        
        roll = -atan2(accCal(2),sqrt(accCal(1)^2 + accCal(3)^2)); % ����
        pitch = atan2(accCal(1),sqrt(accCal(2)^2 + accCal(3)^2)); % ������
        Moh = angles2dcm(roll, pitch, 0); % object -> horizon
        magHeading = Moh*magCal';
        yaw = -atan2(magHeading(2),magHeading(1));
        magDec = -yaw;
    end
end