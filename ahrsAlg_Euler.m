function [ang_, anr_, bias_] = ahrsAlg_Euler(acc, magNorm, anr, accCoefs, magCoefs, gyrCoefs)
    % Основная функция
    % На вход податся ВСЕ данные (полные массивы)
    initP = [1e-1 1e-1 1e-1 1e-3 1e-3 1e-3 1e-4 1e-4 1e-4]; % euler
    initQ = [1e-6 1e-6 1e-6 1e-3 1e-3 1e-3 1e-5 1e-5 1e-5]; % euler
    initR = [2e-2 2e-2 2e-2 1e-2 1e-2 1e-2];                % euler
%     initQuat = [0 0 0 0];
    frequency = 95;
    % Filter step
    dT = 1/single(frequency);
%     accBasic = [-0.15 -0.11 9.73]; % m/s^2
%     accBasic = [0 0 1]; % m/s^2
%     magBasic = [0.367 -0.051 -0.928];
    accelerationThreshold = 0.3; % m/s^2
    fieldThreshold = 0.14; % 0.06
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
%     X = [0.694 0.021 -0.081 0 0 0 0 0 0];
%     angles = zeros(len, 3, 'single');
%     triad = zeros(len, 3, 'single');
    debug = zeros(len, 3, 'single');
    X = zeros(len, 9, 'single');
    T = zeros(len, 1, 'single');
    % accBasic = acc(start,:);
    len = length(acc);
%     magBasic = magNorm(start,:);
    mag_ = zeros(len, 1, 'single');
    
    for k=start:len-1
%         if k > 2200
%             sd = 0;
%         end;
        T(k+1) = T(k) + dT;
        [X(k+1,:), P, Q, R, debug(k,:)] = LKF_Euler(acc(k,:), magNorm(k,:), anr(k,:),...
        accCoefs, magCoefs, gyrCoefs,...
        dT,... 
        X(k,:), P, Q, R, accelerationThreshold, fieldThreshold, magDeclination);
%         triad(k,:) = triadOrientation(acc(k,:), magNorm(k,:), accBasic, magBasic);
    mag_(k) = norm(magNorm(k,:));
        
    end;
    
    ang_ = X(:,1:3);
    anr_ = X(:,4:6);
    bias_ = X(:,7:9);
    %euler
    figure
    plot(T, ang_*180/pi);
    title('angles');
    figure
    plot(T, anr_);
    title('angular rate');
    figure
    plot(T, bias_);
    title('gyro bias');
%     figure
%     plot(T, anr_(:,1),'b');
%     hold on
%     plot(T, anr(:,1),'k');
%     title('anr(x) comparison')
    figure
    plot(T, debug);
    title('debug');
    legend accCorrection magCorrection gyrCorrection
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
        
        roll = -atan2(accCal(2),sqrt(accCal(1)^2 + accCal(3)^2)); % крен
        pitch = atan2(accCal(1),sqrt(accCal(2)^2 + accCal(3)^2)); % тангаж
        Moh = angle2dcm_eml(0, pitch, roll); % object -> horizon
        magHeading = Moh*magCal';
        yaw = -atan2(magHeading(2),magHeading(1));
        magDec = -yaw;
    end
end