function [quaternion, anr_, bias_] = ahrsAlg_QuatExt(acc, magNorm, anr, accCoefs, magCoefs, gyrCoefs)
    % Основная функция
    % На вход податся ВСЕ данные (полные массивы)
    
    initP = [1e-1 1e-1 1e-1 1e-1 1e-1 1e-1 1e-1 1e-1 1e-1]; % quatExt
    initQ = [1e-4 1e-4 1e-4 1e-2 1e-2 1e-2 1e-6 1e-6 1e-6]; % quatExt
    initR = [1e-1 1e-1 1e-1 1e0 1e0 1e0 1e-2 1e-2 1e-2]; % quatExt
    initQuat = [1 0 0 0]; % corresponds to zero rotation
    frequency = 95;
    % Filter step
    dT = 1/single(frequency);
    accelerationThreshold = 0.2; % m/s^2
    fieldThreshold = 0.2; % 0.06

    % Initial Covariance matrix
    P = single(diag(initP));
    % System noises
    Q = single(diag(initQ));
    % Measurement noises
    R = single(diag(initR));

    len = length(acc);
    quaternion = zeros(len, 4, 'single');
    quaternion(1,:) = initQuat;
%     angles = zeros(len, 3, 'single');
    debug = zeros(len, 3, 'single');
    bias_ = zeros(len, 3, 'single');
    anr_ = zeros(len, 3, 'single');
    X = zeros(len, 9, 'single');
    T = zeros(len, 1, 'single');
    % accBasic и magBasic берут 2е значение из массива на случай если
    % первое значение массива нулевое
    accBasic = acc(2,:);
    accBasic = accBasic./norm(accBasic);
    magBasic = magNorm(2,:);   

    for k=1:len-1
%         if k > 3400
%             sd = 0; % this is used for debug only
%         end;
        T(k+1) = T(k) + dT;

        [X(k+1,:), P, Q, R, debug(k+1,:), quaternion(k+1,:)] = LKF_QuatExt(acc(k,:), magNorm(k,:), anr(k,:),...
            accCoefs, magCoefs, gyrCoefs,...
            dT, quaternion(k,:), accBasic, magBasic,... 
            X(k,:), P, Q, R, accelerationThreshold, fieldThreshold);
        
        % convert quaternions into angles
%         [angles(k+1,1), angles(k+1,2), angles(k+1,3)] = quat2angle(quaternion(k+1,:));
    end;
    
    anr_ = X(:,4:6);
    bias_ = X(:,7:9);

%     figure
%     plot(T, ang_*180/pi);
%     title('angles');
%     figure
%     plot(T, anr_);
%     title('angular rate');
%     figure
%     plot(T, bias_);
%     title('gyro bias');
    figure
    plot(T, quaternion);
    title('quaternion');

end