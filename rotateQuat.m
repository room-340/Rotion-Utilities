function [q] = rotateQuat(q, w, dt)
    Fx = w(1)*dt;
    Fy = w(2)*dt;
    Fz = w(3)*dt;
    Fm = sqrt(Fx*Fx+Fy*Fy+Fz*Fz);
    sinFm2 = sin(Fm/2);
    cosFm2 = cos(Fm/2); 
    if Fm~=single(0)
        Nd=[cosFm2, Fx/Fm*sinFm2, Fy/Fm*sinFm2, Fz/Fm*sinFm2];
    else
        Nd=single([ 1, 0, 0, 0]);
    end
    q = quat_mult( q, Nd );
end

