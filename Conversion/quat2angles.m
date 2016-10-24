function [angles] = quat2angles(q)
%Матрица направляющих косинусов по углам, углы в последовательность XYZ
    
    %Normalize quaternion 
    qin = q./(sqrt(sum(q.^2,2))* ones(1,4));

    %DCM elements
    c21 = 2*(qin(2)*qin(3) + qin(1)*qin(4));
    c11 = qin(1)^2 + qin(2)^2 - qin(3)^2 - qin(4)^2;
    c31 = 2*(qin(2)*qin(4) - qin(1)*qin(3));
    c32 = 2*(qin(3)*qin(4) + qin(1)*qin(2));
    c33 = qin(1)^2 - qin(2)^2 - qin(3)^2 + qin(4)^2;

    %Angles from DCM
    angles(1) = atan2(c21,c11);
    angles(2) = asin(-c31);
    angles(3) = atan2(c32,c33);
    
    [angles(1), angles(2), angles(3)]=quat2angle(q);
end


