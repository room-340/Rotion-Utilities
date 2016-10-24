function [v] = rotateVQ(v,q)
%Поворот вектора на кватернион
qv(1) = 0;
qv(2) = v(1);
qv(3) = v(2);
qv(4) = v(3);

qv2 = multQQ(multQQ(conjQ(q),qv),q);

v(1) = qv2(2);
v(2) = qv2(3);
v(3) = qv2(4);

end

