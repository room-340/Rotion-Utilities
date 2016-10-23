function [ q ] = quat_mult(q1,q2)
%функция перемножения кватернионов

q = [q2(1)*q1(1)-q2(2:4)*q1(2:4)' ...
        q2(1)*q1(2:4)+q1(1)*q2(2:4)+cross_multiplication(q1(2:4),q2(2:4))];

end

