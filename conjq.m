function [q] = conjq(q)
%обратный кватернион
    q(2)=-q(2);
    q(3)=-q(3);
    q(4)=-q(4);
end

