function [ cross ] = multVV(vect1, vect2)
%Функция векторного умножения
cross(1) = vect1(2)*vect2(3)-vect1(3)*vect2(2);
cross(2) = vect1(3)*vect2(1)-vect1(1)*vect2(3);
cross(3) = vect1(1)*vect2(2)-vect1(2)*vect2(1);
end

