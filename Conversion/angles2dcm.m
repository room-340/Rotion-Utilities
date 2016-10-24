function [DCM] = angles2dcm(anglex,angley,anglez)
%Матрица направляющих косинусов по углам, углы в последовательность XYZ
    DCM(1,1) = cos(angley)*cos(anglez);
    DCM(1,2) = -sin(anglex)*sin(angley)*cos(anglez) + cos(anglex)*sin(anglez);
    DCM(1,3) = -cos(anglex)*sin(angley)*cos(anglez) - sin(anglex)*sin(anglez);
    DCM(2,1) = -cos(angley)*sin(anglez);
    DCM(2,2) = sin(anglex)*sin(angley)*sin(anglez) + cos(anglex)*cos(anglez);
    DCM(2,3) = cos(anglex)*sin(angley)*sin(anglez) - sin(anglex)*cos(anglez);        
    DCM(3,1) = sin(angley);
    DCM(3,2) = sin(anglex)*cos(angley);
    DCM(3,3) = cos(anglex)*cos(angley);
end

