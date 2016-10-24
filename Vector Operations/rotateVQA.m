function [vector] = rotateVQA(vector,orientation)
%Поворот вектора, при помощи кватерниона или углов в последовательности xyz

    if (length(orientation)==3)
        DCM = angles2dcm(orientation(1),orientation(2),orientation(3));
        vector = DCM*vector'; 
    end
    if (length(orientation)==4)
        rotateVQ(vector,orientation);
    end

end

