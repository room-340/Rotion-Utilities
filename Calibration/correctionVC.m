function [a_corr] = correctionVC(a, coefs)
B = eye(3)-[coefs(1) coefs(4) coefs(5); coefs(6) coefs(2) coefs(7); coefs(8) coefs(9) coefs(3)];
B0 = [coefs(10); coefs(11); coefs(12)];
a_corr = B*a'-B0;
end

