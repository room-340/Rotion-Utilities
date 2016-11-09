function [ value ] = corr_value_Single_RefReal(ref, test)
dist=0;
norm_ref = sum(ref)/length(ref);
norm_test = sum(test)/length(ref);
ref=ref-norm_ref;
test=test-norm_test;
for i=1:length(ref)
    dist=dist+norm_value(ref(i),test(i),0.3);
end
value = dist/(length(ref));


end

