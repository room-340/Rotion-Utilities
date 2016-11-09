function [ value ] = corr_value_Acc_ArefAreal(ref, test)
dist=0;
for i=1:3
    norm_ref = sum(ref(:,i))/length(ref);
    norm_test = sum(test(:,i))/length(ref);
    ref(:,i)=ref(:,i)-norm_ref;
    test(:,i)=test(:,i)-norm_test;
end
for i=1:length(ref)
    dist=dist+norm_value(ref(i,1),test(i,1),0.3)+norm_value(ref(i,1),test(i,1),0.3)+norm_value(ref(i,1),test(i,1),0.3);
end
value = dist/(length(ref)*3);


end

