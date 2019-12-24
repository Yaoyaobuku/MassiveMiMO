function [pop] = functioncreatpop(l,tcf,num,Beta)
ibm = zeros(1,tcf);
B = Beta;
sum_Beta = sum(B);
for i=1:tcf
   ibm(i) = find(sum_Beta == max(sum_Beta));
   sum_Beta(ibm(i)) = 0;
end
pop = [];
for k =1:num
array_tcf = randperm(tcf,tcf);
pop_start = zeros(1,l);
j=1;
for i=1:l
   if(find(ibm == i) ~= 0)
       pop_start(i) = array_tcf(find(ibm == i));
       j = j+1;
   else
       pop_start(i) = randi(tcf);
   end  
end
   pop = [pop;pop_start];
end
