%S_CF = Rate_CF*20*(1-20/200)/2;
%S_CF_g= Rate_CF_genetic*20*(1-20/200)/2;
%S_SC = Rate_SC*20*(1-(20+20)/200)/2;
%avr_cfgenetic_maxrate = zeros(1,100);
%avr_cfgenetic_minpc = zeros(1,100);
avr_cfgreedy_ =zeros(1,100);
avr_random =zeros(1,100);
avr_cluster = zeros(1,100);
Rate_CF_cluster;
%Rate_CF_greedy = reshape(Rate_CF(:,7,:),30,100); 
%Rate_SC_greedy = reshape(Rate_SC(:,7,:),30,100); 
for i=1:100
 %   avr_cfgenetic_maxrate(1,i) = mean(Rate_CF_genetic_maxrate(:,i));
  %  avr_cfgenetic_minpc(1,i) = mean(Rate_CF_genetic_minpc(:,i));
    avr_cfgreedy_(1,i) = mean(Rate_CF_genetic_greedy(:,i));
    avr_random(1,i) = mean(Rate_CF_genetic_random(:,i));
    avr_cluster(1,i) = mean(Rate_CF_cluster(:,i));
end
%S = reshape(Rate_CF_genetic_maxrate,1,2000);
%S1 = reshape(Rate_CF_genetic_minpc,1,2000);
S2 = reshape(Rate_CF_genetic_greedy,1,1000);
S3 = reshape(Rate_CF_genetic_random,1,1000);
S1 = reshape(Rate_CF_cluster,1,1000);
%S = sort(S);
S1 = sort(S1);
S2 = sort(S2);
S3 = sort(S3);
%S_rate = sort(avr_cfgenetic_maxrate);
S1_rate = sort(avr_cluster);
S2_rate = sort(avr_cfgreedy_);
S3_rate = sort(avr_random);
%hold on
ecdf(S1_rate);
hold on 
ecdf(S2_rate);
hold on
ecdf(S3_rate);
title('M=60, K=20, tcf=5');
legend('cluster','greedy','random');
xlabel('Average Users Downlink Rate');
ylabel('CDF');
figure(2)
%ecdf(S);
hold on
ecdf(S1);
hold on
ecdf(S2);
hold on
ecdf(S3);
hold on
title('M=40, K=20, tcf=5');
legend('cluster','greedy','random');
xlabel('Per Users Downlink Rate');
ylabel('CDF');

