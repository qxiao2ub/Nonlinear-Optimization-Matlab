


%-------case 1.3-------

'-------Case 1.3-------'

IniP=100;

%initially price Pijt all the same
P=zeros(20,20,40);%create 20*20*40 matrix
for i=1:i_up 
    for j=1:j_up
        for t=1:t_up
          P(i,j,t)=IniP; %0608-initial P should be infinity, I use 2 for test
        end
    end    
end

% Pijt=2;

Iter_count = 0;

% 0517-M is c16 time window, set at 1 gives it relaxation
% 0424-M is c16 time window, set at 3 gives it relaxation
M=1; 

%0515-change all var to continours instead of MIP
X = sdpvar(i_up,j_up,t_up+M); %0424-Xi,j,t+M, cuz time window constraint c16 Xijt is beyond t
%X = intvar(i_up,j_up,t_up+M); %0424-Xi,j,t+M, cuz time window constraint c16 Xijt is beyond t
X_t0 = sdpvar(i_up,j_up,1); %0424-Xi,j,0 means moves initially for c18: inventory balance eq
%X_t0 = intvar(i_up,j_up,1); %0424-Xi,j,0 means moves initially for c18: inventory balance eq
U = sdpvar(i_up,j_up,t_up);
%U = intvar(i_up,j_up,t_up); 
U_t0 = sdpvar(i_up,j_up,1); %0424-this is for Ui,j,0 
%U_t0 = intvar(i_up,j_up,1); %0424-this is for Ui,j,0 
V = sdpvar(i_up,t_up); 
%V = intvar(i_up,t_up); 
V_t0 = sdpvar(i_up,1); %0424-this is especially for Vi,0
%V_t0 = intvar(i_up,1); %0424-this is especially for Vi,0
D_ij_SAV = sdpvar(i_up,j_up,t_up); %0425-this is to show D_ij_SAV, for futher plot, then find upper bound of Pijt
%D_ij_SAV = intvar(i_up,j_up,t_up); %0425-this is to show D_ij_SAV, for futher plot, then find upper bound of Pijt
% 0604-k ratio
k = sdpvar(i_up,j_up,t_up);

% 070220-set place for Wijt
W=zeros(20,20,40);%create 20*20*40 matrix

%1.1-parameters
% a=b=m=n=1,T=40
a=1;
b=1;
m=1;
n=1;
T=8; %070120-total 40 time periods, 5 time block, so every block is 8
T_abs=8; %070120-T_abs should be same as T
% cijt=40, pijt=10, hit=10, alphaijt=0.1 (all c,p,h,alpha are subjectively to be same)
cijt=20;
pijt=10; % 070120-increase for penalty of unmet demand
% 0425-intentionally make pijt this small to avoid all 0 in Uijt solution
hit=10;
alpha=1;

%%%-----0607----
%%% How to add loop?
%%%----

% 0610-here starts the loop

P_sav = zeros(20,20,40);
U_sav = zeros(20,20,40);
min_k = 0;
P_bench = zeros(20,20,40);

% tik = 2;

satis_P = 0;
feasib_P = 0;
consec_P = 0;

% while (satis_P<i_up*j_up*t_up) || (feasib_P<i_up*j_up*t_up)
% while Iter_count < 2
while (consec_P<3)

Iter_count = Iter_count+1; 


for i=1:i_up 
    for j=1:j_up
        for t=1:t_up
          P_sav(i,j,t)=P(i,j,t); %0608-initial P should be infinity, I use 2 for test
        end
    end    
end  


for i=1:i_up 
    for j=1:j_up
        for t=1:t_up
            U_sav(i,j,t)=U(i,j,t); 
        end
    end    
end


min_k_sav = min_k;


%1.2-add operational level constraints

C = [];

% c15
for i=1:i_up %c15
   for j=1:j_up
     for t=1:t_up
       if (t==1)
         C = [C, U(i,j,1) >= U_t0(i,j,1)+(b*d(i,j,1)+cijt*d(i,j,1)*m-n*T_abs*U(i,j,1))./(a+b+m*(cijt+P(i,j,t)))-X(i,j,1)]; %c15: t=1
       else       
         C = [C, U(i,j,t) >= U(i,j,t-1)+(b*d(i,j,t)+cijt*d(i,j,t)*m-n*T_abs*U(i,j,t))./(a+b+m*(cijt+P(i,j,t)))-X(i,j,t)]; %c15: t from 2
       end
     end
   end    
end

% 0424-c16-I suppose time window is M period for c16 - relax this one, hope get solution
% this also considers Xijt beyond T
for i=1:i_up %c16
    for j=1:j_up
      for t=1:t_up
        if (t==1)
          C = [C, sum(X(i,j,t:(t+M))) >= U_t0(i,j,1)+(b*d(i,j,t)+cijt*d(i,j,t)*m-n*T_abs*U(i,j,t))./(a+b+m*(cijt+P(i,j,t)))]; %c16: t=1
        else
          C = [C, sum(X(i,j,t:(t+M))) >= U(i,j,t-1)+(b*d(i,j,t)+cijt*d(i,j,t)*m-n*T_abs*U(i,j,t))./(a+b+m*(cijt+P(i,j,t)))]; %c16: t from 2
        end
      end  
    end     
end

% 0424 - c17 is in upper level, no need in lower level
% for i=1:i_up %c17
%   for j=1:j_up
%     for t=1:t_up
%       C = [C, ((b-1)*d(i,j,t)+cijt*d(i,j,t)*m)./(n*T_abs) <= U(i,j,t)];
%       C = [C, U(i,j,t) <= (b*d(i,j,t)+cijt*d(i,j,t)*m)./(n*T_abs)];
%     end
%   end
% end

% 0424-c18-here makes <= for Vi,t cuz alpha all set to 1, too relaxed, should be variaty
for i=1:i_up %c18
  for t=1:t_up
    if (t==1)
      C = [C, V(i,t) <= V_t0(i,1) + sum(alpha*X_t0(1:j_up,i,1)) - sum(X(i,1:j_up,t)) ]; %c18: t=1
    else
      C = [C, V(i,t) <= V(i,t-1) + sum(sum(alpha*X(1:j_up,i,1:(t-1)))) - sum(X(i,1:j_up,t)) ]; %c18: t from 2 (remember to use 2 sum for 2nd)
    end
  end
end

% 0508-c19-not added, cuz it looks already relaxed in math function, hidden property of our formula

for i=1:i_up %c20: Xijt
    for j=1:j_up
        for t=1:(t_up+M)
          C = [C, X(i,j,t) >= 0];
        end
    end    
end    

for i=1:i_up %c20: X_t0(i,j,1) is Xijt=0
    for j=1:j_up
        C = [C, X_t0(i,j,1) >= 0];
    end    
end    

for i=1:i_up %c20: Uijt
    for j=1:j_up
        for t=1:t_up
            C = [C, U(i,j,t) >= 0];
        end
    end    
end    

for i=1:i_up %c20: U_t0(i,j,1) is Uijt=0
    for j=1:j_up
        C = [C, U_t0(i,j,1) >= 0];
    end    
end 

for i=1:i_up %c20: Vit
    for t=1:t_up
      C = [C, V(i,t) >= 0];
    end    
end    

for i=1:i_up %c20: V_t0(i,1) is V(i,0)
  C = [C, V_t0(i,1) >= 0];
end

%1.3-solve operational level: settings: use cplex as solver

%0515-following is set optimality for MIP, but kang said make it LNContinuous-https://www.ibm.com/support/knowledgecenter/SSSA5P_12.9.0/ilog.odms.cplex.help/CPLEX/Parameters/topics/OptimalityTarget.html
%0515-'cplex.optimalitytarget'=1: Searches for a globally optimal solution to a convex model.
ops = sdpsettings ('solver','cplex','verbose',2,'cplex.optimalitytarget',1);%0508-verbose 2 shows iterations, 1 means shown normal message, 0 shows silent message
%0515-MIP limit tolerance gap to 1-https://www.ibm.com/support/knowledgecenter/SSSA5P_12.9.0/ilog.odms.cplex.help/CPLEX/Parameters/topics/EpGap.html
%ops = sdpsettings ('solver','cplex','verbose',2,'cplex.mip.tolerances.mipgap',1);%0508-verbose 2 shows iterations, 1 means shown normal message, 0 shows silent message

%0306-limit barrier iteration-https://www.ibm.com/support/knowledgecenter/SSSA5P_12.9.0/ilog.odms.cplex.help/CPLEX/Parameters/topics/BarItLim.html
%Cplex.Param.barrier.limits.iteration = 4;[not work]

%0306-limit network simplex iteration-https://www.ibm.com/support/knowledgecenter/SSSA5P_12.9.0/ilog.odms.cplex.help/CPLEX/Parameters/topics/NetItLim.html
%Cplex.Param.network.iterations = 2;[not work]

%ops.mip.limits.nodes = 4;[not work]

%0304-use the online method for out-of-memory: https://www.ilovematlab.cn/thread-266244-1-1.html
%ops.mip.strategy.file = 3; %[not work: should put in sdpsetting] this has same results 1159110 node, but reduces time for 3 min.

%0305-use online method for out-of-memory: https://www.ibm.com/developerworks/community/forums/html/topic?id=e3ea344c-7249-4edd-a47d-29e1f80fa480
%Cplex.Param.mip.limits.auxrootthreads = 2; %[not work: should put in sdpsetting] -1 or 1

%obj-operational obj z_op
z_op = sum(sum(sum(cijt.*X(1:i_up,1:j_up,1:t_up)))) + sum(sum(cijt.*X_t0(1:i_up,1:j_up,1)))...
  + sum(sum(sum(pijt.*U(1:i_up,1:j_up,1:t_up)))) + sum(sum(pijt.*U_t0(1:i_up,1:j_up,1)))...
  + sum(sum(hit.*V(1:i_up,1:t_up))) + sum(hit.*V_t0(1:i_up,1));

%1.35-express D_ij_SAV
for i=1:i_up
  for j=1:j_up
    for t=1:t_up
      D_ij_SAV(i,j,t) = (b*d(i,j,t)+cijt*d(i,j,t)*m-n*T_abs*U(i,j,t))./(a+b+m*(cijt+P(i,j,t)));
    end
  end
end

% 0604-calculate minimum k ratio covered by SAVs
for i=1:i_up
  for j=1:j_up
    for t=1:t_up
      k(i,j,t) = D_ij_SAV(i,j,t)./d(i,j,t);
    end
  end
end

% k_min = min(k(1,1,1),k(1,2,1),k(1,3,1))

% k_min = 0.99;
% % 
% if ((k_min) >= (k(1,1,1)))
%   k_min = k(1,1,1)                
% end  

% 0604-find minimum covered ratio k
% for i=1:i_up
%   for j=1:j_up
%     for t=1:t_up
%       if (k_min >= k(i,j,t))
%         k_min = k(i,j,t)                
%       end  
%     end
%   end
% end

%1.4-check if solved and output
result = optimize(C,z_op,ops); % 0508-origin
% result = solvesdp(C,z_op,ops); % 0508-this works

%result = optimize(C,z);
% if result.problem == 0 % problem =0 means solve successfully,only print Uijt here

%   fprintf('\nk value: ') %ratio k covered by SAVs
%   value(k)

  
%   fprintf('\nP_Ini_(1,1,1) = ')
%   value(P(1,1,1))

%   fprintf('\nMin cost by Lower = ')
%   value(z_op)

  all_k = [];

  for i=1:i_up
    for j=1:j_up
      for t=1:t_up
        all_k(end+1) = k(i,j,t);
      end
    end
  end
  
  min_k = min(all_k);
  
  k_avg = mean(all_k);
  
%   fprintf('\n min k: ')
%   value(min_k)
  
  max_k = max(all_k);
  
%   fprintf('\n max k: ')
%   value(max_k)
  
  % 0615-save max_k in txt

%   fid = fopen(['IniP=',num2str(2),'max_k.txt'],'a'); %save to last 'a'
%   fprintf(fid,'%.4f\t',double(max_k)); %0516-keep 4 decimals
%   fclose(fid);
  
%   fprintf('\nU(1,1,1) value: ')
%   value(U(1,1,1))
        
%   fprintf('\nX(1,1,1) value: ') 
%   value(X(1,1,1))

fprintf('\n Iteration No. = ')
value(Iter_count)



% fprintf('\n D_ij_SAV(1,1,1) value: ') %0515-as for now hide D_ij_SAV
% value(D_ij_SAV(1,1,1))
  
%   fprintf('\n k(1,1,1) ')
%   value(k(1,1,1))

%   fprintf('\n')
%     fprintf('\nV value: ')
%     value(V)
%     
%     fprintf('\nU_t0 value: ')
%     value(U_t0)
%     
%     fprintf('\nX_t0 value: ') 
%     value(X_t0)
%     
%     fprintf('\nV_t0 value: ')
%     value(V_t0)

% else
%   disp('Not solved!')
% 
% end

% 0519-save total cost z_op in txt
% fid = fopen(['Iter=',num2str(Iter_count),'_OperCost.txt'],'wt'); %creat & write z_op to txt
% fprintf(fid,'%.4f\n',double(z_op)); %0516-keep 4 decimals
% fclose(fid);

% 0605-save min ratio k covered by SAVs in txt
% fid_ratiok = fopen(['Iter=',num2str(Iter_count),'_MinRatiok.txt'],'wt'); %creat & write z_op to txt
% fprintf(fid_ratiok,'%.4f\n',double(min_k)); %0604-keep 4 decimals
% fclose(fid_ratiok);

% 0517-following is save Uijt in txt
%1.6-save Uijt as txt for upper use
% for U_t=1:t_up 
%      fid = fopen(['IniP=',num2str(IniP),'Iter=',num2str(Iter_count),'_Uijt=',num2str(U_t),'.txt'],'wt'); %creat & write Uijt to txt
%      for i=1:i_up
%          for j=1:j_up
%              if j == j_up
%                  fprintf(fid,'%.4f\n',double(U(i,j,U_t))); %0516-keep 4 decimals
%              else
%                  fprintf(fid,'%.4f\t',double(U(i,j,U_t))); %0516-keep 4 decimals
%              end     
%          end    
%      end    
%      fclose(fid);
% end

% 0607-calculate Upper Pijt upper bound
for i=1:i_up
  for j=1:j_up
    for t=1:t_up
      P(i,j,t)=(b*d(i,j,t)+cijt*d(i,j,t)*m-n*T_abs*U(i,j,t)-(a+b)*min_k*d(i,j,t))./(m*min_k*d(i,j,t))-cijt;
    end
  end
end

% fprintf('\nUpper bound P(1,1,1) = %.4f\n',double(P(1,1,1)));

%0608-----
%calculate Upper revenue
%0608-----

revenue=0;

for i=1:i_up
  for j=1:j_up
    for t=1:t_up
      if (t==1)||(t==2)||(t==9)||(t==10)||(t==17)||(t==18)||(t==31)||(t==32)||(t==39)||(t==40)
%         T_abs=2;
        revenue = revenue-(P(i,j,t).*(b.*d(i,j,t)+cijt.*d(i,j,t)*m-n*T_abs*U(i,j,t)))./(a+b+m.*(cijt+P(i,j,t))); 
      else
%         T_abs=2; %here should be 3, but temporarily set to 2, aligh with Oper level
        revenue = revenue-(P(i,j,t).*(b.*d(i,j,t)+cijt.*d(i,j,t)*m-n*T_abs*U(i,j,t)))./(a+b+m.*(cijt+P(i,j,t)));
      end
    end
  end
end

% 070320-avg. all P

 all_P = [];

  for i=1:i_up
    for j=1:j_up
      for t=1:t_up
        all_P = P(i,j,t);
      end
    end
  end

avg_P = all_P./(i_up*j_up*t_up);



% 070320-avg. all D

 all_Dijt = [];

  for i=1:i_up
    for j=1:j_up
      for t=1:t_up
        all_Dijt = D_ij_SAV(i,j,t);
      end
    end
  end

avg_Dijt = all_Dijt./(i_up*j_up*t_up);

% 070320-print results

fprintf('\n Upper bound P(1,1,1) = %.4f\n',double(P(1,1,1))); % particular P each iteration

fprintf('\n');

fprintf('\n D(1,1,1) = %.4f\n',double(D_ij_SAV(1,1,1))); % particular Dijt each iteration

fprintf('\n');

fprintf('\n Avg_P = %.4f\n',double(avg_P)); % particular avg_P each iteration

fprintf('\n');

fprintf('\n Avg_Dijt = %.4f\n',double(avg_Dijt)); % particular avg_Dijt each iteration

fprintf('\n');

fprintf('\n Max Revenue by Upper= %.4f\n',double(-revenue));

fprintf('\n');

fprintf('\n Min Cost by Lower= %.4f\n',double(z_op));

fprintf('\n');

fprintf('\n Profit = %.4f\n',double(-revenue-z_op));

fprintf('\n');

% test how many new P > old P

satis_P = 0;

for i=1:i_up
  for j=1:j_up
    for t=1:t_up
%       if (double(P_sav(i,j,t))<=double(P(i,j,t))) % 0610-check new P greater than old ones
        %0610-if new P no more than 1 cent than old, than add to satis_P terminate condition
      %if ((double(P_sav(i,j,t)))*1.0001 >= double(P(i,j,t))) 
%       if ((double(P_sav(i,j,t)))+0.01 >= double(P(i,j,t))) 
      if ((P_sav(i,j,t))+0.01 > P(i,j,t)) 
        satis_P = (satis_P)+1;
      end  
    end
  end
end

satis_P;

fprintf('\n');

% test how many new P all feasible that new_P <= new P upper bound condition


if (Iter_count==1)
  feasib_P = 0;
else
  
  feasib_P = 0;
  
  for i=1:i_up
    for j=1:j_up
      for t=1:t_up
        P_bench(i,j,t) = (b*d(i,j,t)+cijt*d(i,j,t)*m-n*T_abs*U(i,j,t)-(a+b)*min_k*d(i,j,t))./(m*min_k*d(i,j,t))-cijt;        
      end
    end
  end
  
    
  for i=1:i_up
    for j=1:j_up
      for t=1:t_up
  %       if (double(P_sav(i,j,t))<=double(P(i,j,t))) % 0610-check new P greater than old ones
          %0610-if new P no more than 1 cent than old, than add to satis_P terminate condition
        %if ((double(P_sav(i,j,t)))*1.0001 >= double(P(i,j,t))) 
        
%         if (double(P_sav(i,j,t)) <= double(P_bench(i,j,t)))
        if (P_bench(i,j,t) >= P_sav(i,j,t))
%           i,j,t
%           fprintf("\n");
          feasib_P = (feasib_P)+1;
        end  
      end
    end
  end

end

feasib_P;

fprintf('\n');

if (satis_P==i_up*j_up*t_up) && (feasib_P==i_up*j_up*t_up)
  consec_P = 1+consec_P;
else
  consec_P = 0;
end

% tik = tik-1;

end

fprintf('\n');
fprintf('\n');
fprintf('\n');


%-------case 2.1-change b-------

'-------Case 2.1-------'

IniP=1;

%initially price Pijt all the same
P=zeros(20,20,40);%create 20*20*40 matrix
for i=1:i_up 
    for j=1:j_up
        for t=1:t_up
          P(i,j,t)=IniP; %0608-initial P should be infinity, I use 2 for test
        end
    end    
end

% Pijt=2;

Iter_count = 0;

% 0517-M is c16 time window, set at 1 gives it relaxation
% 0424-M is c16 time window, set at 3 gives it relaxation
M=1; 

%0515-change all var to continours instead of MIP
X = sdpvar(i_up,j_up,t_up+M); %0424-Xi,j,t+M, cuz time window constraint c16 Xijt is beyond t
%X = intvar(i_up,j_up,t_up+M); %0424-Xi,j,t+M, cuz time window constraint c16 Xijt is beyond t
X_t0 = sdpvar(i_up,j_up,1); %0424-Xi,j,0 means moves initially for c18: inventory balance eq
%X_t0 = intvar(i_up,j_up,1); %0424-Xi,j,0 means moves initially for c18: inventory balance eq
U = sdpvar(i_up,j_up,t_up);
%U = intvar(i_up,j_up,t_up); 
U_t0 = sdpvar(i_up,j_up,1); %0424-this is for Ui,j,0 
%U_t0 = intvar(i_up,j_up,1); %0424-this is for Ui,j,0 
V = sdpvar(i_up,t_up); 
%V = intvar(i_up,t_up); 
V_t0 = sdpvar(i_up,1); %0424-this is especially for Vi,0
%V_t0 = intvar(i_up,1); %0424-this is especially for Vi,0
D_ij_SAV = sdpvar(i_up,j_up,t_up); %0425-this is to show D_ij_SAV, for futher plot, then find upper bound of Pijt
%D_ij_SAV = intvar(i_up,j_up,t_up); %0425-this is to show D_ij_SAV, for futher plot, then find upper bound of Pijt
% 0604-k ratio
k = sdpvar(i_up,j_up,t_up);

% 070220-set place for Wijt
W=zeros(20,20,40);%create 20*20*40 matrix

%1.1-parameters
% a=b=m=n=1,T=40
a=1;
b=10;
m=1;
n=1;
T=8; %070120-total 40 time periods, 5 time block, so every block is 8
T_abs=8; %070120-T_abs should be same as T
% cijt=40, pijt=10, hit=10, alphaijt=0.1 (all c,p,h,alpha are subjectively to be same)
cijt=20;
pijt=10; % 070120-increase for penalty of unmet demand
% 0425-intentionally make pijt this small to avoid all 0 in Uijt solution
hit=10;
alpha=1;

%%%-----0607----
%%% How to add loop?
%%%----

% 0610-here starts the loop

P_sav = zeros(20,20,40);
U_sav = zeros(20,20,40);
min_k = 0;
P_bench = zeros(20,20,40);

% tik = 2;

satis_P = 0;
feasib_P = 0;
consec_P = 0;

% while (satis_P<i_up*j_up*t_up) || (feasib_P<i_up*j_up*t_up)
% while Iter_count < 2
while (consec_P<3)

Iter_count = Iter_count+1; 


for i=1:i_up 
    for j=1:j_up
        for t=1:t_up
          P_sav(i,j,t)=P(i,j,t); %0608-initial P should be infinity, I use 2 for test
        end
    end    
end  


for i=1:i_up 
    for j=1:j_up
        for t=1:t_up
            U_sav(i,j,t)=U(i,j,t); 
        end
    end    
end


min_k_sav = min_k;


%1.2-add operational level constraints

C = [];

% c15
for i=1:i_up %c15
   for j=1:j_up
     for t=1:t_up
       if (t==1)
         C = [C, U(i,j,1) >= U_t0(i,j,1)+(b*d(i,j,1)+cijt*d(i,j,1)*m-n*T_abs*U(i,j,1))./(a+b+m*(cijt+P(i,j,t)))-X(i,j,1)]; %c15: t=1
       else       
         C = [C, U(i,j,t) >= U(i,j,t-1)+(b*d(i,j,t)+cijt*d(i,j,t)*m-n*T_abs*U(i,j,t))./(a+b+m*(cijt+P(i,j,t)))-X(i,j,t)]; %c15: t from 2
       end
     end
   end    
end

% 0424-c16-I suppose time window is M period for c16 - relax this one, hope get solution
% this also considers Xijt beyond T
for i=1:i_up %c16
    for j=1:j_up
      for t=1:t_up
        if (t==1)
          C = [C, sum(X(i,j,t:(t+M))) >= U_t0(i,j,1)+(b*d(i,j,t)+cijt*d(i,j,t)*m-n*T_abs*U(i,j,t))./(a+b+m*(cijt+P(i,j,t)))]; %c16: t=1
        else
          C = [C, sum(X(i,j,t:(t+M))) >= U(i,j,t-1)+(b*d(i,j,t)+cijt*d(i,j,t)*m-n*T_abs*U(i,j,t))./(a+b+m*(cijt+P(i,j,t)))]; %c16: t from 2
        end
      end  
    end     
end

% 0424 - c17 is in upper level, no need in lower level
% for i=1:i_up %c17
%   for j=1:j_up
%     for t=1:t_up
%       C = [C, ((b-1)*d(i,j,t)+cijt*d(i,j,t)*m)./(n*T_abs) <= U(i,j,t)];
%       C = [C, U(i,j,t) <= (b*d(i,j,t)+cijt*d(i,j,t)*m)./(n*T_abs)];
%     end
%   end
% end

% 0424-c18-here makes <= for Vi,t cuz alpha all set to 1, too relaxed, should be variaty
for i=1:i_up %c18
  for t=1:t_up
    if (t==1)
      C = [C, V(i,t) <= V_t0(i,1) + sum(alpha*X_t0(1:j_up,i,1)) - sum(X(i,1:j_up,t)) ]; %c18: t=1
    else
      C = [C, V(i,t) <= V(i,t-1) + sum(sum(alpha*X(1:j_up,i,1:(t-1)))) - sum(X(i,1:j_up,t)) ]; %c18: t from 2 (remember to use 2 sum for 2nd)
    end
  end
end

% 0508-c19-not added, cuz it looks already relaxed in math function, hidden property of our formula

for i=1:i_up %c20: Xijt
    for j=1:j_up
        for t=1:(t_up+M)
          C = [C, X(i,j,t) >= 0];
        end
    end    
end    

for i=1:i_up %c20: X_t0(i,j,1) is Xijt=0
    for j=1:j_up
        C = [C, X_t0(i,j,1) >= 0];
    end    
end    

for i=1:i_up %c20: Uijt
    for j=1:j_up
        for t=1:t_up
            C = [C, U(i,j,t) >= 0];
        end
    end    
end    

for i=1:i_up %c20: U_t0(i,j,1) is Uijt=0
    for j=1:j_up
        C = [C, U_t0(i,j,1) >= 0];
    end    
end 

for i=1:i_up %c20: Vit
    for t=1:t_up
      C = [C, V(i,t) >= 0];
    end    
end    

for i=1:i_up %c20: V_t0(i,1) is V(i,0)
  C = [C, V_t0(i,1) >= 0];
end

%1.3-solve operational level: settings: use cplex as solver

%0515-following is set optimality for MIP, but kang said make it LNContinuous-https://www.ibm.com/support/knowledgecenter/SSSA5P_12.9.0/ilog.odms.cplex.help/CPLEX/Parameters/topics/OptimalityTarget.html
%0515-'cplex.optimalitytarget'=1: Searches for a globally optimal solution to a convex model.
ops = sdpsettings ('solver','cplex','verbose',2,'cplex.optimalitytarget',1);%0508-verbose 2 shows iterations, 1 means shown normal message, 0 shows silent message
%0515-MIP limit tolerance gap to 1-https://www.ibm.com/support/knowledgecenter/SSSA5P_12.9.0/ilog.odms.cplex.help/CPLEX/Parameters/topics/EpGap.html
%ops = sdpsettings ('solver','cplex','verbose',2,'cplex.mip.tolerances.mipgap',1);%0508-verbose 2 shows iterations, 1 means shown normal message, 0 shows silent message

%0306-limit barrier iteration-https://www.ibm.com/support/knowledgecenter/SSSA5P_12.9.0/ilog.odms.cplex.help/CPLEX/Parameters/topics/BarItLim.html
%Cplex.Param.barrier.limits.iteration = 4;[not work]

%0306-limit network simplex iteration-https://www.ibm.com/support/knowledgecenter/SSSA5P_12.9.0/ilog.odms.cplex.help/CPLEX/Parameters/topics/NetItLim.html
%Cplex.Param.network.iterations = 2;[not work]

%ops.mip.limits.nodes = 4;[not work]

%0304-use the online method for out-of-memory: https://www.ilovematlab.cn/thread-266244-1-1.html
%ops.mip.strategy.file = 3; %[not work: should put in sdpsetting] this has same results 1159110 node, but reduces time for 3 min.

%0305-use online method for out-of-memory: https://www.ibm.com/developerworks/community/forums/html/topic?id=e3ea344c-7249-4edd-a47d-29e1f80fa480
%Cplex.Param.mip.limits.auxrootthreads = 2; %[not work: should put in sdpsetting] -1 or 1

%obj-operational obj z_op
z_op = sum(sum(sum(cijt.*X(1:i_up,1:j_up,1:t_up)))) + sum(sum(cijt.*X_t0(1:i_up,1:j_up,1)))...
  + sum(sum(sum(pijt.*U(1:i_up,1:j_up,1:t_up)))) + sum(sum(pijt.*U_t0(1:i_up,1:j_up,1)))...
  + sum(sum(hit.*V(1:i_up,1:t_up))) + sum(hit.*V_t0(1:i_up,1));

%1.35-express D_ij_SAV
for i=1:i_up
  for j=1:j_up
    for t=1:t_up
      D_ij_SAV(i,j,t) = (b*d(i,j,t)+cijt*d(i,j,t)*m-n*T_abs*U(i,j,t))./(a+b+m*(cijt+P(i,j,t)));
    end
  end
end

% 0604-calculate minimum k ratio covered by SAVs
for i=1:i_up
  for j=1:j_up
    for t=1:t_up
      k(i,j,t) = D_ij_SAV(i,j,t)./d(i,j,t);
    end
  end
end

% k_min = min(k(1,1,1),k(1,2,1),k(1,3,1))

% k_min = 0.99;
% % 
% if ((k_min) >= (k(1,1,1)))
%   k_min = k(1,1,1)                
% end  

% 0604-find minimum covered ratio k
% for i=1:i_up
%   for j=1:j_up
%     for t=1:t_up
%       if (k_min >= k(i,j,t))
%         k_min = k(i,j,t)                
%       end  
%     end
%   end
% end

%1.4-check if solved and output
result = optimize(C,z_op,ops); % 0508-origin
% result = solvesdp(C,z_op,ops); % 0508-this works

%result = optimize(C,z);
% if result.problem == 0 % problem =0 means solve successfully,only print Uijt here

%   fprintf('\nk value: ') %ratio k covered by SAVs
%   value(k)

  
%   fprintf('\nP_Ini_(1,1,1) = ')
%   value(P(1,1,1))

%   fprintf('\nMin cost by Lower = ')
%   value(z_op)

  all_k = [];

  for i=1:i_up
    for j=1:j_up
      for t=1:t_up
        all_k(end+1) = k(i,j,t);
      end
    end
  end
  
  min_k = min(all_k);
  
  k_avg = mean(all_k);
  
%   fprintf('\n min k: ')
%   value(min_k)
  
  max_k = max(all_k);
  
%   fprintf('\n max k: ')
%   value(max_k)
  
  % 0615-save max_k in txt

%   fid = fopen(['IniP=',num2str(2),'max_k.txt'],'a'); %save to last 'a'
%   fprintf(fid,'%.4f\t',double(max_k)); %0516-keep 4 decimals
%   fclose(fid);
  
%   fprintf('\nU(1,1,1) value: ')
%   value(U(1,1,1))
        
%   fprintf('\nX(1,1,1) value: ') 
%   value(X(1,1,1))

fprintf('\n Iteration No. = ')
value(Iter_count)



% fprintf('\n D_ij_SAV(1,1,1) value: ') %0515-as for now hide D_ij_SAV
% value(D_ij_SAV(1,1,1))
  
%   fprintf('\n k(1,1,1) ')
%   value(k(1,1,1))

%   fprintf('\n')
%     fprintf('\nV value: ')
%     value(V)
%     
%     fprintf('\nU_t0 value: ')
%     value(U_t0)
%     
%     fprintf('\nX_t0 value: ') 
%     value(X_t0)
%     
%     fprintf('\nV_t0 value: ')
%     value(V_t0)

% else
%   disp('Not solved!')
% 
% end

% 0519-save total cost z_op in txt
% fid = fopen(['Iter=',num2str(Iter_count),'_OperCost.txt'],'wt'); %creat & write z_op to txt
% fprintf(fid,'%.4f\n',double(z_op)); %0516-keep 4 decimals
% fclose(fid);

% 0605-save min ratio k covered by SAVs in txt
% fid_ratiok = fopen(['Iter=',num2str(Iter_count),'_MinRatiok.txt'],'wt'); %creat & write z_op to txt
% fprintf(fid_ratiok,'%.4f\n',double(min_k)); %0604-keep 4 decimals
% fclose(fid_ratiok);

% 0517-following is save Uijt in txt
%1.6-save Uijt as txt for upper use
% for U_t=1:t_up 
%      fid = fopen(['IniP=',num2str(IniP),'Iter=',num2str(Iter_count),'_Uijt=',num2str(U_t),'.txt'],'wt'); %creat & write Uijt to txt
%      for i=1:i_up
%          for j=1:j_up
%              if j == j_up
%                  fprintf(fid,'%.4f\n',double(U(i,j,U_t))); %0516-keep 4 decimals
%              else
%                  fprintf(fid,'%.4f\t',double(U(i,j,U_t))); %0516-keep 4 decimals
%              end     
%          end    
%      end    
%      fclose(fid);
% end

% 0607-calculate Upper Pijt upper bound
for i=1:i_up
  for j=1:j_up
    for t=1:t_up
      P(i,j,t)=(b*d(i,j,t)+cijt*d(i,j,t)*m-n*T_abs*U(i,j,t)-(a+b)*min_k*d(i,j,t))./(m*min_k*d(i,j,t))-cijt;
    end
  end
end

% fprintf('\nUpper bound P(1,1,1) = %.4f\n',double(P(1,1,1)));

%0608-----
%calculate Upper revenue
%0608-----

revenue=0;

for i=1:i_up
  for j=1:j_up
    for t=1:t_up
      if (t==1)||(t==2)||(t==9)||(t==10)||(t==17)||(t==18)||(t==31)||(t==32)||(t==39)||(t==40)
%         T_abs=2;
        revenue = revenue-(P(i,j,t).*(b.*d(i,j,t)+cijt.*d(i,j,t)*m-n*T_abs*U(i,j,t)))./(a+b+m.*(cijt+P(i,j,t))); 
      else
%         T_abs=2; %here should be 3, but temporarily set to 2, aligh with Oper level
        revenue = revenue-(P(i,j,t).*(b.*d(i,j,t)+cijt.*d(i,j,t)*m-n*T_abs*U(i,j,t)))./(a+b+m.*(cijt+P(i,j,t)));
      end
    end
  end
end

% 070320-avg. all P

 all_P = [];

  for i=1:i_up
    for j=1:j_up
      for t=1:t_up
        all_P = P(i,j,t);
      end
    end
  end

avg_P = all_P./(i_up*j_up*t_up);



% 070320-avg. all D

 all_Dijt = [];

  for i=1:i_up
    for j=1:j_up
      for t=1:t_up
        all_Dijt = D_ij_SAV(i,j,t);
      end
    end
  end

avg_Dijt = all_Dijt./(i_up*j_up*t_up);

% 070320-print results

fprintf('\n Upper bound P(1,1,1) = %.4f\n',double(P(1,1,1))); % particular P each iteration

fprintf('\n');

fprintf('\n D(1,1,1) = %.4f\n',double(D_ij_SAV(1,1,1))); % particular Dijt each iteration

fprintf('\n');

fprintf('\n Avg_P = %.4f\n',double(avg_P)); % particular avg_P each iteration

fprintf('\n');

fprintf('\n Avg_Dijt = %.4f\n',double(avg_Dijt)); % particular avg_Dijt each iteration

fprintf('\n');

fprintf('\n Max Revenue by Upper= %.4f\n',double(-revenue));

fprintf('\n');

fprintf('\n Min Cost by Lower= %.4f\n',double(z_op));

fprintf('\n');

fprintf('\n Profit = %.4f\n',double(-revenue-z_op));

fprintf('\n');

% test how many new P > old P

satis_P = 0;

for i=1:i_up
  for j=1:j_up
    for t=1:t_up
%       if (double(P_sav(i,j,t))<=double(P(i,j,t))) % 0610-check new P greater than old ones
        %0610-if new P no more than 1 cent than old, than add to satis_P terminate condition
      %if ((double(P_sav(i,j,t)))*1.0001 >= double(P(i,j,t))) 
%       if ((double(P_sav(i,j,t)))+0.01 >= double(P(i,j,t))) 
      if ((P_sav(i,j,t))+0.01 > P(i,j,t)) 
        satis_P = (satis_P)+1;
      end  
    end
  end
end

satis_P;

fprintf('\n');

% test how many new P all feasible that new_P <= new P upper bound condition


if (Iter_count==1)
  feasib_P = 0;
else
  
  feasib_P = 0;
  
  for i=1:i_up
    for j=1:j_up
      for t=1:t_up
        P_bench(i,j,t) = (b*d(i,j,t)+cijt*d(i,j,t)*m-n*T_abs*U(i,j,t)-(a+b)*min_k*d(i,j,t))./(m*min_k*d(i,j,t))-cijt;        
      end
    end
  end
  
    
  for i=1:i_up
    for j=1:j_up
      for t=1:t_up
  %       if (double(P_sav(i,j,t))<=double(P(i,j,t))) % 0610-check new P greater than old ones
          %0610-if new P no more than 1 cent than old, than add to satis_P terminate condition
        %if ((double(P_sav(i,j,t)))*1.0001 >= double(P(i,j,t))) 
        
%         if (double(P_sav(i,j,t)) <= double(P_bench(i,j,t)))
        if (P_bench(i,j,t) >= P_sav(i,j,t))
%           i,j,t
%           fprintf("\n");
          feasib_P = (feasib_P)+1;
        end  
      end
    end
  end

end

feasib_P;

fprintf('\n');

if (satis_P==i_up*j_up*t_up) && (feasib_P==i_up*j_up*t_up)
  consec_P = 1+consec_P;
else
  consec_P = 0;
end

% tik = tik-1;

end

fprintf('\n');
fprintf('\n');
fprintf('\n');



%-------case 2.2-b=10------

'-------Case 2.2------'

IniP=10;

%initially price Pijt all the same
P=zeros(20,20,40);%create 20*20*40 matrix
for i=1:i_up 
    for j=1:j_up
        for t=1:t_up
          P(i,j,t)=IniP; %0608-initial P should be infinity, I use 2 for test
        end
    end    
end

% Pijt=2;

Iter_count = 0;

% 0517-M is c16 time window, set at 1 gives it relaxation
% 0424-M is c16 time window, set at 3 gives it relaxation
M=1; 

%0515-change all var to continours instead of MIP
X = sdpvar(i_up,j_up,t_up+M); %0424-Xi,j,t+M, cuz time window constraint c16 Xijt is beyond t
%X = intvar(i_up,j_up,t_up+M); %0424-Xi,j,t+M, cuz time window constraint c16 Xijt is beyond t
X_t0 = sdpvar(i_up,j_up,1); %0424-Xi,j,0 means moves initially for c18: inventory balance eq
%X_t0 = intvar(i_up,j_up,1); %0424-Xi,j,0 means moves initially for c18: inventory balance eq
U = sdpvar(i_up,j_up,t_up);
%U = intvar(i_up,j_up,t_up); 
U_t0 = sdpvar(i_up,j_up,1); %0424-this is for Ui,j,0 
%U_t0 = intvar(i_up,j_up,1); %0424-this is for Ui,j,0 
V = sdpvar(i_up,t_up); 
%V = intvar(i_up,t_up); 
V_t0 = sdpvar(i_up,1); %0424-this is especially for Vi,0
%V_t0 = intvar(i_up,1); %0424-this is especially for Vi,0
D_ij_SAV = sdpvar(i_up,j_up,t_up); %0425-this is to show D_ij_SAV, for futher plot, then find upper bound of Pijt
%D_ij_SAV = intvar(i_up,j_up,t_up); %0425-this is to show D_ij_SAV, for futher plot, then find upper bound of Pijt
% 0604-k ratio
k = sdpvar(i_up,j_up,t_up);

% 070220-set place for Wijt
W=zeros(20,20,40);%create 20*20*40 matrix

%1.1-parameters
% a=b=m=n=1,T=40
a=1;
b=10;
m=1;
n=1;
T=8; %070120-total 40 time periods, 5 time block, so every block is 8
T_abs=8; %070120-T_abs should be same as T
% cijt=40, pijt=10, hit=10, alphaijt=0.1 (all c,p,h,alpha are subjectively to be same)
cijt=20;
pijt=10; % 070120-increase for penalty of unmet demand
% 0425-intentionally make pijt this small to avoid all 0 in Uijt solution
hit=10;
alpha=1;

%%%-----0607----
%%% How to add loop?
%%%----

% 0610-here starts the loop

P_sav = zeros(20,20,40);
U_sav = zeros(20,20,40);
min_k = 0;
P_bench = zeros(20,20,40);

% tik = 2;

satis_P = 0;
feasib_P = 0;
consec_P = 0;

% while (satis_P<i_up*j_up*t_up) || (feasib_P<i_up*j_up*t_up)
% while Iter_count < 2
while (consec_P<3)

Iter_count = Iter_count+1; 


for i=1:i_up 
    for j=1:j_up
        for t=1:t_up
          P_sav(i,j,t)=P(i,j,t); %0608-initial P should be infinity, I use 2 for test
        end
    end    
end  


for i=1:i_up 
    for j=1:j_up
        for t=1:t_up
            U_sav(i,j,t)=U(i,j,t); 
        end
    end    
end


min_k_sav = min_k;


%1.2-add operational level constraints

C = [];

% c15
for i=1:i_up %c15
   for j=1:j_up
     for t=1:t_up
       if (t==1)
         C = [C, U(i,j,1) >= U_t0(i,j,1)+(b*d(i,j,1)+cijt*d(i,j,1)*m-n*T_abs*U(i,j,1))./(a+b+m*(cijt+P(i,j,t)))-X(i,j,1)]; %c15: t=1
       else       
         C = [C, U(i,j,t) >= U(i,j,t-1)+(b*d(i,j,t)+cijt*d(i,j,t)*m-n*T_abs*U(i,j,t))./(a+b+m*(cijt+P(i,j,t)))-X(i,j,t)]; %c15: t from 2
       end
     end
   end    
end

% 0424-c16-I suppose time window is M period for c16 - relax this one, hope get solution
% this also considers Xijt beyond T
for i=1:i_up %c16
    for j=1:j_up
      for t=1:t_up
        if (t==1)
          C = [C, sum(X(i,j,t:(t+M))) >= U_t0(i,j,1)+(b*d(i,j,t)+cijt*d(i,j,t)*m-n*T_abs*U(i,j,t))./(a+b+m*(cijt+P(i,j,t)))]; %c16: t=1
        else
          C = [C, sum(X(i,j,t:(t+M))) >= U(i,j,t-1)+(b*d(i,j,t)+cijt*d(i,j,t)*m-n*T_abs*U(i,j,t))./(a+b+m*(cijt+P(i,j,t)))]; %c16: t from 2
        end
      end  
    end     
end

% 0424 - c17 is in upper level, no need in lower level
% for i=1:i_up %c17
%   for j=1:j_up
%     for t=1:t_up
%       C = [C, ((b-1)*d(i,j,t)+cijt*d(i,j,t)*m)./(n*T_abs) <= U(i,j,t)];
%       C = [C, U(i,j,t) <= (b*d(i,j,t)+cijt*d(i,j,t)*m)./(n*T_abs)];
%     end
%   end
% end

% 0424-c18-here makes <= for Vi,t cuz alpha all set to 1, too relaxed, should be variaty
for i=1:i_up %c18
  for t=1:t_up
    if (t==1)
      C = [C, V(i,t) <= V_t0(i,1) + sum(alpha*X_t0(1:j_up,i,1)) - sum(X(i,1:j_up,t)) ]; %c18: t=1
    else
      C = [C, V(i,t) <= V(i,t-1) + sum(sum(alpha*X(1:j_up,i,1:(t-1)))) - sum(X(i,1:j_up,t)) ]; %c18: t from 2 (remember to use 2 sum for 2nd)
    end
  end
end

% 0508-c19-not added, cuz it looks already relaxed in math function, hidden property of our formula

for i=1:i_up %c20: Xijt
    for j=1:j_up
        for t=1:(t_up+M)
          C = [C, X(i,j,t) >= 0];
        end
    end    
end    

for i=1:i_up %c20: X_t0(i,j,1) is Xijt=0
    for j=1:j_up
        C = [C, X_t0(i,j,1) >= 0];
    end    
end    

for i=1:i_up %c20: Uijt
    for j=1:j_up
        for t=1:t_up
            C = [C, U(i,j,t) >= 0];
        end
    end    
end    

for i=1:i_up %c20: U_t0(i,j,1) is Uijt=0
    for j=1:j_up
        C = [C, U_t0(i,j,1) >= 0];
    end    
end 

for i=1:i_up %c20: Vit
    for t=1:t_up
      C = [C, V(i,t) >= 0];
    end    
end    

for i=1:i_up %c20: V_t0(i,1) is V(i,0)
  C = [C, V_t0(i,1) >= 0];
end

%1.3-solve operational level: settings: use cplex as solver

%0515-following is set optimality for MIP, but kang said make it LNContinuous-https://www.ibm.com/support/knowledgecenter/SSSA5P_12.9.0/ilog.odms.cplex.help/CPLEX/Parameters/topics/OptimalityTarget.html
%0515-'cplex.optimalitytarget'=1: Searches for a globally optimal solution to a convex model.
ops = sdpsettings ('solver','cplex','verbose',2,'cplex.optimalitytarget',1);%0508-verbose 2 shows iterations, 1 means shown normal message, 0 shows silent message
%0515-MIP limit tolerance gap to 1-https://www.ibm.com/support/knowledgecenter/SSSA5P_12.9.0/ilog.odms.cplex.help/CPLEX/Parameters/topics/EpGap.html
%ops = sdpsettings ('solver','cplex','verbose',2,'cplex.mip.tolerances.mipgap',1);%0508-verbose 2 shows iterations, 1 means shown normal message, 0 shows silent message

%0306-limit barrier iteration-https://www.ibm.com/support/knowledgecenter/SSSA5P_12.9.0/ilog.odms.cplex.help/CPLEX/Parameters/topics/BarItLim.html
%Cplex.Param.barrier.limits.iteration = 4;[not work]

%0306-limit network simplex iteration-https://www.ibm.com/support/knowledgecenter/SSSA5P_12.9.0/ilog.odms.cplex.help/CPLEX/Parameters/topics/NetItLim.html
%Cplex.Param.network.iterations = 2;[not work]

%ops.mip.limits.nodes = 4;[not work]

%0304-use the online method for out-of-memory: https://www.ilovematlab.cn/thread-266244-1-1.html
%ops.mip.strategy.file = 3; %[not work: should put in sdpsetting] this has same results 1159110 node, but reduces time for 3 min.

%0305-use online method for out-of-memory: https://www.ibm.com/developerworks/community/forums/html/topic?id=e3ea344c-7249-4edd-a47d-29e1f80fa480
%Cplex.Param.mip.limits.auxrootthreads = 2; %[not work: should put in sdpsetting] -1 or 1

%obj-operational obj z_op
z_op = sum(sum(sum(cijt.*X(1:i_up,1:j_up,1:t_up)))) + sum(sum(cijt.*X_t0(1:i_up,1:j_up,1)))...
  + sum(sum(sum(pijt.*U(1:i_up,1:j_up,1:t_up)))) + sum(sum(pijt.*U_t0(1:i_up,1:j_up,1)))...
  + sum(sum(hit.*V(1:i_up,1:t_up))) + sum(hit.*V_t0(1:i_up,1));

%1.35-express D_ij_SAV
for i=1:i_up
  for j=1:j_up
    for t=1:t_up
      D_ij_SAV(i,j,t) = (b*d(i,j,t)+cijt*d(i,j,t)*m-n*T_abs*U(i,j,t))./(a+b+m*(cijt+P(i,j,t)));
    end
  end
end

% 0604-calculate minimum k ratio covered by SAVs
for i=1:i_up
  for j=1:j_up
    for t=1:t_up
      k(i,j,t) = D_ij_SAV(i,j,t)./d(i,j,t);
    end
  end
end

% k_min = min(k(1,1,1),k(1,2,1),k(1,3,1))

% k_min = 0.99;
% % 
% if ((k_min) >= (k(1,1,1)))
%   k_min = k(1,1,1)                
% end  

% 0604-find minimum covered ratio k
% for i=1:i_up
%   for j=1:j_up
%     for t=1:t_up
%       if (k_min >= k(i,j,t))
%         k_min = k(i,j,t)                
%       end  
%     end
%   end
% end

%1.4-check if solved and output
result = optimize(C,z_op,ops); % 0508-origin
% result = solvesdp(C,z_op,ops); % 0508-this works

%result = optimize(C,z);
% if result.problem == 0 % problem =0 means solve successfully,only print Uijt here

%   fprintf('\nk value: ') %ratio k covered by SAVs
%   value(k)

  
%   fprintf('\nP_Ini_(1,1,1) = ')
%   value(P(1,1,1))

%   fprintf('\nMin cost by Lower = ')
%   value(z_op)

  all_k = [];

  for i=1:i_up
    for j=1:j_up
      for t=1:t_up
        all_k(end+1) = k(i,j,t);
      end
    end
  end
  
  min_k = min(all_k);
  
  k_avg = mean(all_k);
  
%   fprintf('\n min k: ')
%   value(min_k)
  
  max_k = max(all_k);
  
%   fprintf('\n max k: ')
%   value(max_k)
  
  % 0615-save max_k in txt

%   fid = fopen(['IniP=',num2str(2),'max_k.txt'],'a'); %save to last 'a'
%   fprintf(fid,'%.4f\t',double(max_k)); %0516-keep 4 decimals
%   fclose(fid);
  
%   fprintf('\nU(1,1,1) value: ')
%   value(U(1,1,1))
        
%   fprintf('\nX(1,1,1) value: ') 
%   value(X(1,1,1))

fprintf('\n Iteration No. = ')
value(Iter_count)



% fprintf('\n D_ij_SAV(1,1,1) value: ') %0515-as for now hide D_ij_SAV
% value(D_ij_SAV(1,1,1))
  
%   fprintf('\n k(1,1,1) ')
%   value(k(1,1,1))

%   fprintf('\n')
%     fprintf('\nV value: ')
%     value(V)
%     
%     fprintf('\nU_t0 value: ')
%     value(U_t0)
%     
%     fprintf('\nX_t0 value: ') 
%     value(X_t0)
%     
%     fprintf('\nV_t0 value: ')
%     value(V_t0)

% else
%   disp('Not solved!')
% 
% end

% 0519-save total cost z_op in txt
% fid = fopen(['Iter=',num2str(Iter_count),'_OperCost.txt'],'wt'); %creat & write z_op to txt
% fprintf(fid,'%.4f\n',double(z_op)); %0516-keep 4 decimals
% fclose(fid);

% 0605-save min ratio k covered by SAVs in txt
% fid_ratiok = fopen(['Iter=',num2str(Iter_count),'_MinRatiok.txt'],'wt'); %creat & write z_op to txt
% fprintf(fid_ratiok,'%.4f\n',double(min_k)); %0604-keep 4 decimals
% fclose(fid_ratiok);

% 0517-following is save Uijt in txt
%1.6-save Uijt as txt for upper use
% for U_t=1:t_up 
%      fid = fopen(['IniP=',num2str(IniP),'Iter=',num2str(Iter_count),'_Uijt=',num2str(U_t),'.txt'],'wt'); %creat & write Uijt to txt
%      for i=1:i_up
%          for j=1:j_up
%              if j == j_up
%                  fprintf(fid,'%.4f\n',double(U(i,j,U_t))); %0516-keep 4 decimals
%              else
%                  fprintf(fid,'%.4f\t',double(U(i,j,U_t))); %0516-keep 4 decimals
%              end     
%          end    
%      end    
%      fclose(fid);
% end

% 0607-calculate Upper Pijt upper bound
for i=1:i_up
  for j=1:j_up
    for t=1:t_up
      P(i,j,t)=(b*d(i,j,t)+cijt*d(i,j,t)*m-n*T_abs*U(i,j,t)-(a+b)*min_k*d(i,j,t))./(m*min_k*d(i,j,t))-cijt;
    end
  end
end

% fprintf('\nUpper bound P(1,1,1) = %.4f\n',double(P(1,1,1)));

%0608-----
%calculate Upper revenue
%0608-----

revenue=0;

for i=1:i_up
  for j=1:j_up
    for t=1:t_up
      if (t==1)||(t==2)||(t==9)||(t==10)||(t==17)||(t==18)||(t==31)||(t==32)||(t==39)||(t==40)
%         T_abs=2;
        revenue = revenue-(P(i,j,t).*(b.*d(i,j,t)+cijt.*d(i,j,t)*m-n*T_abs*U(i,j,t)))./(a+b+m.*(cijt+P(i,j,t))); 
      else
%         T_abs=2; %here should be 3, but temporarily set to 2, aligh with Oper level
        revenue = revenue-(P(i,j,t).*(b.*d(i,j,t)+cijt.*d(i,j,t)*m-n*T_abs*U(i,j,t)))./(a+b+m.*(cijt+P(i,j,t)));
      end
    end
  end
end

% 070320-avg. all P

 all_P = [];

  for i=1:i_up
    for j=1:j_up
      for t=1:t_up
        all_P = P(i,j,t);
      end
    end
  end

avg_P = all_P./(i_up*j_up*t_up);



% 070320-avg. all D

 all_Dijt = [];

  for i=1:i_up
    for j=1:j_up
      for t=1:t_up
        all_Dijt = D_ij_SAV(i,j,t);
      end
    end
  end

avg_Dijt = all_Dijt./(i_up*j_up*t_up);

% 070320-print results

fprintf('\n Upper bound P(1,1,1) = %.4f\n',double(P(1,1,1))); % particular P each iteration

fprintf('\n');

fprintf('\n D(1,1,1) = %.4f\n',double(D_ij_SAV(1,1,1))); % particular Dijt each iteration

fprintf('\n');

fprintf('\n Avg_P = %.4f\n',double(avg_P)); % particular avg_P each iteration

fprintf('\n');

fprintf('\n Avg_Dijt = %.4f\n',double(avg_Dijt)); % particular avg_Dijt each iteration

fprintf('\n');

fprintf('\n Max Revenue by Upper= %.4f\n',double(-revenue));

fprintf('\n');

fprintf('\n Min Cost by Lower= %.4f\n',double(z_op));

fprintf('\n');

fprintf('\n Profit = %.4f\n',double(-revenue-z_op));

fprintf('\n');

% test how many new P > old P

satis_P = 0;

for i=1:i_up
  for j=1:j_up
    for t=1:t_up
%       if (double(P_sav(i,j,t))<=double(P(i,j,t))) % 0610-check new P greater than old ones
        %0610-if new P no more than 1 cent than old, than add to satis_P terminate condition
      %if ((double(P_sav(i,j,t)))*1.0001 >= double(P(i,j,t))) 
%       if ((double(P_sav(i,j,t)))+0.01 >= double(P(i,j,t))) 
      if ((P_sav(i,j,t))+0.01 > P(i,j,t)) 
        satis_P = (satis_P)+1;
      end  
    end
  end
end

satis_P;

fprintf('\n');

% test how many new P all feasible that new_P <= new P upper bound condition


if (Iter_count==1)
  feasib_P = 0;
else
  
  feasib_P = 0;
  
  for i=1:i_up
    for j=1:j_up
      for t=1:t_up
        P_bench(i,j,t) = (b*d(i,j,t)+cijt*d(i,j,t)*m-n*T_abs*U(i,j,t)-(a+b)*min_k*d(i,j,t))./(m*min_k*d(i,j,t))-cijt;        
      end
    end
  end
  
    
  for i=1:i_up
    for j=1:j_up
      for t=1:t_up
  %       if (double(P_sav(i,j,t))<=double(P(i,j,t))) % 0610-check new P greater than old ones
          %0610-if new P no more than 1 cent than old, than add to satis_P terminate condition
        %if ((double(P_sav(i,j,t)))*1.0001 >= double(P(i,j,t))) 
        
%         if (double(P_sav(i,j,t)) <= double(P_bench(i,j,t)))
        if (P_bench(i,j,t) >= P_sav(i,j,t))
%           i,j,t
%           fprintf("\n");
          feasib_P = (feasib_P)+1;
        end  
      end
    end
  end

end

feasib_P;

fprintf('\n');

if (satis_P==i_up*j_up*t_up) && (feasib_P==i_up*j_up*t_up)
  consec_P = 1+consec_P;
else
  consec_P = 0;
end

% tik = tik-1;

end

fprintf('\n');
fprintf('\n');
fprintf('\n');


%-------case 2.3-b=10------

'-------Case 2.3-------'

IniP=100;

%initially price Pijt all the same
P=zeros(20,20,40);%create 20*20*40 matrix
for i=1:i_up 
    for j=1:j_up
        for t=1:t_up
          P(i,j,t)=IniP; %0608-initial P should be infinity, I use 2 for test
        end
    end    
end

% Pijt=2;

Iter_count = 0;

% 0517-M is c16 time window, set at 1 gives it relaxation
% 0424-M is c16 time window, set at 3 gives it relaxation
M=1; 

%0515-change all var to continours instead of MIP
X = sdpvar(i_up,j_up,t_up+M); %0424-Xi,j,t+M, cuz time window constraint c16 Xijt is beyond t
%X = intvar(i_up,j_up,t_up+M); %0424-Xi,j,t+M, cuz time window constraint c16 Xijt is beyond t
X_t0 = sdpvar(i_up,j_up,1); %0424-Xi,j,0 means moves initially for c18: inventory balance eq
%X_t0 = intvar(i_up,j_up,1); %0424-Xi,j,0 means moves initially for c18: inventory balance eq
U = sdpvar(i_up,j_up,t_up);
%U = intvar(i_up,j_up,t_up); 
U_t0 = sdpvar(i_up,j_up,1); %0424-this is for Ui,j,0 
%U_t0 = intvar(i_up,j_up,1); %0424-this is for Ui,j,0 
V = sdpvar(i_up,t_up); 
%V = intvar(i_up,t_up); 
V_t0 = sdpvar(i_up,1); %0424-this is especially for Vi,0
%V_t0 = intvar(i_up,1); %0424-this is especially for Vi,0
D_ij_SAV = sdpvar(i_up,j_up,t_up); %0425-this is to show D_ij_SAV, for futher plot, then find upper bound of Pijt
%D_ij_SAV = intvar(i_up,j_up,t_up); %0425-this is to show D_ij_SAV, for futher plot, then find upper bound of Pijt
% 0604-k ratio
k = sdpvar(i_up,j_up,t_up);

% 070220-set place for Wijt
W=zeros(20,20,40);%create 20*20*40 matrix

%1.1-parameters
% a=b=m=n=1,T=40
a=1;
b=10;
m=1;
n=1;
T=8; %070120-total 40 time periods, 5 time block, so every block is 8
T_abs=8; %070120-T_abs should be same as T
% cijt=40, pijt=10, hit=10, alphaijt=0.1 (all c,p,h,alpha are subjectively to be same)
cijt=20;
pijt=10; % 070120-increase for penalty of unmet demand
% 0425-intentionally make pijt this small to avoid all 0 in Uijt solution
hit=10;
alpha=1;

%%%-----0607----
%%% How to add loop?
%%%----

% 0610-here starts the loop

P_sav = zeros(20,20,40);
U_sav = zeros(20,20,40);
min_k = 0;
P_bench = zeros(20,20,40);

% tik = 2;

satis_P = 0;
feasib_P = 0;
consec_P = 0;

% while (satis_P<i_up*j_up*t_up) || (feasib_P<i_up*j_up*t_up)
% while Iter_count < 2
while (consec_P<3)

Iter_count = Iter_count+1; 


for i=1:i_up 
    for j=1:j_up
        for t=1:t_up
          P_sav(i,j,t)=P(i,j,t); %0608-initial P should be infinity, I use 2 for test
        end
    end    
end  


for i=1:i_up 
    for j=1:j_up
        for t=1:t_up
            U_sav(i,j,t)=U(i,j,t); 
        end
    end    
end


min_k_sav = min_k;


%1.2-add operational level constraints

C = [];

% c15
for i=1:i_up %c15
   for j=1:j_up
     for t=1:t_up
       if (t==1)
         C = [C, U(i,j,1) >= U_t0(i,j,1)+(b*d(i,j,1)+cijt*d(i,j,1)*m-n*T_abs*U(i,j,1))./(a+b+m*(cijt+P(i,j,t)))-X(i,j,1)]; %c15: t=1
       else       
         C = [C, U(i,j,t) >= U(i,j,t-1)+(b*d(i,j,t)+cijt*d(i,j,t)*m-n*T_abs*U(i,j,t))./(a+b+m*(cijt+P(i,j,t)))-X(i,j,t)]; %c15: t from 2
       end
     end
   end    
end

% 0424-c16-I suppose time window is M period for c16 - relax this one, hope get solution
% this also considers Xijt beyond T
for i=1:i_up %c16
    for j=1:j_up
      for t=1:t_up
        if (t==1)
          C = [C, sum(X(i,j,t:(t+M))) >= U_t0(i,j,1)+(b*d(i,j,t)+cijt*d(i,j,t)*m-n*T_abs*U(i,j,t))./(a+b+m*(cijt+P(i,j,t)))]; %c16: t=1
        else
          C = [C, sum(X(i,j,t:(t+M))) >= U(i,j,t-1)+(b*d(i,j,t)+cijt*d(i,j,t)*m-n*T_abs*U(i,j,t))./(a+b+m*(cijt+P(i,j,t)))]; %c16: t from 2
        end
      end  
    end     
end

% 0424 - c17 is in upper level, no need in lower level
% for i=1:i_up %c17
%   for j=1:j_up
%     for t=1:t_up
%       C = [C, ((b-1)*d(i,j,t)+cijt*d(i,j,t)*m)./(n*T_abs) <= U(i,j,t)];
%       C = [C, U(i,j,t) <= (b*d(i,j,t)+cijt*d(i,j,t)*m)./(n*T_abs)];
%     end
%   end
% end

% 0424-c18-here makes <= for Vi,t cuz alpha all set to 1, too relaxed, should be variaty
for i=1:i_up %c18
  for t=1:t_up
    if (t==1)
      C = [C, V(i,t) <= V_t0(i,1) + sum(alpha*X_t0(1:j_up,i,1)) - sum(X(i,1:j_up,t)) ]; %c18: t=1
    else
      C = [C, V(i,t) <= V(i,t-1) + sum(sum(alpha*X(1:j_up,i,1:(t-1)))) - sum(X(i,1:j_up,t)) ]; %c18: t from 2 (remember to use 2 sum for 2nd)
    end
  end
end

% 0508-c19-not added, cuz it looks already relaxed in math function, hidden property of our formula

for i=1:i_up %c20: Xijt
    for j=1:j_up
        for t=1:(t_up+M)
          C = [C, X(i,j,t) >= 0];
        end
    end    
end    

for i=1:i_up %c20: X_t0(i,j,1) is Xijt=0
    for j=1:j_up
        C = [C, X_t0(i,j,1) >= 0];
    end    
end    

for i=1:i_up %c20: Uijt
    for j=1:j_up
        for t=1:t_up
            C = [C, U(i,j,t) >= 0];
        end
    end    
end    

for i=1:i_up %c20: U_t0(i,j,1) is Uijt=0
    for j=1:j_up
        C = [C, U_t0(i,j,1) >= 0];
    end    
end 

for i=1:i_up %c20: Vit
    for t=1:t_up
      C = [C, V(i,t) >= 0];
    end    
end    

for i=1:i_up %c20: V_t0(i,1) is V(i,0)
  C = [C, V_t0(i,1) >= 0];
end

%1.3-solve operational level: settings: use cplex as solver

%0515-following is set optimality for MIP, but kang said make it LNContinuous-https://www.ibm.com/support/knowledgecenter/SSSA5P_12.9.0/ilog.odms.cplex.help/CPLEX/Parameters/topics/OptimalityTarget.html
%0515-'cplex.optimalitytarget'=1: Searches for a globally optimal solution to a convex model.
ops = sdpsettings ('solver','cplex','verbose',2,'cplex.optimalitytarget',1);%0508-verbose 2 shows iterations, 1 means shown normal message, 0 shows silent message
%0515-MIP limit tolerance gap to 1-https://www.ibm.com/support/knowledgecenter/SSSA5P_12.9.0/ilog.odms.cplex.help/CPLEX/Parameters/topics/EpGap.html
%ops = sdpsettings ('solver','cplex','verbose',2,'cplex.mip.tolerances.mipgap',1);%0508-verbose 2 shows iterations, 1 means shown normal message, 0 shows silent message

%0306-limit barrier iteration-https://www.ibm.com/support/knowledgecenter/SSSA5P_12.9.0/ilog.odms.cplex.help/CPLEX/Parameters/topics/BarItLim.html
%Cplex.Param.barrier.limits.iteration = 4;[not work]

%0306-limit network simplex iteration-https://www.ibm.com/support/knowledgecenter/SSSA5P_12.9.0/ilog.odms.cplex.help/CPLEX/Parameters/topics/NetItLim.html
%Cplex.Param.network.iterations = 2;[not work]

%ops.mip.limits.nodes = 4;[not work]

%0304-use the online method for out-of-memory: https://www.ilovematlab.cn/thread-266244-1-1.html
%ops.mip.strategy.file = 3; %[not work: should put in sdpsetting] this has same results 1159110 node, but reduces time for 3 min.

%0305-use online method for out-of-memory: https://www.ibm.com/developerworks/community/forums/html/topic?id=e3ea344c-7249-4edd-a47d-29e1f80fa480
%Cplex.Param.mip.limits.auxrootthreads = 2; %[not work: should put in sdpsetting] -1 or 1

%obj-operational obj z_op
z_op = sum(sum(sum(cijt.*X(1:i_up,1:j_up,1:t_up)))) + sum(sum(cijt.*X_t0(1:i_up,1:j_up,1)))...
  + sum(sum(sum(pijt.*U(1:i_up,1:j_up,1:t_up)))) + sum(sum(pijt.*U_t0(1:i_up,1:j_up,1)))...
  + sum(sum(hit.*V(1:i_up,1:t_up))) + sum(hit.*V_t0(1:i_up,1));

%1.35-express D_ij_SAV
for i=1:i_up
  for j=1:j_up
    for t=1:t_up
      D_ij_SAV(i,j,t) = (b*d(i,j,t)+cijt*d(i,j,t)*m-n*T_abs*U(i,j,t))./(a+b+m*(cijt+P(i,j,t)));
    end
  end
end

% 0604-calculate minimum k ratio covered by SAVs
for i=1:i_up
  for j=1:j_up
    for t=1:t_up
      k(i,j,t) = D_ij_SAV(i,j,t)./d(i,j,t);
    end
  end
end

% k_min = min(k(1,1,1),k(1,2,1),k(1,3,1))

% k_min = 0.99;
% % 
% if ((k_min) >= (k(1,1,1)))
%   k_min = k(1,1,1)                
% end  

% 0604-find minimum covered ratio k
% for i=1:i_up
%   for j=1:j_up
%     for t=1:t_up
%       if (k_min >= k(i,j,t))
%         k_min = k(i,j,t)                
%       end  
%     end
%   end
% end

%1.4-check if solved and output
result = optimize(C,z_op,ops); % 0508-origin
% result = solvesdp(C,z_op,ops); % 0508-this works

%result = optimize(C,z);
% if result.problem == 0 % problem =0 means solve successfully,only print Uijt here

%   fprintf('\nk value: ') %ratio k covered by SAVs
%   value(k)

  
%   fprintf('\nP_Ini_(1,1,1) = ')
%   value(P(1,1,1))

%   fprintf('\nMin cost by Lower = ')
%   value(z_op)

  all_k = [];

  for i=1:i_up
    for j=1:j_up
      for t=1:t_up
        all_k(end+1) = k(i,j,t);
      end
    end
  end
  
  min_k = min(all_k);
  
  k_avg = mean(all_k);
  
%   fprintf('\n min k: ')
%   value(min_k)
  
  max_k = max(all_k);
  
%   fprintf('\n max k: ')
%   value(max_k)
  
  % 0615-save max_k in txt

%   fid = fopen(['IniP=',num2str(2),'max_k.txt'],'a'); %save to last 'a'
%   fprintf(fid,'%.4f\t',double(max_k)); %0516-keep 4 decimals
%   fclose(fid);
  
%   fprintf('\nU(1,1,1) value: ')
%   value(U(1,1,1))
        
%   fprintf('\nX(1,1,1) value: ') 
%   value(X(1,1,1))

fprintf('\n Iteration No. = ')
value(Iter_count)



% fprintf('\n D_ij_SAV(1,1,1) value: ') %0515-as for now hide D_ij_SAV
% value(D_ij_SAV(1,1,1))
  
%   fprintf('\n k(1,1,1) ')
%   value(k(1,1,1))

%   fprintf('\n')
%     fprintf('\nV value: ')
%     value(V)
%     
%     fprintf('\nU_t0 value: ')
%     value(U_t0)
%     
%     fprintf('\nX_t0 value: ') 
%     value(X_t0)
%     
%     fprintf('\nV_t0 value: ')
%     value(V_t0)

% else
%   disp('Not solved!')
% 
% end

% 0519-save total cost z_op in txt
% fid = fopen(['Iter=',num2str(Iter_count),'_OperCost.txt'],'wt'); %creat & write z_op to txt
% fprintf(fid,'%.4f\n',double(z_op)); %0516-keep 4 decimals
% fclose(fid);

% 0605-save min ratio k covered by SAVs in txt
% fid_ratiok = fopen(['Iter=',num2str(Iter_count),'_MinRatiok.txt'],'wt'); %creat & write z_op to txt
% fprintf(fid_ratiok,'%.4f\n',double(min_k)); %0604-keep 4 decimals
% fclose(fid_ratiok);

% 0517-following is save Uijt in txt
%1.6-save Uijt as txt for upper use
% for U_t=1:t_up 
%      fid = fopen(['IniP=',num2str(IniP),'Iter=',num2str(Iter_count),'_Uijt=',num2str(U_t),'.txt'],'wt'); %creat & write Uijt to txt
%      for i=1:i_up
%          for j=1:j_up
%              if j == j_up
%                  fprintf(fid,'%.4f\n',double(U(i,j,U_t))); %0516-keep 4 decimals
%              else
%                  fprintf(fid,'%.4f\t',double(U(i,j,U_t))); %0516-keep 4 decimals
%              end     
%          end    
%      end    
%      fclose(fid);
% end

% 0607-calculate Upper Pijt upper bound
for i=1:i_up
  for j=1:j_up
    for t=1:t_up
      P(i,j,t)=(b*d(i,j,t)+cijt*d(i,j,t)*m-n*T_abs*U(i,j,t)-(a+b)*min_k*d(i,j,t))./(m*min_k*d(i,j,t))-cijt;
    end
  end
end

% fprintf('\nUpper bound P(1,1,1) = %.4f\n',double(P(1,1,1)));

%0608-----
%calculate Upper revenue
%0608-----

revenue=0;

for i=1:i_up
  for j=1:j_up
    for t=1:t_up
      if (t==1)||(t==2)||(t==9)||(t==10)||(t==17)||(t==18)||(t==31)||(t==32)||(t==39)||(t==40)
%         T_abs=2;
        revenue = revenue-(P(i,j,t).*(b.*d(i,j,t)+cijt.*d(i,j,t)*m-n*T_abs*U(i,j,t)))./(a+b+m.*(cijt+P(i,j,t))); 
      else
%         T_abs=2; %here should be 3, but temporarily set to 2, aligh with Oper level
        revenue = revenue-(P(i,j,t).*(b.*d(i,j,t)+cijt.*d(i,j,t)*m-n*T_abs*U(i,j,t)))./(a+b+m.*(cijt+P(i,j,t)));
      end
    end
  end
end

% 070320-avg. all P

 all_P = [];

  for i=1:i_up
    for j=1:j_up
      for t=1:t_up
        all_P = P(i,j,t);
      end
    end
  end

avg_P = all_P./(i_up*j_up*t_up);



% 070320-avg. all D

 all_Dijt = [];

  for i=1:i_up
    for j=1:j_up
      for t=1:t_up
        all_Dijt = D_ij_SAV(i,j,t);
      end
    end
  end

avg_Dijt = all_Dijt./(i_up*j_up*t_up);

% 070320-print results

fprintf('\n Upper bound P(1,1,1) = %.4f\n',double(P(1,1,1))); % particular P each iteration

fprintf('\n');

fprintf('\n D(1,1,1) = %.4f\n',double(D_ij_SAV(1,1,1))); % particular Dijt each iteration

fprintf('\n');

fprintf('\n Avg_P = %.4f\n',double(avg_P)); % particular avg_P each iteration

fprintf('\n');

fprintf('\n Avg_Dijt = %.4f\n',double(avg_Dijt)); % particular avg_Dijt each iteration

fprintf('\n');

fprintf('\n Max Revenue by Upper= %.4f\n',double(-revenue));

fprintf('\n');

fprintf('\n Min Cost by Lower= %.4f\n',double(z_op));

fprintf('\n');

fprintf('\n Profit = %.4f\n',double(-revenue-z_op));

fprintf('\n');

% test how many new P > old P

satis_P = 0;

for i=1:i_up
  for j=1:j_up
    for t=1:t_up
%       if (double(P_sav(i,j,t))<=double(P(i,j,t))) % 0610-check new P greater than old ones
        %0610-if new P no more than 1 cent than old, than add to satis_P terminate condition
      %if ((double(P_sav(i,j,t)))*1.0001 >= double(P(i,j,t))) 
%       if ((double(P_sav(i,j,t)))+0.01 >= double(P(i,j,t))) 
      if ((P_sav(i,j,t))+0.01 > P(i,j,t)) 
        satis_P = (satis_P)+1;
      end  
    end
  end
end

satis_P;

fprintf('\n');

% test how many new P all feasible that new_P <= new P upper bound condition


if (Iter_count==1)
  feasib_P = 0;
else
  
  feasib_P = 0;
  
  for i=1:i_up
    for j=1:j_up
      for t=1:t_up
        P_bench(i,j,t) = (b*d(i,j,t)+cijt*d(i,j,t)*m-n*T_abs*U(i,j,t)-(a+b)*min_k*d(i,j,t))./(m*min_k*d(i,j,t))-cijt;        
      end
    end
  end
  
    
  for i=1:i_up
    for j=1:j_up
      for t=1:t_up
  %       if (double(P_sav(i,j,t))<=double(P(i,j,t))) % 0610-check new P greater than old ones
          %0610-if new P no more than 1 cent than old, than add to satis_P terminate condition
        %if ((double(P_sav(i,j,t)))*1.0001 >= double(P(i,j,t))) 
        
%         if (double(P_sav(i,j,t)) <= double(P_bench(i,j,t)))
        if (P_bench(i,j,t) >= P_sav(i,j,t))
%           i,j,t
%           fprintf("\n");
          feasib_P = (feasib_P)+1;
        end  
      end
    end
  end

end

feasib_P;

fprintf('\n');

if (satis_P==i_up*j_up*t_up) && (feasib_P==i_up*j_up*t_up)
  consec_P = 1+consec_P;
else
  consec_P = 0;
end

% tik = tik-1;

end

fprintf('\n');
fprintf('\n');
fprintf('\n');


%-------case 3.1------

'------Case 3.1------'

IniP=1;

%initially price Pijt all the same
P=zeros(20,20,40);%create 20*20*40 matrix
for i=1:i_up 
    for j=1:j_up
        for t=1:t_up
          P(i,j,t)=IniP; %0608-initial P should be infinity, I use 2 for test
        end
    end    
end

% Pijt=2;

Iter_count = 0;

% 0517-M is c16 time window, set at 1 gives it relaxation
% 0424-M is c16 time window, set at 3 gives it relaxation
M=1; 

%0515-change all var to continours instead of MIP
X = sdpvar(i_up,j_up,t_up+M); %0424-Xi,j,t+M, cuz time window constraint c16 Xijt is beyond t
%X = intvar(i_up,j_up,t_up+M); %0424-Xi,j,t+M, cuz time window constraint c16 Xijt is beyond t
X_t0 = sdpvar(i_up,j_up,1); %0424-Xi,j,0 means moves initially for c18: inventory balance eq
%X_t0 = intvar(i_up,j_up,1); %0424-Xi,j,0 means moves initially for c18: inventory balance eq
U = sdpvar(i_up,j_up,t_up);
%U = intvar(i_up,j_up,t_up); 
U_t0 = sdpvar(i_up,j_up,1); %0424-this is for Ui,j,0 
%U_t0 = intvar(i_up,j_up,1); %0424-this is for Ui,j,0 
V = sdpvar(i_up,t_up); 
%V = intvar(i_up,t_up); 
V_t0 = sdpvar(i_up,1); %0424-this is especially for Vi,0
%V_t0 = intvar(i_up,1); %0424-this is especially for Vi,0
D_ij_SAV = sdpvar(i_up,j_up,t_up); %0425-this is to show D_ij_SAV, for futher plot, then find upper bound of Pijt
%D_ij_SAV = intvar(i_up,j_up,t_up); %0425-this is to show D_ij_SAV, for futher plot, then find upper bound of Pijt
% 0604-k ratio
k = sdpvar(i_up,j_up,t_up);

% 070220-set place for Wijt
W=zeros(20,20,40);%create 20*20*40 matrix

%1.1-parameters
% a=b=m=n=1,T=40
a=1;
b=1;
m=10;
n=1;
T=8; %070120-total 40 time periods, 5 time block, so every block is 8
T_abs=8; %070120-T_abs should be same as T
% cijt=40, pijt=10, hit=10, alphaijt=0.1 (all c,p,h,alpha are subjectively to be same)
cijt=20;
pijt=10; % 070120-increase for penalty of unmet demand
% 0425-intentionally make pijt this small to avoid all 0 in Uijt solution
hit=10;
alpha=1;

%%%-----0607----
%%% How to add loop?
%%%----

% 0610-here starts the loop

P_sav = zeros(20,20,40);
U_sav = zeros(20,20,40);
min_k = 0;
P_bench = zeros(20,20,40);

% tik = 2;

satis_P = 0;
feasib_P = 0;
consec_P = 0;

% while (satis_P<i_up*j_up*t_up) || (feasib_P<i_up*j_up*t_up)
% while Iter_count < 2
while (consec_P<3)

Iter_count = Iter_count+1; 


for i=1:i_up 
    for j=1:j_up
        for t=1:t_up
          P_sav(i,j,t)=P(i,j,t); %0608-initial P should be infinity, I use 2 for test
        end
    end    
end  


for i=1:i_up 
    for j=1:j_up
        for t=1:t_up
            U_sav(i,j,t)=U(i,j,t); 
        end
    end    
end


min_k_sav = min_k;


%1.2-add operational level constraints

C = [];

% c15
for i=1:i_up %c15
   for j=1:j_up
     for t=1:t_up
       if (t==1)
         C = [C, U(i,j,1) >= U_t0(i,j,1)+(b*d(i,j,1)+cijt*d(i,j,1)*m-n*T_abs*U(i,j,1))./(a+b+m*(cijt+P(i,j,t)))-X(i,j,1)]; %c15: t=1
       else       
         C = [C, U(i,j,t) >= U(i,j,t-1)+(b*d(i,j,t)+cijt*d(i,j,t)*m-n*T_abs*U(i,j,t))./(a+b+m*(cijt+P(i,j,t)))-X(i,j,t)]; %c15: t from 2
       end
     end
   end    
end

% 0424-c16-I suppose time window is M period for c16 - relax this one, hope get solution
% this also considers Xijt beyond T
for i=1:i_up %c16
    for j=1:j_up
      for t=1:t_up
        if (t==1)
          C = [C, sum(X(i,j,t:(t+M))) >= U_t0(i,j,1)+(b*d(i,j,t)+cijt*d(i,j,t)*m-n*T_abs*U(i,j,t))./(a+b+m*(cijt+P(i,j,t)))]; %c16: t=1
        else
          C = [C, sum(X(i,j,t:(t+M))) >= U(i,j,t-1)+(b*d(i,j,t)+cijt*d(i,j,t)*m-n*T_abs*U(i,j,t))./(a+b+m*(cijt+P(i,j,t)))]; %c16: t from 2
        end
      end  
    end     
end

% 0424 - c17 is in upper level, no need in lower level
% for i=1:i_up %c17
%   for j=1:j_up
%     for t=1:t_up
%       C = [C, ((b-1)*d(i,j,t)+cijt*d(i,j,t)*m)./(n*T_abs) <= U(i,j,t)];
%       C = [C, U(i,j,t) <= (b*d(i,j,t)+cijt*d(i,j,t)*m)./(n*T_abs)];
%     end
%   end
% end

% 0424-c18-here makes <= for Vi,t cuz alpha all set to 1, too relaxed, should be variaty
for i=1:i_up %c18
  for t=1:t_up
    if (t==1)
      C = [C, V(i,t) <= V_t0(i,1) + sum(alpha*X_t0(1:j_up,i,1)) - sum(X(i,1:j_up,t)) ]; %c18: t=1
    else
      C = [C, V(i,t) <= V(i,t-1) + sum(sum(alpha*X(1:j_up,i,1:(t-1)))) - sum(X(i,1:j_up,t)) ]; %c18: t from 2 (remember to use 2 sum for 2nd)
    end
  end
end

% 0508-c19-not added, cuz it looks already relaxed in math function, hidden property of our formula

for i=1:i_up %c20: Xijt
    for j=1:j_up
        for t=1:(t_up+M)
          C = [C, X(i,j,t) >= 0];
        end
    end    
end    

for i=1:i_up %c20: X_t0(i,j,1) is Xijt=0
    for j=1:j_up
        C = [C, X_t0(i,j,1) >= 0];
    end    
end    

for i=1:i_up %c20: Uijt
    for j=1:j_up
        for t=1:t_up
            C = [C, U(i,j,t) >= 0];
        end
    end    
end    

for i=1:i_up %c20: U_t0(i,j,1) is Uijt=0
    for j=1:j_up
        C = [C, U_t0(i,j,1) >= 0];
    end    
end 

for i=1:i_up %c20: Vit
    for t=1:t_up
      C = [C, V(i,t) >= 0];
    end    
end    

for i=1:i_up %c20: V_t0(i,1) is V(i,0)
  C = [C, V_t0(i,1) >= 0];
end

%1.3-solve operational level: settings: use cplex as solver

%0515-following is set optimality for MIP, but kang said make it LNContinuous-https://www.ibm.com/support/knowledgecenter/SSSA5P_12.9.0/ilog.odms.cplex.help/CPLEX/Parameters/topics/OptimalityTarget.html
%0515-'cplex.optimalitytarget'=1: Searches for a globally optimal solution to a convex model.
ops = sdpsettings ('solver','cplex','verbose',2,'cplex.optimalitytarget',1);%0508-verbose 2 shows iterations, 1 means shown normal message, 0 shows silent message
%0515-MIP limit tolerance gap to 1-https://www.ibm.com/support/knowledgecenter/SSSA5P_12.9.0/ilog.odms.cplex.help/CPLEX/Parameters/topics/EpGap.html
%ops = sdpsettings ('solver','cplex','verbose',2,'cplex.mip.tolerances.mipgap',1);%0508-verbose 2 shows iterations, 1 means shown normal message, 0 shows silent message

%0306-limit barrier iteration-https://www.ibm.com/support/knowledgecenter/SSSA5P_12.9.0/ilog.odms.cplex.help/CPLEX/Parameters/topics/BarItLim.html
%Cplex.Param.barrier.limits.iteration = 4;[not work]

%0306-limit network simplex iteration-https://www.ibm.com/support/knowledgecenter/SSSA5P_12.9.0/ilog.odms.cplex.help/CPLEX/Parameters/topics/NetItLim.html
%Cplex.Param.network.iterations = 2;[not work]

%ops.mip.limits.nodes = 4;[not work]

%0304-use the online method for out-of-memory: https://www.ilovematlab.cn/thread-266244-1-1.html
%ops.mip.strategy.file = 3; %[not work: should put in sdpsetting] this has same results 1159110 node, but reduces time for 3 min.

%0305-use online method for out-of-memory: https://www.ibm.com/developerworks/community/forums/html/topic?id=e3ea344c-7249-4edd-a47d-29e1f80fa480
%Cplex.Param.mip.limits.auxrootthreads = 2; %[not work: should put in sdpsetting] -1 or 1

%obj-operational obj z_op
z_op = sum(sum(sum(cijt.*X(1:i_up,1:j_up,1:t_up)))) + sum(sum(cijt.*X_t0(1:i_up,1:j_up,1)))...
  + sum(sum(sum(pijt.*U(1:i_up,1:j_up,1:t_up)))) + sum(sum(pijt.*U_t0(1:i_up,1:j_up,1)))...
  + sum(sum(hit.*V(1:i_up,1:t_up))) + sum(hit.*V_t0(1:i_up,1));

%1.35-express D_ij_SAV
for i=1:i_up
  for j=1:j_up
    for t=1:t_up
      D_ij_SAV(i,j,t) = (b*d(i,j,t)+cijt*d(i,j,t)*m-n*T_abs*U(i,j,t))./(a+b+m*(cijt+P(i,j,t)));
    end
  end
end

% 0604-calculate minimum k ratio covered by SAVs
for i=1:i_up
  for j=1:j_up
    for t=1:t_up
      k(i,j,t) = D_ij_SAV(i,j,t)./d(i,j,t);
    end
  end
end

% k_min = min(k(1,1,1),k(1,2,1),k(1,3,1))

% k_min = 0.99;
% % 
% if ((k_min) >= (k(1,1,1)))
%   k_min = k(1,1,1)                
% end  

% 0604-find minimum covered ratio k
% for i=1:i_up
%   for j=1:j_up
%     for t=1:t_up
%       if (k_min >= k(i,j,t))
%         k_min = k(i,j,t)                
%       end  
%     end
%   end
% end

%1.4-check if solved and output
result = optimize(C,z_op,ops); % 0508-origin
% result = solvesdp(C,z_op,ops); % 0508-this works

%result = optimize(C,z);
% if result.problem == 0 % problem =0 means solve successfully,only print Uijt here

%   fprintf('\nk value: ') %ratio k covered by SAVs
%   value(k)

  
%   fprintf('\nP_Ini_(1,1,1) = ')
%   value(P(1,1,1))

%   fprintf('\nMin cost by Lower = ')
%   value(z_op)

  all_k = [];

  for i=1:i_up
    for j=1:j_up
      for t=1:t_up
        all_k(end+1) = k(i,j,t);
      end
    end
  end
  
  min_k = min(all_k);
  
  k_avg = mean(all_k);
  
%   fprintf('\n min k: ')
%   value(min_k)
  
  max_k = max(all_k);
  
%   fprintf('\n max k: ')
%   value(max_k)
  
  % 0615-save max_k in txt

%   fid = fopen(['IniP=',num2str(2),'max_k.txt'],'a'); %save to last 'a'
%   fprintf(fid,'%.4f\t',double(max_k)); %0516-keep 4 decimals
%   fclose(fid);
  
%   fprintf('\nU(1,1,1) value: ')
%   value(U(1,1,1))
        
%   fprintf('\nX(1,1,1) value: ') 
%   value(X(1,1,1))

fprintf('\n Iteration No. = ')
value(Iter_count)



% fprintf('\n D_ij_SAV(1,1,1) value: ') %0515-as for now hide D_ij_SAV
% value(D_ij_SAV(1,1,1))
  
%   fprintf('\n k(1,1,1) ')
%   value(k(1,1,1))

%   fprintf('\n')
%     fprintf('\nV value: ')
%     value(V)
%     
%     fprintf('\nU_t0 value: ')
%     value(U_t0)
%     
%     fprintf('\nX_t0 value: ') 
%     value(X_t0)
%     
%     fprintf('\nV_t0 value: ')
%     value(V_t0)

% else
%   disp('Not solved!')
% 
% end

% 0519-save total cost z_op in txt
% fid = fopen(['Iter=',num2str(Iter_count),'_OperCost.txt'],'wt'); %creat & write z_op to txt
% fprintf(fid,'%.4f\n',double(z_op)); %0516-keep 4 decimals
% fclose(fid);

% 0605-save min ratio k covered by SAVs in txt
% fid_ratiok = fopen(['Iter=',num2str(Iter_count),'_MinRatiok.txt'],'wt'); %creat & write z_op to txt
% fprintf(fid_ratiok,'%.4f\n',double(min_k)); %0604-keep 4 decimals
% fclose(fid_ratiok);

% 0517-following is save Uijt in txt
%1.6-save Uijt as txt for upper use
% for U_t=1:t_up 
%      fid = fopen(['IniP=',num2str(IniP),'Iter=',num2str(Iter_count),'_Uijt=',num2str(U_t),'.txt'],'wt'); %creat & write Uijt to txt
%      for i=1:i_up
%          for j=1:j_up
%              if j == j_up
%                  fprintf(fid,'%.4f\n',double(U(i,j,U_t))); %0516-keep 4 decimals
%              else
%                  fprintf(fid,'%.4f\t',double(U(i,j,U_t))); %0516-keep 4 decimals
%              end     
%          end    
%      end    
%      fclose(fid);
% end

% 0607-calculate Upper Pijt upper bound
for i=1:i_up
  for j=1:j_up
    for t=1:t_up
      P(i,j,t)=(b*d(i,j,t)+cijt*d(i,j,t)*m-n*T_abs*U(i,j,t)-(a+b)*min_k*d(i,j,t))./(m*min_k*d(i,j,t))-cijt;
    end
  end
end

% fprintf('\nUpper bound P(1,1,1) = %.4f\n',double(P(1,1,1)));

%0608-----
%calculate Upper revenue
%0608-----

revenue=0;

for i=1:i_up
  for j=1:j_up
    for t=1:t_up
      if (t==1)||(t==2)||(t==9)||(t==10)||(t==17)||(t==18)||(t==31)||(t==32)||(t==39)||(t==40)
%         T_abs=2;
        revenue = revenue-(P(i,j,t).*(b.*d(i,j,t)+cijt.*d(i,j,t)*m-n*T_abs*U(i,j,t)))./(a+b+m.*(cijt+P(i,j,t))); 
      else
%         T_abs=2; %here should be 3, but temporarily set to 2, aligh with Oper level
        revenue = revenue-(P(i,j,t).*(b.*d(i,j,t)+cijt.*d(i,j,t)*m-n*T_abs*U(i,j,t)))./(a+b+m.*(cijt+P(i,j,t)));
      end
    end
  end
end

% 070320-avg. all P

 all_P = [];

  for i=1:i_up
    for j=1:j_up
      for t=1:t_up
        all_P = P(i,j,t);
      end
    end
  end

avg_P = all_P./(i_up*j_up*t_up);



% 070320-avg. all D

 all_Dijt = [];

  for i=1:i_up
    for j=1:j_up
      for t=1:t_up
        all_Dijt = D_ij_SAV(i,j,t);
      end
    end
  end

avg_Dijt = all_Dijt./(i_up*j_up*t_up);

% 070320-print results

fprintf('\n Upper bound P(1,1,1) = %.4f\n',double(P(1,1,1))); % particular P each iteration

fprintf('\n');

fprintf('\n D(1,1,1) = %.4f\n',double(D_ij_SAV(1,1,1))); % particular Dijt each iteration

fprintf('\n');

fprintf('\n Avg_P = %.4f\n',double(avg_P)); % particular avg_P each iteration

fprintf('\n');

fprintf('\n Avg_Dijt = %.4f\n',double(avg_Dijt)); % particular avg_Dijt each iteration

fprintf('\n');

fprintf('\n Max Revenue by Upper= %.4f\n',double(-revenue));

fprintf('\n');

fprintf('\n Min Cost by Lower= %.4f\n',double(z_op));

fprintf('\n');

fprintf('\n Profit = %.4f\n',double(-revenue-z_op));

fprintf('\n');

% test how many new P > old P

satis_P = 0;

for i=1:i_up
  for j=1:j_up
    for t=1:t_up
%       if (double(P_sav(i,j,t))<=double(P(i,j,t))) % 0610-check new P greater than old ones
        %0610-if new P no more than 1 cent than old, than add to satis_P terminate condition
      %if ((double(P_sav(i,j,t)))*1.0001 >= double(P(i,j,t))) 
%       if ((double(P_sav(i,j,t)))+0.01 >= double(P(i,j,t))) 
      if ((P_sav(i,j,t))+0.01 > P(i,j,t)) 
        satis_P = (satis_P)+1;
      end  
    end
  end
end

satis_P;

fprintf('\n');

% test how many new P all feasible that new_P <= new P upper bound condition


if (Iter_count==1)
  feasib_P = 0;
else
  
  feasib_P = 0;
  
  for i=1:i_up
    for j=1:j_up
      for t=1:t_up
        P_bench(i,j,t) = (b*d(i,j,t)+cijt*d(i,j,t)*m-n*T_abs*U(i,j,t)-(a+b)*min_k*d(i,j,t))./(m*min_k*d(i,j,t))-cijt;        
      end
    end
  end
  
    
  for i=1:i_up
    for j=1:j_up
      for t=1:t_up
  %       if (double(P_sav(i,j,t))<=double(P(i,j,t))) % 0610-check new P greater than old ones
          %0610-if new P no more than 1 cent than old, than add to satis_P terminate condition
        %if ((double(P_sav(i,j,t)))*1.0001 >= double(P(i,j,t))) 
        
%         if (double(P_sav(i,j,t)) <= double(P_bench(i,j,t)))
        if (P_bench(i,j,t) >= P_sav(i,j,t))
%           i,j,t
%           fprintf("\n");
          feasib_P = (feasib_P)+1;
        end  
      end
    end
  end

end

feasib_P;

fprintf('\n');

if (satis_P==i_up*j_up*t_up) && (feasib_P==i_up*j_up*t_up)
  consec_P = 1+consec_P;
else
  consec_P = 0;
end

% tik = tik-1;

end

fprintf('\n');
fprintf('\n');
fprintf('\n');


%-------case 3.2------

'-------Case 3.2-------'

IniP=10;

%initially price Pijt all the same
P=zeros(20,20,40);%create 20*20*40 matrix
for i=1:i_up 
    for j=1:j_up
        for t=1:t_up
          P(i,j,t)=IniP; %0608-initial P should be infinity, I use 2 for test
        end
    end    
end

% Pijt=2;

Iter_count = 0;

% 0517-M is c16 time window, set at 1 gives it relaxation
% 0424-M is c16 time window, set at 3 gives it relaxation
M=1; 

%0515-change all var to continours instead of MIP
X = sdpvar(i_up,j_up,t_up+M); %0424-Xi,j,t+M, cuz time window constraint c16 Xijt is beyond t
%X = intvar(i_up,j_up,t_up+M); %0424-Xi,j,t+M, cuz time window constraint c16 Xijt is beyond t
X_t0 = sdpvar(i_up,j_up,1); %0424-Xi,j,0 means moves initially for c18: inventory balance eq
%X_t0 = intvar(i_up,j_up,1); %0424-Xi,j,0 means moves initially for c18: inventory balance eq
U = sdpvar(i_up,j_up,t_up);
%U = intvar(i_up,j_up,t_up); 
U_t0 = sdpvar(i_up,j_up,1); %0424-this is for Ui,j,0 
%U_t0 = intvar(i_up,j_up,1); %0424-this is for Ui,j,0 
V = sdpvar(i_up,t_up); 
%V = intvar(i_up,t_up); 
V_t0 = sdpvar(i_up,1); %0424-this is especially for Vi,0
%V_t0 = intvar(i_up,1); %0424-this is especially for Vi,0
D_ij_SAV = sdpvar(i_up,j_up,t_up); %0425-this is to show D_ij_SAV, for futher plot, then find upper bound of Pijt
%D_ij_SAV = intvar(i_up,j_up,t_up); %0425-this is to show D_ij_SAV, for futher plot, then find upper bound of Pijt
% 0604-k ratio
k = sdpvar(i_up,j_up,t_up);

% 070220-set place for Wijt
W=zeros(20,20,40);%create 20*20*40 matrix

%1.1-parameters
% a=b=m=n=1,T=40
a=1;
b=1;
m=10;
n=1;
T=8; %070120-total 40 time periods, 5 time block, so every block is 8
T_abs=8; %070120-T_abs should be same as T
% cijt=40, pijt=10, hit=10, alphaijt=0.1 (all c,p,h,alpha are subjectively to be same)
cijt=20;
pijt=10; % 070120-increase for penalty of unmet demand
% 0425-intentionally make pijt this small to avoid all 0 in Uijt solution
hit=10;
alpha=1;

%%%-----0607----
%%% How to add loop?
%%%----

% 0610-here starts the loop

P_sav = zeros(20,20,40);
U_sav = zeros(20,20,40);
min_k = 0;
P_bench = zeros(20,20,40);

% tik = 2;

satis_P = 0;
feasib_P = 0;
consec_P = 0;

% while (satis_P<i_up*j_up*t_up) || (feasib_P<i_up*j_up*t_up)
% while Iter_count < 2
while (consec_P<3)

Iter_count = Iter_count+1; 


for i=1:i_up 
    for j=1:j_up
        for t=1:t_up
          P_sav(i,j,t)=P(i,j,t); %0608-initial P should be infinity, I use 2 for test
        end
    end    
end  


for i=1:i_up 
    for j=1:j_up
        for t=1:t_up
            U_sav(i,j,t)=U(i,j,t); 
        end
    end    
end


min_k_sav = min_k;


%1.2-add operational level constraints

C = [];

% c15
for i=1:i_up %c15
   for j=1:j_up
     for t=1:t_up
       if (t==1)
         C = [C, U(i,j,1) >= U_t0(i,j,1)+(b*d(i,j,1)+cijt*d(i,j,1)*m-n*T_abs*U(i,j,1))./(a+b+m*(cijt+P(i,j,t)))-X(i,j,1)]; %c15: t=1
       else       
         C = [C, U(i,j,t) >= U(i,j,t-1)+(b*d(i,j,t)+cijt*d(i,j,t)*m-n*T_abs*U(i,j,t))./(a+b+m*(cijt+P(i,j,t)))-X(i,j,t)]; %c15: t from 2
       end
     end
   end    
end

% 0424-c16-I suppose time window is M period for c16 - relax this one, hope get solution
% this also considers Xijt beyond T
for i=1:i_up %c16
    for j=1:j_up
      for t=1:t_up
        if (t==1)
          C = [C, sum(X(i,j,t:(t+M))) >= U_t0(i,j,1)+(b*d(i,j,t)+cijt*d(i,j,t)*m-n*T_abs*U(i,j,t))./(a+b+m*(cijt+P(i,j,t)))]; %c16: t=1
        else
          C = [C, sum(X(i,j,t:(t+M))) >= U(i,j,t-1)+(b*d(i,j,t)+cijt*d(i,j,t)*m-n*T_abs*U(i,j,t))./(a+b+m*(cijt+P(i,j,t)))]; %c16: t from 2
        end
      end  
    end     
end

% 0424 - c17 is in upper level, no need in lower level
% for i=1:i_up %c17
%   for j=1:j_up
%     for t=1:t_up
%       C = [C, ((b-1)*d(i,j,t)+cijt*d(i,j,t)*m)./(n*T_abs) <= U(i,j,t)];
%       C = [C, U(i,j,t) <= (b*d(i,j,t)+cijt*d(i,j,t)*m)./(n*T_abs)];
%     end
%   end
% end

% 0424-c18-here makes <= for Vi,t cuz alpha all set to 1, too relaxed, should be variaty
for i=1:i_up %c18
  for t=1:t_up
    if (t==1)
      C = [C, V(i,t) <= V_t0(i,1) + sum(alpha*X_t0(1:j_up,i,1)) - sum(X(i,1:j_up,t)) ]; %c18: t=1
    else
      C = [C, V(i,t) <= V(i,t-1) + sum(sum(alpha*X(1:j_up,i,1:(t-1)))) - sum(X(i,1:j_up,t)) ]; %c18: t from 2 (remember to use 2 sum for 2nd)
    end
  end
end

% 0508-c19-not added, cuz it looks already relaxed in math function, hidden property of our formula

for i=1:i_up %c20: Xijt
    for j=1:j_up
        for t=1:(t_up+M)
          C = [C, X(i,j,t) >= 0];
        end
    end    
end    

for i=1:i_up %c20: X_t0(i,j,1) is Xijt=0
    for j=1:j_up
        C = [C, X_t0(i,j,1) >= 0];
    end    
end    

for i=1:i_up %c20: Uijt
    for j=1:j_up
        for t=1:t_up
            C = [C, U(i,j,t) >= 0];
        end
    end    
end    

for i=1:i_up %c20: U_t0(i,j,1) is Uijt=0
    for j=1:j_up
        C = [C, U_t0(i,j,1) >= 0];
    end    
end 

for i=1:i_up %c20: Vit
    for t=1:t_up
      C = [C, V(i,t) >= 0];
    end    
end    

for i=1:i_up %c20: V_t0(i,1) is V(i,0)
  C = [C, V_t0(i,1) >= 0];
end

%1.3-solve operational level: settings: use cplex as solver

%0515-following is set optimality for MIP, but kang said make it LNContinuous-https://www.ibm.com/support/knowledgecenter/SSSA5P_12.9.0/ilog.odms.cplex.help/CPLEX/Parameters/topics/OptimalityTarget.html
%0515-'cplex.optimalitytarget'=1: Searches for a globally optimal solution to a convex model.
ops = sdpsettings ('solver','cplex','verbose',2,'cplex.optimalitytarget',1);%0508-verbose 2 shows iterations, 1 means shown normal message, 0 shows silent message
%0515-MIP limit tolerance gap to 1-https://www.ibm.com/support/knowledgecenter/SSSA5P_12.9.0/ilog.odms.cplex.help/CPLEX/Parameters/topics/EpGap.html
%ops = sdpsettings ('solver','cplex','verbose',2,'cplex.mip.tolerances.mipgap',1);%0508-verbose 2 shows iterations, 1 means shown normal message, 0 shows silent message

%0306-limit barrier iteration-https://www.ibm.com/support/knowledgecenter/SSSA5P_12.9.0/ilog.odms.cplex.help/CPLEX/Parameters/topics/BarItLim.html
%Cplex.Param.barrier.limits.iteration = 4;[not work]

%0306-limit network simplex iteration-https://www.ibm.com/support/knowledgecenter/SSSA5P_12.9.0/ilog.odms.cplex.help/CPLEX/Parameters/topics/NetItLim.html
%Cplex.Param.network.iterations = 2;[not work]

%ops.mip.limits.nodes = 4;[not work]

%0304-use the online method for out-of-memory: https://www.ilovematlab.cn/thread-266244-1-1.html
%ops.mip.strategy.file = 3; %[not work: should put in sdpsetting] this has same results 1159110 node, but reduces time for 3 min.

%0305-use online method for out-of-memory: https://www.ibm.com/developerworks/community/forums/html/topic?id=e3ea344c-7249-4edd-a47d-29e1f80fa480
%Cplex.Param.mip.limits.auxrootthreads = 2; %[not work: should put in sdpsetting] -1 or 1

%obj-operational obj z_op
z_op = sum(sum(sum(cijt.*X(1:i_up,1:j_up,1:t_up)))) + sum(sum(cijt.*X_t0(1:i_up,1:j_up,1)))...
  + sum(sum(sum(pijt.*U(1:i_up,1:j_up,1:t_up)))) + sum(sum(pijt.*U_t0(1:i_up,1:j_up,1)))...
  + sum(sum(hit.*V(1:i_up,1:t_up))) + sum(hit.*V_t0(1:i_up,1));

%1.35-express D_ij_SAV
for i=1:i_up
  for j=1:j_up
    for t=1:t_up
      D_ij_SAV(i,j,t) = (b*d(i,j,t)+cijt*d(i,j,t)*m-n*T_abs*U(i,j,t))./(a+b+m*(cijt+P(i,j,t)));
    end
  end
end

% 0604-calculate minimum k ratio covered by SAVs
for i=1:i_up
  for j=1:j_up
    for t=1:t_up
      k(i,j,t) = D_ij_SAV(i,j,t)./d(i,j,t);
    end
  end
end

% k_min = min(k(1,1,1),k(1,2,1),k(1,3,1))

% k_min = 0.99;
% % 
% if ((k_min) >= (k(1,1,1)))
%   k_min = k(1,1,1)                
% end  

% 0604-find minimum covered ratio k
% for i=1:i_up
%   for j=1:j_up
%     for t=1:t_up
%       if (k_min >= k(i,j,t))
%         k_min = k(i,j,t)                
%       end  
%     end
%   end
% end

%1.4-check if solved and output
result = optimize(C,z_op,ops); % 0508-origin
% result = solvesdp(C,z_op,ops); % 0508-this works

%result = optimize(C,z);
% if result.problem == 0 % problem =0 means solve successfully,only print Uijt here

%   fprintf('\nk value: ') %ratio k covered by SAVs
%   value(k)

  
%   fprintf('\nP_Ini_(1,1,1) = ')
%   value(P(1,1,1))

%   fprintf('\nMin cost by Lower = ')
%   value(z_op)

  all_k = [];

  for i=1:i_up
    for j=1:j_up
      for t=1:t_up
        all_k(end+1) = k(i,j,t);
      end
    end
  end
  
  min_k = min(all_k);
  
  k_avg = mean(all_k);
  
%   fprintf('\n min k: ')
%   value(min_k)
  
  max_k = max(all_k);
  
%   fprintf('\n max k: ')
%   value(max_k)
  
  % 0615-save max_k in txt

%   fid = fopen(['IniP=',num2str(2),'max_k.txt'],'a'); %save to last 'a'
%   fprintf(fid,'%.4f\t',double(max_k)); %0516-keep 4 decimals
%   fclose(fid);
  
%   fprintf('\nU(1,1,1) value: ')
%   value(U(1,1,1))
        
%   fprintf('\nX(1,1,1) value: ') 
%   value(X(1,1,1))

fprintf('\n Iteration No. = ')
value(Iter_count)



% fprintf('\n D_ij_SAV(1,1,1) value: ') %0515-as for now hide D_ij_SAV
% value(D_ij_SAV(1,1,1))
  
%   fprintf('\n k(1,1,1) ')
%   value(k(1,1,1))

%   fprintf('\n')
%     fprintf('\nV value: ')
%     value(V)
%     
%     fprintf('\nU_t0 value: ')
%     value(U_t0)
%     
%     fprintf('\nX_t0 value: ') 
%     value(X_t0)
%     
%     fprintf('\nV_t0 value: ')
%     value(V_t0)

% else
%   disp('Not solved!')
% 
% end

% 0519-save total cost z_op in txt
% fid = fopen(['Iter=',num2str(Iter_count),'_OperCost.txt'],'wt'); %creat & write z_op to txt
% fprintf(fid,'%.4f\n',double(z_op)); %0516-keep 4 decimals
% fclose(fid);

% 0605-save min ratio k covered by SAVs in txt
% fid_ratiok = fopen(['Iter=',num2str(Iter_count),'_MinRatiok.txt'],'wt'); %creat & write z_op to txt
% fprintf(fid_ratiok,'%.4f\n',double(min_k)); %0604-keep 4 decimals
% fclose(fid_ratiok);

% 0517-following is save Uijt in txt
%1.6-save Uijt as txt for upper use
% for U_t=1:t_up 
%      fid = fopen(['IniP=',num2str(IniP),'Iter=',num2str(Iter_count),'_Uijt=',num2str(U_t),'.txt'],'wt'); %creat & write Uijt to txt
%      for i=1:i_up
%          for j=1:j_up
%              if j == j_up
%                  fprintf(fid,'%.4f\n',double(U(i,j,U_t))); %0516-keep 4 decimals
%              else
%                  fprintf(fid,'%.4f\t',double(U(i,j,U_t))); %0516-keep 4 decimals
%              end     
%          end    
%      end    
%      fclose(fid);
% end

% 0607-calculate Upper Pijt upper bound
for i=1:i_up
  for j=1:j_up
    for t=1:t_up
      P(i,j,t)=(b*d(i,j,t)+cijt*d(i,j,t)*m-n*T_abs*U(i,j,t)-(a+b)*min_k*d(i,j,t))./(m*min_k*d(i,j,t))-cijt;
    end
  end
end

% fprintf('\nUpper bound P(1,1,1) = %.4f\n',double(P(1,1,1)));

%0608-----
%calculate Upper revenue
%0608-----

revenue=0;

for i=1:i_up
  for j=1:j_up
    for t=1:t_up
      if (t==1)||(t==2)||(t==9)||(t==10)||(t==17)||(t==18)||(t==31)||(t==32)||(t==39)||(t==40)
%         T_abs=2;
        revenue = revenue-(P(i,j,t).*(b.*d(i,j,t)+cijt.*d(i,j,t)*m-n*T_abs*U(i,j,t)))./(a+b+m.*(cijt+P(i,j,t))); 
      else
%         T_abs=2; %here should be 3, but temporarily set to 2, aligh with Oper level
        revenue = revenue-(P(i,j,t).*(b.*d(i,j,t)+cijt.*d(i,j,t)*m-n*T_abs*U(i,j,t)))./(a+b+m.*(cijt+P(i,j,t)));
      end
    end
  end
end

% 070320-avg. all P

 all_P = [];

  for i=1:i_up
    for j=1:j_up
      for t=1:t_up
        all_P = P(i,j,t);
      end
    end
  end

avg_P = all_P./(i_up*j_up*t_up);



% 070320-avg. all D

 all_Dijt = [];

  for i=1:i_up
    for j=1:j_up
      for t=1:t_up
        all_Dijt = D_ij_SAV(i,j,t);
      end
    end
  end

avg_Dijt = all_Dijt./(i_up*j_up*t_up);

% 070320-print results

fprintf('\n Upper bound P(1,1,1) = %.4f\n',double(P(1,1,1))); % particular P each iteration

fprintf('\n');

fprintf('\n D(1,1,1) = %.4f\n',double(D_ij_SAV(1,1,1))); % particular Dijt each iteration

fprintf('\n');

fprintf('\n Avg_P = %.4f\n',double(avg_P)); % particular avg_P each iteration

fprintf('\n');

fprintf('\n Avg_Dijt = %.4f\n',double(avg_Dijt)); % particular avg_Dijt each iteration

fprintf('\n');

fprintf('\n Max Revenue by Upper= %.4f\n',double(-revenue));

fprintf('\n');

fprintf('\n Min Cost by Lower= %.4f\n',double(z_op));

fprintf('\n');

fprintf('\n Profit = %.4f\n',double(-revenue-z_op));

fprintf('\n');

% test how many new P > old P

satis_P = 0;

for i=1:i_up
  for j=1:j_up
    for t=1:t_up
%       if (double(P_sav(i,j,t))<=double(P(i,j,t))) % 0610-check new P greater than old ones
        %0610-if new P no more than 1 cent than old, than add to satis_P terminate condition
      %if ((double(P_sav(i,j,t)))*1.0001 >= double(P(i,j,t))) 
%       if ((double(P_sav(i,j,t)))+0.01 >= double(P(i,j,t))) 
      if ((P_sav(i,j,t))+0.01 > P(i,j,t)) 
        satis_P = (satis_P)+1;
      end  
    end
  end
end

satis_P;

fprintf('\n');

% test how many new P all feasible that new_P <= new P upper bound condition


if (Iter_count==1)
  feasib_P = 0;
else
  
  feasib_P = 0;
  
  for i=1:i_up
    for j=1:j_up
      for t=1:t_up
        P_bench(i,j,t) = (b*d(i,j,t)+cijt*d(i,j,t)*m-n*T_abs*U(i,j,t)-(a+b)*min_k*d(i,j,t))./(m*min_k*d(i,j,t))-cijt;        
      end
    end
  end
  
    
  for i=1:i_up
    for j=1:j_up
      for t=1:t_up
  %       if (double(P_sav(i,j,t))<=double(P(i,j,t))) % 0610-check new P greater than old ones
          %0610-if new P no more than 1 cent than old, than add to satis_P terminate condition
        %if ((double(P_sav(i,j,t)))*1.0001 >= double(P(i,j,t))) 
        
%         if (double(P_sav(i,j,t)) <= double(P_bench(i,j,t)))
        if (P_bench(i,j,t) >= P_sav(i,j,t))
%           i,j,t
%           fprintf("\n");
          feasib_P = (feasib_P)+1;
        end  
      end
    end
  end

end

feasib_P;

fprintf('\n');

if (satis_P==i_up*j_up*t_up) && (feasib_P==i_up*j_up*t_up)
  consec_P = 1+consec_P;
else
  consec_P = 0;
end

% tik = tik-1;

end

fprintf('\n');
fprintf('\n');
fprintf('\n');


%-------case 3.3------

'-------Case 3.3-----'

IniP=100;

%initially price Pijt all the same
P=zeros(20,20,40);%create 20*20*40 matrix
for i=1:i_up 
    for j=1:j_up
        for t=1:t_up
          P(i,j,t)=IniP; %0608-initial P should be infinity, I use 2 for test
        end
    end    
end

% Pijt=2;

Iter_count = 0;

% 0517-M is c16 time window, set at 1 gives it relaxation
% 0424-M is c16 time window, set at 3 gives it relaxation
M=1; 

%0515-change all var to continours instead of MIP
X = sdpvar(i_up,j_up,t_up+M); %0424-Xi,j,t+M, cuz time window constraint c16 Xijt is beyond t
%X = intvar(i_up,j_up,t_up+M); %0424-Xi,j,t+M, cuz time window constraint c16 Xijt is beyond t
X_t0 = sdpvar(i_up,j_up,1); %0424-Xi,j,0 means moves initially for c18: inventory balance eq
%X_t0 = intvar(i_up,j_up,1); %0424-Xi,j,0 means moves initially for c18: inventory balance eq
U = sdpvar(i_up,j_up,t_up);
%U = intvar(i_up,j_up,t_up); 
U_t0 = sdpvar(i_up,j_up,1); %0424-this is for Ui,j,0 
%U_t0 = intvar(i_up,j_up,1); %0424-this is for Ui,j,0 
V = sdpvar(i_up,t_up); 
%V = intvar(i_up,t_up); 
V_t0 = sdpvar(i_up,1); %0424-this is especially for Vi,0
%V_t0 = intvar(i_up,1); %0424-this is especially for Vi,0
D_ij_SAV = sdpvar(i_up,j_up,t_up); %0425-this is to show D_ij_SAV, for futher plot, then find upper bound of Pijt
%D_ij_SAV = intvar(i_up,j_up,t_up); %0425-this is to show D_ij_SAV, for futher plot, then find upper bound of Pijt
% 0604-k ratio
k = sdpvar(i_up,j_up,t_up);

% 070220-set place for Wijt
W=zeros(20,20,40);%create 20*20*40 matrix

%1.1-parameters
% a=b=m=n=1,T=40
a=1;
b=1;
m=10;
n=1;
T=8; %070120-total 40 time periods, 5 time block, so every block is 8
T_abs=8; %070120-T_abs should be same as T
% cijt=40, pijt=10, hit=10, alphaijt=0.1 (all c,p,h,alpha are subjectively to be same)
cijt=20;
pijt=10; % 070120-increase for penalty of unmet demand
% 0425-intentionally make pijt this small to avoid all 0 in Uijt solution
hit=10;
alpha=1;

%%%-----0607----
%%% How to add loop?
%%%----

% 0610-here starts the loop

P_sav = zeros(20,20,40);
U_sav = zeros(20,20,40);
min_k = 0;
P_bench = zeros(20,20,40);

% tik = 2;

satis_P = 0;
feasib_P = 0;
consec_P = 0;

% while (satis_P<i_up*j_up*t_up) || (feasib_P<i_up*j_up*t_up)
% while Iter_count < 2
while (consec_P<3)

Iter_count = Iter_count+1; 


for i=1:i_up 
    for j=1:j_up
        for t=1:t_up
          P_sav(i,j,t)=P(i,j,t); %0608-initial P should be infinity, I use 2 for test
        end
    end    
end  


for i=1:i_up 
    for j=1:j_up
        for t=1:t_up
            U_sav(i,j,t)=U(i,j,t); 
        end
    end    
end


min_k_sav = min_k;


%1.2-add operational level constraints

C = [];

% c15
for i=1:i_up %c15
   for j=1:j_up
     for t=1:t_up
       if (t==1)
         C = [C, U(i,j,1) >= U_t0(i,j,1)+(b*d(i,j,1)+cijt*d(i,j,1)*m-n*T_abs*U(i,j,1))./(a+b+m*(cijt+P(i,j,t)))-X(i,j,1)]; %c15: t=1
       else       
         C = [C, U(i,j,t) >= U(i,j,t-1)+(b*d(i,j,t)+cijt*d(i,j,t)*m-n*T_abs*U(i,j,t))./(a+b+m*(cijt+P(i,j,t)))-X(i,j,t)]; %c15: t from 2
       end
     end
   end    
end

% 0424-c16-I suppose time window is M period for c16 - relax this one, hope get solution
% this also considers Xijt beyond T
for i=1:i_up %c16
    for j=1:j_up
      for t=1:t_up
        if (t==1)
          C = [C, sum(X(i,j,t:(t+M))) >= U_t0(i,j,1)+(b*d(i,j,t)+cijt*d(i,j,t)*m-n*T_abs*U(i,j,t))./(a+b+m*(cijt+P(i,j,t)))]; %c16: t=1
        else
          C = [C, sum(X(i,j,t:(t+M))) >= U(i,j,t-1)+(b*d(i,j,t)+cijt*d(i,j,t)*m-n*T_abs*U(i,j,t))./(a+b+m*(cijt+P(i,j,t)))]; %c16: t from 2
        end
      end  
    end     
end

% 0424 - c17 is in upper level, no need in lower level
% for i=1:i_up %c17
%   for j=1:j_up
%     for t=1:t_up
%       C = [C, ((b-1)*d(i,j,t)+cijt*d(i,j,t)*m)./(n*T_abs) <= U(i,j,t)];
%       C = [C, U(i,j,t) <= (b*d(i,j,t)+cijt*d(i,j,t)*m)./(n*T_abs)];
%     end
%   end
% end

% 0424-c18-here makes <= for Vi,t cuz alpha all set to 1, too relaxed, should be variaty
for i=1:i_up %c18
  for t=1:t_up
    if (t==1)
      C = [C, V(i,t) <= V_t0(i,1) + sum(alpha*X_t0(1:j_up,i,1)) - sum(X(i,1:j_up,t)) ]; %c18: t=1
    else
      C = [C, V(i,t) <= V(i,t-1) + sum(sum(alpha*X(1:j_up,i,1:(t-1)))) - sum(X(i,1:j_up,t)) ]; %c18: t from 2 (remember to use 2 sum for 2nd)
    end
  end
end

% 0508-c19-not added, cuz it looks already relaxed in math function, hidden property of our formula

for i=1:i_up %c20: Xijt
    for j=1:j_up
        for t=1:(t_up+M)
          C = [C, X(i,j,t) >= 0];
        end
    end    
end    

for i=1:i_up %c20: X_t0(i,j,1) is Xijt=0
    for j=1:j_up
        C = [C, X_t0(i,j,1) >= 0];
    end    
end    

for i=1:i_up %c20: Uijt
    for j=1:j_up
        for t=1:t_up
            C = [C, U(i,j,t) >= 0];
        end
    end    
end    

for i=1:i_up %c20: U_t0(i,j,1) is Uijt=0
    for j=1:j_up
        C = [C, U_t0(i,j,1) >= 0];
    end    
end 

for i=1:i_up %c20: Vit
    for t=1:t_up
      C = [C, V(i,t) >= 0];
    end    
end    

for i=1:i_up %c20: V_t0(i,1) is V(i,0)
  C = [C, V_t0(i,1) >= 0];
end

%1.3-solve operational level: settings: use cplex as solver

%0515-following is set optimality for MIP, but kang said make it LNContinuous-https://www.ibm.com/support/knowledgecenter/SSSA5P_12.9.0/ilog.odms.cplex.help/CPLEX/Parameters/topics/OptimalityTarget.html
%0515-'cplex.optimalitytarget'=1: Searches for a globally optimal solution to a convex model.
ops = sdpsettings ('solver','cplex','verbose',2,'cplex.optimalitytarget',1);%0508-verbose 2 shows iterations, 1 means shown normal message, 0 shows silent message
%0515-MIP limit tolerance gap to 1-https://www.ibm.com/support/knowledgecenter/SSSA5P_12.9.0/ilog.odms.cplex.help/CPLEX/Parameters/topics/EpGap.html
%ops = sdpsettings ('solver','cplex','verbose',2,'cplex.mip.tolerances.mipgap',1);%0508-verbose 2 shows iterations, 1 means shown normal message, 0 shows silent message

%0306-limit barrier iteration-https://www.ibm.com/support/knowledgecenter/SSSA5P_12.9.0/ilog.odms.cplex.help/CPLEX/Parameters/topics/BarItLim.html
%Cplex.Param.barrier.limits.iteration = 4;[not work]

%0306-limit network simplex iteration-https://www.ibm.com/support/knowledgecenter/SSSA5P_12.9.0/ilog.odms.cplex.help/CPLEX/Parameters/topics/NetItLim.html
%Cplex.Param.network.iterations = 2;[not work]

%ops.mip.limits.nodes = 4;[not work]

%0304-use the online method for out-of-memory: https://www.ilovematlab.cn/thread-266244-1-1.html
%ops.mip.strategy.file = 3; %[not work: should put in sdpsetting] this has same results 1159110 node, but reduces time for 3 min.

%0305-use online method for out-of-memory: https://www.ibm.com/developerworks/community/forums/html/topic?id=e3ea344c-7249-4edd-a47d-29e1f80fa480
%Cplex.Param.mip.limits.auxrootthreads = 2; %[not work: should put in sdpsetting] -1 or 1

%obj-operational obj z_op
z_op = sum(sum(sum(cijt.*X(1:i_up,1:j_up,1:t_up)))) + sum(sum(cijt.*X_t0(1:i_up,1:j_up,1)))...
  + sum(sum(sum(pijt.*U(1:i_up,1:j_up,1:t_up)))) + sum(sum(pijt.*U_t0(1:i_up,1:j_up,1)))...
  + sum(sum(hit.*V(1:i_up,1:t_up))) + sum(hit.*V_t0(1:i_up,1));

%1.35-express D_ij_SAV
for i=1:i_up
  for j=1:j_up
    for t=1:t_up
      D_ij_SAV(i,j,t) = (b*d(i,j,t)+cijt*d(i,j,t)*m-n*T_abs*U(i,j,t))./(a+b+m*(cijt+P(i,j,t)));
    end
  end
end

% 0604-calculate minimum k ratio covered by SAVs
for i=1:i_up
  for j=1:j_up
    for t=1:t_up
      k(i,j,t) = D_ij_SAV(i,j,t)./d(i,j,t);
    end
  end
end

% k_min = min(k(1,1,1),k(1,2,1),k(1,3,1))

% k_min = 0.99;
% % 
% if ((k_min) >= (k(1,1,1)))
%   k_min = k(1,1,1)                
% end  

% 0604-find minimum covered ratio k
% for i=1:i_up
%   for j=1:j_up
%     for t=1:t_up
%       if (k_min >= k(i,j,t))
%         k_min = k(i,j,t)                
%       end  
%     end
%   end
% end

%1.4-check if solved and output
result = optimize(C,z_op,ops); % 0508-origin
% result = solvesdp(C,z_op,ops); % 0508-this works

%result = optimize(C,z);
% if result.problem == 0 % problem =0 means solve successfully,only print Uijt here

%   fprintf('\nk value: ') %ratio k covered by SAVs
%   value(k)

  
%   fprintf('\nP_Ini_(1,1,1) = ')
%   value(P(1,1,1))

%   fprintf('\nMin cost by Lower = ')
%   value(z_op)

  all_k = [];

  for i=1:i_up
    for j=1:j_up
      for t=1:t_up
        all_k(end+1) = k(i,j,t);
      end
    end
  end
  
  min_k = min(all_k);
  
  k_avg = mean(all_k);
  
%   fprintf('\n min k: ')
%   value(min_k)
  
  max_k = max(all_k);
  
%   fprintf('\n max k: ')
%   value(max_k)
  
  % 0615-save max_k in txt

%   fid = fopen(['IniP=',num2str(2),'max_k.txt'],'a'); %save to last 'a'
%   fprintf(fid,'%.4f\t',double(max_k)); %0516-keep 4 decimals
%   fclose(fid);
  
%   fprintf('\nU(1,1,1) value: ')
%   value(U(1,1,1))
        
%   fprintf('\nX(1,1,1) value: ') 
%   value(X(1,1,1))

fprintf('\n Iteration No. = ')
value(Iter_count)



% fprintf('\n D_ij_SAV(1,1,1) value: ') %0515-as for now hide D_ij_SAV
% value(D_ij_SAV(1,1,1))
  
%   fprintf('\n k(1,1,1) ')
%   value(k(1,1,1))

%   fprintf('\n')
%     fprintf('\nV value: ')
%     value(V)
%     
%     fprintf('\nU_t0 value: ')
%     value(U_t0)
%     
%     fprintf('\nX_t0 value: ') 
%     value(X_t0)
%     
%     fprintf('\nV_t0 value: ')
%     value(V_t0)

% else
%   disp('Not solved!')
% 
% end

% 0519-save total cost z_op in txt
% fid = fopen(['Iter=',num2str(Iter_count),'_OperCost.txt'],'wt'); %creat & write z_op to txt
% fprintf(fid,'%.4f\n',double(z_op)); %0516-keep 4 decimals
% fclose(fid);

% 0605-save min ratio k covered by SAVs in txt
% fid_ratiok = fopen(['Iter=',num2str(Iter_count),'_MinRatiok.txt'],'wt'); %creat & write z_op to txt
% fprintf(fid_ratiok,'%.4f\n',double(min_k)); %0604-keep 4 decimals
% fclose(fid_ratiok);

% 0517-following is save Uijt in txt
%1.6-save Uijt as txt for upper use
% for U_t=1:t_up 
%      fid = fopen(['IniP=',num2str(IniP),'Iter=',num2str(Iter_count),'_Uijt=',num2str(U_t),'.txt'],'wt'); %creat & write Uijt to txt
%      for i=1:i_up
%          for j=1:j_up
%              if j == j_up
%                  fprintf(fid,'%.4f\n',double(U(i,j,U_t))); %0516-keep 4 decimals
%              else
%                  fprintf(fid,'%.4f\t',double(U(i,j,U_t))); %0516-keep 4 decimals
%              end     
%          end    
%      end    
%      fclose(fid);
% end

% 0607-calculate Upper Pijt upper bound
for i=1:i_up
  for j=1:j_up
    for t=1:t_up
      P(i,j,t)=(b*d(i,j,t)+cijt*d(i,j,t)*m-n*T_abs*U(i,j,t)-(a+b)*min_k*d(i,j,t))./(m*min_k*d(i,j,t))-cijt;
    end
  end
end

% fprintf('\nUpper bound P(1,1,1) = %.4f\n',double(P(1,1,1)));

%0608-----
%calculate Upper revenue
%0608-----

revenue=0;

for i=1:i_up
  for j=1:j_up
    for t=1:t_up
      if (t==1)||(t==2)||(t==9)||(t==10)||(t==17)||(t==18)||(t==31)||(t==32)||(t==39)||(t==40)
%         T_abs=2;
        revenue = revenue-(P(i,j,t).*(b.*d(i,j,t)+cijt.*d(i,j,t)*m-n*T_abs*U(i,j,t)))./(a+b+m.*(cijt+P(i,j,t))); 
      else
%         T_abs=2; %here should be 3, but temporarily set to 2, aligh with Oper level
        revenue = revenue-(P(i,j,t).*(b.*d(i,j,t)+cijt.*d(i,j,t)*m-n*T_abs*U(i,j,t)))./(a+b+m.*(cijt+P(i,j,t)));
      end
    end
  end
end

% 070320-avg. all P

 all_P = [];

  for i=1:i_up
    for j=1:j_up
      for t=1:t_up
        all_P = P(i,j,t);
      end
    end
  end

avg_P = all_P./(i_up*j_up*t_up);



% 070320-avg. all D

 all_Dijt = [];

  for i=1:i_up
    for j=1:j_up
      for t=1:t_up
        all_Dijt = D_ij_SAV(i,j,t);
      end
    end
  end

avg_Dijt = all_Dijt./(i_up*j_up*t_up);

% 070320-print results

fprintf('\n Upper bound P(1,1,1) = %.4f\n',double(P(1,1,1))); % particular P each iteration

fprintf('\n');

fprintf('\n D(1,1,1) = %.4f\n',double(D_ij_SAV(1,1,1))); % particular Dijt each iteration

fprintf('\n');

fprintf('\n Avg_P = %.4f\n',double(avg_P)); % particular avg_P each iteration

fprintf('\n');

fprintf('\n Avg_Dijt = %.4f\n',double(avg_Dijt)); % particular avg_Dijt each iteration

fprintf('\n');

fprintf('\n Max Revenue by Upper= %.4f\n',double(-revenue));

fprintf('\n');

fprintf('\n Min Cost by Lower= %.4f\n',double(z_op));

fprintf('\n');

fprintf('\n Profit = %.4f\n',double(-revenue-z_op));

fprintf('\n');

% test how many new P > old P

satis_P = 0;

for i=1:i_up
  for j=1:j_up
    for t=1:t_up
%       if (double(P_sav(i,j,t))<=double(P(i,j,t))) % 0610-check new P greater than old ones
        %0610-if new P no more than 1 cent than old, than add to satis_P terminate condition
      %if ((double(P_sav(i,j,t)))*1.0001 >= double(P(i,j,t))) 
%       if ((double(P_sav(i,j,t)))+0.01 >= double(P(i,j,t))) 
      if ((P_sav(i,j,t))+0.01 > P(i,j,t)) 
        satis_P = (satis_P)+1;
      end  
    end
  end
end

satis_P;

fprintf('\n');

% test how many new P all feasible that new_P <= new P upper bound condition


if (Iter_count==1)
  feasib_P = 0;
else
  
  feasib_P = 0;
  
  for i=1:i_up
    for j=1:j_up
      for t=1:t_up
        P_bench(i,j,t) = (b*d(i,j,t)+cijt*d(i,j,t)*m-n*T_abs*U(i,j,t)-(a+b)*min_k*d(i,j,t))./(m*min_k*d(i,j,t))-cijt;        
      end
    end
  end
  
    
  for i=1:i_up
    for j=1:j_up
      for t=1:t_up
  %       if (double(P_sav(i,j,t))<=double(P(i,j,t))) % 0610-check new P greater than old ones
          %0610-if new P no more than 1 cent than old, than add to satis_P terminate condition
        %if ((double(P_sav(i,j,t)))*1.0001 >= double(P(i,j,t))) 
        
%         if (double(P_sav(i,j,t)) <= double(P_bench(i,j,t)))
        if (P_bench(i,j,t) >= P_sav(i,j,t))
%           i,j,t
%           fprintf("\n");
          feasib_P = (feasib_P)+1;
        end  
      end
    end
  end

end

feasib_P;

fprintf('\n');

if (satis_P==i_up*j_up*t_up) && (feasib_P==i_up*j_up*t_up)
  consec_P = 1+consec_P;
else
  consec_P = 0;
end

% tik = tik-1;

end

fprintf('\n');
fprintf('\n');
fprintf('\n');



%-------case 4.1------

'------Case 4.1------'

IniP=1;

%initially price Pijt all the same
P=zeros(20,20,40);%create 20*20*40 matrix
for i=1:i_up 
    for j=1:j_up
        for t=1:t_up
          P(i,j,t)=IniP; %0608-initial P should be infinity, I use 2 for test
        end
    end    
end

% Pijt=2;

Iter_count = 0;

% 0517-M is c16 time window, set at 1 gives it relaxation
% 0424-M is c16 time window, set at 3 gives it relaxation
M=1; 

%0515-change all var to continours instead of MIP
X = sdpvar(i_up,j_up,t_up+M); %0424-Xi,j,t+M, cuz time window constraint c16 Xijt is beyond t
%X = intvar(i_up,j_up,t_up+M); %0424-Xi,j,t+M, cuz time window constraint c16 Xijt is beyond t
X_t0 = sdpvar(i_up,j_up,1); %0424-Xi,j,0 means moves initially for c18: inventory balance eq
%X_t0 = intvar(i_up,j_up,1); %0424-Xi,j,0 means moves initially for c18: inventory balance eq
U = sdpvar(i_up,j_up,t_up);
%U = intvar(i_up,j_up,t_up); 
U_t0 = sdpvar(i_up,j_up,1); %0424-this is for Ui,j,0 
%U_t0 = intvar(i_up,j_up,1); %0424-this is for Ui,j,0 
V = sdpvar(i_up,t_up); 
%V = intvar(i_up,t_up); 
V_t0 = sdpvar(i_up,1); %0424-this is especially for Vi,0
%V_t0 = intvar(i_up,1); %0424-this is especially for Vi,0
D_ij_SAV = sdpvar(i_up,j_up,t_up); %0425-this is to show D_ij_SAV, for futher plot, then find upper bound of Pijt
%D_ij_SAV = intvar(i_up,j_up,t_up); %0425-this is to show D_ij_SAV, for futher plot, then find upper bound of Pijt
% 0604-k ratio
k = sdpvar(i_up,j_up,t_up);

% 070220-set place for Wijt
W=zeros(20,20,40);%create 20*20*40 matrix

%1.1-parameters
% a=b=m=n=1,T=40
a=1;
b=1;
m=1;
n=10;
T=8; %070120-total 40 time periods, 5 time block, so every block is 8
T_abs=8; %070120-T_abs should be same as T
% cijt=40, pijt=10, hit=10, alphaijt=0.1 (all c,p,h,alpha are subjectively to be same)
cijt=20;
pijt=10; % 070120-increase for penalty of unmet demand
% 0425-intentionally make pijt this small to avoid all 0 in Uijt solution
hit=10;
alpha=1;

%%%-----0607----
%%% How to add loop?
%%%----

% 0610-here starts the loop

P_sav = zeros(20,20,40);
U_sav = zeros(20,20,40);
min_k = 0;
P_bench = zeros(20,20,40);

% tik = 2;

satis_P = 0;
feasib_P = 0;
consec_P = 0;

% while (satis_P<i_up*j_up*t_up) || (feasib_P<i_up*j_up*t_up)
% while Iter_count < 2
while (consec_P<3)

Iter_count = Iter_count+1; 


for i=1:i_up 
    for j=1:j_up
        for t=1:t_up
          P_sav(i,j,t)=P(i,j,t); %0608-initial P should be infinity, I use 2 for test
        end
    end    
end  


for i=1:i_up 
    for j=1:j_up
        for t=1:t_up
            U_sav(i,j,t)=U(i,j,t); 
        end
    end    
end


min_k_sav = min_k;


%1.2-add operational level constraints

C = [];

% c15
for i=1:i_up %c15
   for j=1:j_up
     for t=1:t_up
       if (t==1)
         C = [C, U(i,j,1) >= U_t0(i,j,1)+(b*d(i,j,1)+cijt*d(i,j,1)*m-n*T_abs*U(i,j,1))./(a+b+m*(cijt+P(i,j,t)))-X(i,j,1)]; %c15: t=1
       else       
         C = [C, U(i,j,t) >= U(i,j,t-1)+(b*d(i,j,t)+cijt*d(i,j,t)*m-n*T_abs*U(i,j,t))./(a+b+m*(cijt+P(i,j,t)))-X(i,j,t)]; %c15: t from 2
       end
     end
   end    
end

% 0424-c16-I suppose time window is M period for c16 - relax this one, hope get solution
% this also considers Xijt beyond T
for i=1:i_up %c16
    for j=1:j_up
      for t=1:t_up
        if (t==1)
          C = [C, sum(X(i,j,t:(t+M))) >= U_t0(i,j,1)+(b*d(i,j,t)+cijt*d(i,j,t)*m-n*T_abs*U(i,j,t))./(a+b+m*(cijt+P(i,j,t)))]; %c16: t=1
        else
          C = [C, sum(X(i,j,t:(t+M))) >= U(i,j,t-1)+(b*d(i,j,t)+cijt*d(i,j,t)*m-n*T_abs*U(i,j,t))./(a+b+m*(cijt+P(i,j,t)))]; %c16: t from 2
        end
      end  
    end     
end

% 0424 - c17 is in upper level, no need in lower level
% for i=1:i_up %c17
%   for j=1:j_up
%     for t=1:t_up
%       C = [C, ((b-1)*d(i,j,t)+cijt*d(i,j,t)*m)./(n*T_abs) <= U(i,j,t)];
%       C = [C, U(i,j,t) <= (b*d(i,j,t)+cijt*d(i,j,t)*m)./(n*T_abs)];
%     end
%   end
% end

% 0424-c18-here makes <= for Vi,t cuz alpha all set to 1, too relaxed, should be variaty
for i=1:i_up %c18
  for t=1:t_up
    if (t==1)
      C = [C, V(i,t) <= V_t0(i,1) + sum(alpha*X_t0(1:j_up,i,1)) - sum(X(i,1:j_up,t)) ]; %c18: t=1
    else
      C = [C, V(i,t) <= V(i,t-1) + sum(sum(alpha*X(1:j_up,i,1:(t-1)))) - sum(X(i,1:j_up,t)) ]; %c18: t from 2 (remember to use 2 sum for 2nd)
    end
  end
end

% 0508-c19-not added, cuz it looks already relaxed in math function, hidden property of our formula

for i=1:i_up %c20: Xijt
    for j=1:j_up
        for t=1:(t_up+M)
          C = [C, X(i,j,t) >= 0];
        end
    end    
end    

for i=1:i_up %c20: X_t0(i,j,1) is Xijt=0
    for j=1:j_up
        C = [C, X_t0(i,j,1) >= 0];
    end    
end    

for i=1:i_up %c20: Uijt
    for j=1:j_up
        for t=1:t_up
            C = [C, U(i,j,t) >= 0];
        end
    end    
end    

for i=1:i_up %c20: U_t0(i,j,1) is Uijt=0
    for j=1:j_up
        C = [C, U_t0(i,j,1) >= 0];
    end    
end 

for i=1:i_up %c20: Vit
    for t=1:t_up
      C = [C, V(i,t) >= 0];
    end    
end    

for i=1:i_up %c20: V_t0(i,1) is V(i,0)
  C = [C, V_t0(i,1) >= 0];
end

%1.3-solve operational level: settings: use cplex as solver

%0515-following is set optimality for MIP, but kang said make it LNContinuous-https://www.ibm.com/support/knowledgecenter/SSSA5P_12.9.0/ilog.odms.cplex.help/CPLEX/Parameters/topics/OptimalityTarget.html
%0515-'cplex.optimalitytarget'=1: Searches for a globally optimal solution to a convex model.
ops = sdpsettings ('solver','cplex','verbose',2,'cplex.optimalitytarget',1);%0508-verbose 2 shows iterations, 1 means shown normal message, 0 shows silent message
%0515-MIP limit tolerance gap to 1-https://www.ibm.com/support/knowledgecenter/SSSA5P_12.9.0/ilog.odms.cplex.help/CPLEX/Parameters/topics/EpGap.html
%ops = sdpsettings ('solver','cplex','verbose',2,'cplex.mip.tolerances.mipgap',1);%0508-verbose 2 shows iterations, 1 means shown normal message, 0 shows silent message

%0306-limit barrier iteration-https://www.ibm.com/support/knowledgecenter/SSSA5P_12.9.0/ilog.odms.cplex.help/CPLEX/Parameters/topics/BarItLim.html
%Cplex.Param.barrier.limits.iteration = 4;[not work]

%0306-limit network simplex iteration-https://www.ibm.com/support/knowledgecenter/SSSA5P_12.9.0/ilog.odms.cplex.help/CPLEX/Parameters/topics/NetItLim.html
%Cplex.Param.network.iterations = 2;[not work]

%ops.mip.limits.nodes = 4;[not work]

%0304-use the online method for out-of-memory: https://www.ilovematlab.cn/thread-266244-1-1.html
%ops.mip.strategy.file = 3; %[not work: should put in sdpsetting] this has same results 1159110 node, but reduces time for 3 min.

%0305-use online method for out-of-memory: https://www.ibm.com/developerworks/community/forums/html/topic?id=e3ea344c-7249-4edd-a47d-29e1f80fa480
%Cplex.Param.mip.limits.auxrootthreads = 2; %[not work: should put in sdpsetting] -1 or 1

%obj-operational obj z_op
z_op = sum(sum(sum(cijt.*X(1:i_up,1:j_up,1:t_up)))) + sum(sum(cijt.*X_t0(1:i_up,1:j_up,1)))...
  + sum(sum(sum(pijt.*U(1:i_up,1:j_up,1:t_up)))) + sum(sum(pijt.*U_t0(1:i_up,1:j_up,1)))...
  + sum(sum(hit.*V(1:i_up,1:t_up))) + sum(hit.*V_t0(1:i_up,1));

%1.35-express D_ij_SAV
for i=1:i_up
  for j=1:j_up
    for t=1:t_up
      D_ij_SAV(i,j,t) = (b*d(i,j,t)+cijt*d(i,j,t)*m-n*T_abs*U(i,j,t))./(a+b+m*(cijt+P(i,j,t)));
    end
  end
end

% 0604-calculate minimum k ratio covered by SAVs
for i=1:i_up
  for j=1:j_up
    for t=1:t_up
      k(i,j,t) = D_ij_SAV(i,j,t)./d(i,j,t);
    end
  end
end

% k_min = min(k(1,1,1),k(1,2,1),k(1,3,1))

% k_min = 0.99;
% % 
% if ((k_min) >= (k(1,1,1)))
%   k_min = k(1,1,1)                
% end  

% 0604-find minimum covered ratio k
% for i=1:i_up
%   for j=1:j_up
%     for t=1:t_up
%       if (k_min >= k(i,j,t))
%         k_min = k(i,j,t)                
%       end  
%     end
%   end
% end

%1.4-check if solved and output
result = optimize(C,z_op,ops); % 0508-origin
% result = solvesdp(C,z_op,ops); % 0508-this works

%result = optimize(C,z);
% if result.problem == 0 % problem =0 means solve successfully,only print Uijt here

%   fprintf('\nk value: ') %ratio k covered by SAVs
%   value(k)

  
%   fprintf('\nP_Ini_(1,1,1) = ')
%   value(P(1,1,1))

%   fprintf('\nMin cost by Lower = ')
%   value(z_op)

  all_k = [];

  for i=1:i_up
    for j=1:j_up
      for t=1:t_up
        all_k(end+1) = k(i,j,t);
      end
    end
  end
  
  min_k = min(all_k);
  
  k_avg = mean(all_k);
  
%   fprintf('\n min k: ')
%   value(min_k)
  
  max_k = max(all_k);
  
%   fprintf('\n max k: ')
%   value(max_k)
  
  % 0615-save max_k in txt

%   fid = fopen(['IniP=',num2str(2),'max_k.txt'],'a'); %save to last 'a'
%   fprintf(fid,'%.4f\t',double(max_k)); %0516-keep 4 decimals
%   fclose(fid);
  
%   fprintf('\nU(1,1,1) value: ')
%   value(U(1,1,1))
        
%   fprintf('\nX(1,1,1) value: ') 
%   value(X(1,1,1))

fprintf('\n Iteration No. = ')
value(Iter_count)



% fprintf('\n D_ij_SAV(1,1,1) value: ') %0515-as for now hide D_ij_SAV
% value(D_ij_SAV(1,1,1))
  
%   fprintf('\n k(1,1,1) ')
%   value(k(1,1,1))

%   fprintf('\n')
%     fprintf('\nV value: ')
%     value(V)
%     
%     fprintf('\nU_t0 value: ')
%     value(U_t0)
%     
%     fprintf('\nX_t0 value: ') 
%     value(X_t0)
%     
%     fprintf('\nV_t0 value: ')
%     value(V_t0)

% else
%   disp('Not solved!')
% 
% end

% 0519-save total cost z_op in txt
% fid = fopen(['Iter=',num2str(Iter_count),'_OperCost.txt'],'wt'); %creat & write z_op to txt
% fprintf(fid,'%.4f\n',double(z_op)); %0516-keep 4 decimals
% fclose(fid);

% 0605-save min ratio k covered by SAVs in txt
% fid_ratiok = fopen(['Iter=',num2str(Iter_count),'_MinRatiok.txt'],'wt'); %creat & write z_op to txt
% fprintf(fid_ratiok,'%.4f\n',double(min_k)); %0604-keep 4 decimals
% fclose(fid_ratiok);

% 0517-following is save Uijt in txt
%1.6-save Uijt as txt for upper use
% for U_t=1:t_up 
%      fid = fopen(['IniP=',num2str(IniP),'Iter=',num2str(Iter_count),'_Uijt=',num2str(U_t),'.txt'],'wt'); %creat & write Uijt to txt
%      for i=1:i_up
%          for j=1:j_up
%              if j == j_up
%                  fprintf(fid,'%.4f\n',double(U(i,j,U_t))); %0516-keep 4 decimals
%              else
%                  fprintf(fid,'%.4f\t',double(U(i,j,U_t))); %0516-keep 4 decimals
%              end     
%          end    
%      end    
%      fclose(fid);
% end

% 0607-calculate Upper Pijt upper bound
for i=1:i_up
  for j=1:j_up
    for t=1:t_up
      P(i,j,t)=(b*d(i,j,t)+cijt*d(i,j,t)*m-n*T_abs*U(i,j,t)-(a+b)*min_k*d(i,j,t))./(m*min_k*d(i,j,t))-cijt;
    end
  end
end

% fprintf('\nUpper bound P(1,1,1) = %.4f\n',double(P(1,1,1)));

%0608-----
%calculate Upper revenue
%0608-----

revenue=0;

for i=1:i_up
  for j=1:j_up
    for t=1:t_up
      if (t==1)||(t==2)||(t==9)||(t==10)||(t==17)||(t==18)||(t==31)||(t==32)||(t==39)||(t==40)
%         T_abs=2;
        revenue = revenue-(P(i,j,t).*(b.*d(i,j,t)+cijt.*d(i,j,t)*m-n*T_abs*U(i,j,t)))./(a+b+m.*(cijt+P(i,j,t))); 
      else
%         T_abs=2; %here should be 3, but temporarily set to 2, aligh with Oper level
        revenue = revenue-(P(i,j,t).*(b.*d(i,j,t)+cijt.*d(i,j,t)*m-n*T_abs*U(i,j,t)))./(a+b+m.*(cijt+P(i,j,t)));
      end
    end
  end
end

% 070320-avg. all P

 all_P = [];

  for i=1:i_up
    for j=1:j_up
      for t=1:t_up
        all_P = P(i,j,t);
      end
    end
  end

avg_P = all_P./(i_up*j_up*t_up);



% 070320-avg. all D

 all_Dijt = [];

  for i=1:i_up
    for j=1:j_up
      for t=1:t_up
        all_Dijt = D_ij_SAV(i,j,t);
      end
    end
  end

avg_Dijt = all_Dijt./(i_up*j_up*t_up);

% 070320-print results

fprintf('\n Upper bound P(1,1,1) = %.4f\n',double(P(1,1,1))); % particular P each iteration

fprintf('\n');

fprintf('\n D(1,1,1) = %.4f\n',double(D_ij_SAV(1,1,1))); % particular Dijt each iteration

fprintf('\n');

fprintf('\n Avg_P = %.4f\n',double(avg_P)); % particular avg_P each iteration

fprintf('\n');

fprintf('\n Avg_Dijt = %.4f\n',double(avg_Dijt)); % particular avg_Dijt each iteration

fprintf('\n');

fprintf('\n Max Revenue by Upper= %.4f\n',double(-revenue));

fprintf('\n');

fprintf('\n Min Cost by Lower= %.4f\n',double(z_op));

fprintf('\n');

fprintf('\n Profit = %.4f\n',double(-revenue-z_op));

fprintf('\n');

% test how many new P > old P

satis_P = 0;

for i=1:i_up
  for j=1:j_up
    for t=1:t_up
%       if (double(P_sav(i,j,t))<=double(P(i,j,t))) % 0610-check new P greater than old ones
        %0610-if new P no more than 1 cent than old, than add to satis_P terminate condition
      %if ((double(P_sav(i,j,t)))*1.0001 >= double(P(i,j,t))) 
%       if ((double(P_sav(i,j,t)))+0.01 >= double(P(i,j,t))) 
      if ((P_sav(i,j,t))+0.01 > P(i,j,t)) 
        satis_P = (satis_P)+1;
      end  
    end
  end
end

satis_P;

fprintf('\n');

% test how many new P all feasible that new_P <= new P upper bound condition


if (Iter_count==1)
  feasib_P = 0;
else
  
  feasib_P = 0;
  
  for i=1:i_up
    for j=1:j_up
      for t=1:t_up
        P_bench(i,j,t) = (b*d(i,j,t)+cijt*d(i,j,t)*m-n*T_abs*U(i,j,t)-(a+b)*min_k*d(i,j,t))./(m*min_k*d(i,j,t))-cijt;        
      end
    end
  end
  
    
  for i=1:i_up
    for j=1:j_up
      for t=1:t_up
  %       if (double(P_sav(i,j,t))<=double(P(i,j,t))) % 0610-check new P greater than old ones
          %0610-if new P no more than 1 cent than old, than add to satis_P terminate condition
        %if ((double(P_sav(i,j,t)))*1.0001 >= double(P(i,j,t))) 
        
%         if (double(P_sav(i,j,t)) <= double(P_bench(i,j,t)))
        if (P_bench(i,j,t) >= P_sav(i,j,t))
%           i,j,t
%           fprintf("\n");
          feasib_P = (feasib_P)+1;
        end  
      end
    end
  end

end

feasib_P;

fprintf('\n');

if (satis_P==i_up*j_up*t_up) && (feasib_P==i_up*j_up*t_up)
  consec_P = 1+consec_P;
else
  consec_P = 0;
end

% tik = tik-1;

end

fprintf('\n');
fprintf('\n');
fprintf('\n');


%-------case 4.2------

'-------Case 4.2-------'

IniP=10;

%initially price Pijt all the same
P=zeros(20,20,40);%create 20*20*40 matrix
for i=1:i_up 
    for j=1:j_up
        for t=1:t_up
          P(i,j,t)=IniP; %0608-initial P should be infinity, I use 2 for test
        end
    end    
end

% Pijt=2;

Iter_count = 0;

% 0517-M is c16 time window, set at 1 gives it relaxation
% 0424-M is c16 time window, set at 3 gives it relaxation
M=1; 

%0515-change all var to continours instead of MIP
X = sdpvar(i_up,j_up,t_up+M); %0424-Xi,j,t+M, cuz time window constraint c16 Xijt is beyond t
%X = intvar(i_up,j_up,t_up+M); %0424-Xi,j,t+M, cuz time window constraint c16 Xijt is beyond t
X_t0 = sdpvar(i_up,j_up,1); %0424-Xi,j,0 means moves initially for c18: inventory balance eq
%X_t0 = intvar(i_up,j_up,1); %0424-Xi,j,0 means moves initially for c18: inventory balance eq
U = sdpvar(i_up,j_up,t_up);
%U = intvar(i_up,j_up,t_up); 
U_t0 = sdpvar(i_up,j_up,1); %0424-this is for Ui,j,0 
%U_t0 = intvar(i_up,j_up,1); %0424-this is for Ui,j,0 
V = sdpvar(i_up,t_up); 
%V = intvar(i_up,t_up); 
V_t0 = sdpvar(i_up,1); %0424-this is especially for Vi,0
%V_t0 = intvar(i_up,1); %0424-this is especially for Vi,0
D_ij_SAV = sdpvar(i_up,j_up,t_up); %0425-this is to show D_ij_SAV, for futher plot, then find upper bound of Pijt
%D_ij_SAV = intvar(i_up,j_up,t_up); %0425-this is to show D_ij_SAV, for futher plot, then find upper bound of Pijt
% 0604-k ratio
k = sdpvar(i_up,j_up,t_up);

% 070220-set place for Wijt
W=zeros(20,20,40);%create 20*20*40 matrix

%1.1-parameters
% a=b=m=n=1,T=40
a=1;
b=1;
m=1;
n=10;
T=8; %070120-total 40 time periods, 5 time block, so every block is 8
T_abs=8; %070120-T_abs should be same as T
% cijt=40, pijt=10, hit=10, alphaijt=0.1 (all c,p,h,alpha are subjectively to be same)
cijt=20;
pijt=10; % 070120-increase for penalty of unmet demand
% 0425-intentionally make pijt this small to avoid all 0 in Uijt solution
hit=10;
alpha=1;

%%%-----0607----
%%% How to add loop?
%%%----

% 0610-here starts the loop

P_sav = zeros(20,20,40);
U_sav = zeros(20,20,40);
min_k = 0;
P_bench = zeros(20,20,40);

% tik = 2;

satis_P = 0;
feasib_P = 0;
consec_P = 0;

% while (satis_P<i_up*j_up*t_up) || (feasib_P<i_up*j_up*t_up)
% while Iter_count < 2
while (consec_P<3)

Iter_count = Iter_count+1; 


for i=1:i_up 
    for j=1:j_up
        for t=1:t_up
          P_sav(i,j,t)=P(i,j,t); %0608-initial P should be infinity, I use 2 for test
        end
    end    
end  


for i=1:i_up 
    for j=1:j_up
        for t=1:t_up
            U_sav(i,j,t)=U(i,j,t); 
        end
    end    
end


min_k_sav = min_k;


%1.2-add operational level constraints

C = [];

% c15
for i=1:i_up %c15
   for j=1:j_up
     for t=1:t_up
       if (t==1)
         C = [C, U(i,j,1) >= U_t0(i,j,1)+(b*d(i,j,1)+cijt*d(i,j,1)*m-n*T_abs*U(i,j,1))./(a+b+m*(cijt+P(i,j,t)))-X(i,j,1)]; %c15: t=1
       else       
         C = [C, U(i,j,t) >= U(i,j,t-1)+(b*d(i,j,t)+cijt*d(i,j,t)*m-n*T_abs*U(i,j,t))./(a+b+m*(cijt+P(i,j,t)))-X(i,j,t)]; %c15: t from 2
       end
     end
   end    
end

% 0424-c16-I suppose time window is M period for c16 - relax this one, hope get solution
% this also considers Xijt beyond T
for i=1:i_up %c16
    for j=1:j_up
      for t=1:t_up
        if (t==1)
          C = [C, sum(X(i,j,t:(t+M))) >= U_t0(i,j,1)+(b*d(i,j,t)+cijt*d(i,j,t)*m-n*T_abs*U(i,j,t))./(a+b+m*(cijt+P(i,j,t)))]; %c16: t=1
        else
          C = [C, sum(X(i,j,t:(t+M))) >= U(i,j,t-1)+(b*d(i,j,t)+cijt*d(i,j,t)*m-n*T_abs*U(i,j,t))./(a+b+m*(cijt+P(i,j,t)))]; %c16: t from 2
        end
      end  
    end     
end

% 0424 - c17 is in upper level, no need in lower level
% for i=1:i_up %c17
%   for j=1:j_up
%     for t=1:t_up
%       C = [C, ((b-1)*d(i,j,t)+cijt*d(i,j,t)*m)./(n*T_abs) <= U(i,j,t)];
%       C = [C, U(i,j,t) <= (b*d(i,j,t)+cijt*d(i,j,t)*m)./(n*T_abs)];
%     end
%   end
% end

% 0424-c18-here makes <= for Vi,t cuz alpha all set to 1, too relaxed, should be variaty
for i=1:i_up %c18
  for t=1:t_up
    if (t==1)
      C = [C, V(i,t) <= V_t0(i,1) + sum(alpha*X_t0(1:j_up,i,1)) - sum(X(i,1:j_up,t)) ]; %c18: t=1
    else
      C = [C, V(i,t) <= V(i,t-1) + sum(sum(alpha*X(1:j_up,i,1:(t-1)))) - sum(X(i,1:j_up,t)) ]; %c18: t from 2 (remember to use 2 sum for 2nd)
    end
  end
end

% 0508-c19-not added, cuz it looks already relaxed in math function, hidden property of our formula

for i=1:i_up %c20: Xijt
    for j=1:j_up
        for t=1:(t_up+M)
          C = [C, X(i,j,t) >= 0];
        end
    end    
end    

for i=1:i_up %c20: X_t0(i,j,1) is Xijt=0
    for j=1:j_up
        C = [C, X_t0(i,j,1) >= 0];
    end    
end    

for i=1:i_up %c20: Uijt
    for j=1:j_up
        for t=1:t_up
            C = [C, U(i,j,t) >= 0];
        end
    end    
end    

for i=1:i_up %c20: U_t0(i,j,1) is Uijt=0
    for j=1:j_up
        C = [C, U_t0(i,j,1) >= 0];
    end    
end 

for i=1:i_up %c20: Vit
    for t=1:t_up
      C = [C, V(i,t) >= 0];
    end    
end    

for i=1:i_up %c20: V_t0(i,1) is V(i,0)
  C = [C, V_t0(i,1) >= 0];
end

%1.3-solve operational level: settings: use cplex as solver

%0515-following is set optimality for MIP, but kang said make it LNContinuous-https://www.ibm.com/support/knowledgecenter/SSSA5P_12.9.0/ilog.odms.cplex.help/CPLEX/Parameters/topics/OptimalityTarget.html
%0515-'cplex.optimalitytarget'=1: Searches for a globally optimal solution to a convex model.
ops = sdpsettings ('solver','cplex','verbose',2,'cplex.optimalitytarget',1);%0508-verbose 2 shows iterations, 1 means shown normal message, 0 shows silent message
%0515-MIP limit tolerance gap to 1-https://www.ibm.com/support/knowledgecenter/SSSA5P_12.9.0/ilog.odms.cplex.help/CPLEX/Parameters/topics/EpGap.html
%ops = sdpsettings ('solver','cplex','verbose',2,'cplex.mip.tolerances.mipgap',1);%0508-verbose 2 shows iterations, 1 means shown normal message, 0 shows silent message

%0306-limit barrier iteration-https://www.ibm.com/support/knowledgecenter/SSSA5P_12.9.0/ilog.odms.cplex.help/CPLEX/Parameters/topics/BarItLim.html
%Cplex.Param.barrier.limits.iteration = 4;[not work]

%0306-limit network simplex iteration-https://www.ibm.com/support/knowledgecenter/SSSA5P_12.9.0/ilog.odms.cplex.help/CPLEX/Parameters/topics/NetItLim.html
%Cplex.Param.network.iterations = 2;[not work]

%ops.mip.limits.nodes = 4;[not work]

%0304-use the online method for out-of-memory: https://www.ilovematlab.cn/thread-266244-1-1.html
%ops.mip.strategy.file = 3; %[not work: should put in sdpsetting] this has same results 1159110 node, but reduces time for 3 min.

%0305-use online method for out-of-memory: https://www.ibm.com/developerworks/community/forums/html/topic?id=e3ea344c-7249-4edd-a47d-29e1f80fa480
%Cplex.Param.mip.limits.auxrootthreads = 2; %[not work: should put in sdpsetting] -1 or 1

%obj-operational obj z_op
z_op = sum(sum(sum(cijt.*X(1:i_up,1:j_up,1:t_up)))) + sum(sum(cijt.*X_t0(1:i_up,1:j_up,1)))...
  + sum(sum(sum(pijt.*U(1:i_up,1:j_up,1:t_up)))) + sum(sum(pijt.*U_t0(1:i_up,1:j_up,1)))...
  + sum(sum(hit.*V(1:i_up,1:t_up))) + sum(hit.*V_t0(1:i_up,1));

%1.35-express D_ij_SAV
for i=1:i_up
  for j=1:j_up
    for t=1:t_up
      D_ij_SAV(i,j,t) = (b*d(i,j,t)+cijt*d(i,j,t)*m-n*T_abs*U(i,j,t))./(a+b+m*(cijt+P(i,j,t)));
    end
  end
end

% 0604-calculate minimum k ratio covered by SAVs
for i=1:i_up
  for j=1:j_up
    for t=1:t_up
      k(i,j,t) = D_ij_SAV(i,j,t)./d(i,j,t);
    end
  end
end

% k_min = min(k(1,1,1),k(1,2,1),k(1,3,1))

% k_min = 0.99;
% % 
% if ((k_min) >= (k(1,1,1)))
%   k_min = k(1,1,1)                
% end  

% 0604-find minimum covered ratio k
% for i=1:i_up
%   for j=1:j_up
%     for t=1:t_up
%       if (k_min >= k(i,j,t))
%         k_min = k(i,j,t)                
%       end  
%     end
%   end
% end

%1.4-check if solved and output
result = optimize(C,z_op,ops); % 0508-origin
% result = solvesdp(C,z_op,ops); % 0508-this works

%result = optimize(C,z);
% if result.problem == 0 % problem =0 means solve successfully,only print Uijt here

%   fprintf('\nk value: ') %ratio k covered by SAVs
%   value(k)

  
%   fprintf('\nP_Ini_(1,1,1) = ')
%   value(P(1,1,1))

%   fprintf('\nMin cost by Lower = ')
%   value(z_op)

  all_k = [];

  for i=1:i_up
    for j=1:j_up
      for t=1:t_up
        all_k(end+1) = k(i,j,t);
      end
    end
  end
  
  min_k = min(all_k);
  
  k_avg = mean(all_k);
  
%   fprintf('\n min k: ')
%   value(min_k)
  
  max_k = max(all_k);
  
%   fprintf('\n max k: ')
%   value(max_k)
  
  % 0615-save max_k in txt

%   fid = fopen(['IniP=',num2str(2),'max_k.txt'],'a'); %save to last 'a'
%   fprintf(fid,'%.4f\t',double(max_k)); %0516-keep 4 decimals
%   fclose(fid);
  
%   fprintf('\nU(1,1,1) value: ')
%   value(U(1,1,1))
        
%   fprintf('\nX(1,1,1) value: ') 
%   value(X(1,1,1))

fprintf('\n Iteration No. = ')
value(Iter_count)



% fprintf('\n D_ij_SAV(1,1,1) value: ') %0515-as for now hide D_ij_SAV
% value(D_ij_SAV(1,1,1))
  
%   fprintf('\n k(1,1,1) ')
%   value(k(1,1,1))

%   fprintf('\n')
%     fprintf('\nV value: ')
%     value(V)
%     
%     fprintf('\nU_t0 value: ')
%     value(U_t0)
%     
%     fprintf('\nX_t0 value: ') 
%     value(X_t0)
%     
%     fprintf('\nV_t0 value: ')
%     value(V_t0)

% else
%   disp('Not solved!')
% 
% end

% 0519-save total cost z_op in txt
% fid = fopen(['Iter=',num2str(Iter_count),'_OperCost.txt'],'wt'); %creat & write z_op to txt
% fprintf(fid,'%.4f\n',double(z_op)); %0516-keep 4 decimals
% fclose(fid);

% 0605-save min ratio k covered by SAVs in txt
% fid_ratiok = fopen(['Iter=',num2str(Iter_count),'_MinRatiok.txt'],'wt'); %creat & write z_op to txt
% fprintf(fid_ratiok,'%.4f\n',double(min_k)); %0604-keep 4 decimals
% fclose(fid_ratiok);

% 0517-following is save Uijt in txt
%1.6-save Uijt as txt for upper use
% for U_t=1:t_up 
%      fid = fopen(['IniP=',num2str(IniP),'Iter=',num2str(Iter_count),'_Uijt=',num2str(U_t),'.txt'],'wt'); %creat & write Uijt to txt
%      for i=1:i_up
%          for j=1:j_up
%              if j == j_up
%                  fprintf(fid,'%.4f\n',double(U(i,j,U_t))); %0516-keep 4 decimals
%              else
%                  fprintf(fid,'%.4f\t',double(U(i,j,U_t))); %0516-keep 4 decimals
%              end     
%          end    
%      end    
%      fclose(fid);
% end

% 0607-calculate Upper Pijt upper bound
for i=1:i_up
  for j=1:j_up
    for t=1:t_up
      P(i,j,t)=(b*d(i,j,t)+cijt*d(i,j,t)*m-n*T_abs*U(i,j,t)-(a+b)*min_k*d(i,j,t))./(m*min_k*d(i,j,t))-cijt;
    end
  end
end

% fprintf('\nUpper bound P(1,1,1) = %.4f\n',double(P(1,1,1)));

%0608-----
%calculate Upper revenue
%0608-----

revenue=0;

for i=1:i_up
  for j=1:j_up
    for t=1:t_up
      if (t==1)||(t==2)||(t==9)||(t==10)||(t==17)||(t==18)||(t==31)||(t==32)||(t==39)||(t==40)
%         T_abs=2;
        revenue = revenue-(P(i,j,t).*(b.*d(i,j,t)+cijt.*d(i,j,t)*m-n*T_abs*U(i,j,t)))./(a+b+m.*(cijt+P(i,j,t))); 
      else
%         T_abs=2; %here should be 3, but temporarily set to 2, aligh with Oper level
        revenue = revenue-(P(i,j,t).*(b.*d(i,j,t)+cijt.*d(i,j,t)*m-n*T_abs*U(i,j,t)))./(a+b+m.*(cijt+P(i,j,t)));
      end
    end
  end
end

% 070320-avg. all P

 all_P = [];

  for i=1:i_up
    for j=1:j_up
      for t=1:t_up
        all_P = P(i,j,t);
      end
    end
  end

avg_P = all_P./(i_up*j_up*t_up);



% 070320-avg. all D

 all_Dijt = [];

  for i=1:i_up
    for j=1:j_up
      for t=1:t_up
        all_Dijt = D_ij_SAV(i,j,t);
      end
    end
  end

avg_Dijt = all_Dijt./(i_up*j_up*t_up);

% 070320-print results

fprintf('\n Upper bound P(1,1,1) = %.4f\n',double(P(1,1,1))); % particular P each iteration

fprintf('\n');

fprintf('\n D(1,1,1) = %.4f\n',double(D_ij_SAV(1,1,1))); % particular Dijt each iteration

fprintf('\n');

fprintf('\n Avg_P = %.4f\n',double(avg_P)); % particular avg_P each iteration

fprintf('\n');

fprintf('\n Avg_Dijt = %.4f\n',double(avg_Dijt)); % particular avg_Dijt each iteration

fprintf('\n');

fprintf('\n Max Revenue by Upper= %.4f\n',double(-revenue));

fprintf('\n');

fprintf('\n Min Cost by Lower= %.4f\n',double(z_op));

fprintf('\n');

fprintf('\n Profit = %.4f\n',double(-revenue-z_op));

fprintf('\n');

% test how many new P > old P

satis_P = 0;

for i=1:i_up
  for j=1:j_up
    for t=1:t_up
%       if (double(P_sav(i,j,t))<=double(P(i,j,t))) % 0610-check new P greater than old ones
        %0610-if new P no more than 1 cent than old, than add to satis_P terminate condition
      %if ((double(P_sav(i,j,t)))*1.0001 >= double(P(i,j,t))) 
%       if ((double(P_sav(i,j,t)))+0.01 >= double(P(i,j,t))) 
      if ((P_sav(i,j,t))+0.01 > P(i,j,t)) 
        satis_P = (satis_P)+1;
      end  
    end
  end
end

satis_P;

fprintf('\n');

% test how many new P all feasible that new_P <= new P upper bound condition


if (Iter_count==1)
  feasib_P = 0;
else
  
  feasib_P = 0;
  
  for i=1:i_up
    for j=1:j_up
      for t=1:t_up
        P_bench(i,j,t) = (b*d(i,j,t)+cijt*d(i,j,t)*m-n*T_abs*U(i,j,t)-(a+b)*min_k*d(i,j,t))./(m*min_k*d(i,j,t))-cijt;        
      end
    end
  end
  
    
  for i=1:i_up
    for j=1:j_up
      for t=1:t_up
  %       if (double(P_sav(i,j,t))<=double(P(i,j,t))) % 0610-check new P greater than old ones
          %0610-if new P no more than 1 cent than old, than add to satis_P terminate condition
        %if ((double(P_sav(i,j,t)))*1.0001 >= double(P(i,j,t))) 
        
%         if (double(P_sav(i,j,t)) <= double(P_bench(i,j,t)))
        if (P_bench(i,j,t) >= P_sav(i,j,t))
%           i,j,t
%           fprintf("\n");
          feasib_P = (feasib_P)+1;
        end  
      end
    end
  end

end

feasib_P;

fprintf('\n');

if (satis_P==i_up*j_up*t_up) && (feasib_P==i_up*j_up*t_up)
  consec_P = 1+consec_P;
else
  consec_P = 0;
end

% tik = tik-1;

end

fprintf('\n');
fprintf('\n');
fprintf('\n');


%-------case 4.3------

'-------Case 4.3-----'

IniP=100;

%initially price Pijt all the same
P=zeros(20,20,40);%create 20*20*40 matrix
for i=1:i_up 
    for j=1:j_up
        for t=1:t_up
          P(i,j,t)=IniP; %0608-initial P should be infinity, I use 2 for test
        end
    end    
end

% Pijt=2;

Iter_count = 0;

% 0517-M is c16 time window, set at 1 gives it relaxation
% 0424-M is c16 time window, set at 3 gives it relaxation
M=1; 

%0515-change all var to continours instead of MIP
X = sdpvar(i_up,j_up,t_up+M); %0424-Xi,j,t+M, cuz time window constraint c16 Xijt is beyond t
%X = intvar(i_up,j_up,t_up+M); %0424-Xi,j,t+M, cuz time window constraint c16 Xijt is beyond t
X_t0 = sdpvar(i_up,j_up,1); %0424-Xi,j,0 means moves initially for c18: inventory balance eq
%X_t0 = intvar(i_up,j_up,1); %0424-Xi,j,0 means moves initially for c18: inventory balance eq
U = sdpvar(i_up,j_up,t_up);
%U = intvar(i_up,j_up,t_up); 
U_t0 = sdpvar(i_up,j_up,1); %0424-this is for Ui,j,0 
%U_t0 = intvar(i_up,j_up,1); %0424-this is for Ui,j,0 
V = sdpvar(i_up,t_up); 
%V = intvar(i_up,t_up); 
V_t0 = sdpvar(i_up,1); %0424-this is especially for Vi,0
%V_t0 = intvar(i_up,1); %0424-this is especially for Vi,0
D_ij_SAV = sdpvar(i_up,j_up,t_up); %0425-this is to show D_ij_SAV, for futher plot, then find upper bound of Pijt
%D_ij_SAV = intvar(i_up,j_up,t_up); %0425-this is to show D_ij_SAV, for futher plot, then find upper bound of Pijt
% 0604-k ratio
k = sdpvar(i_up,j_up,t_up);

% 070220-set place for Wijt
W=zeros(20,20,40);%create 20*20*40 matrix

%1.1-parameters
% a=b=m=n=1,T=40
a=1;
b=1;
m=1;
n=10;
T=8; %070120-total 40 time periods, 5 time block, so every block is 8
T_abs=8; %070120-T_abs should be same as T
% cijt=40, pijt=10, hit=10, alphaijt=0.1 (all c,p,h,alpha are subjectively to be same)
cijt=20;
pijt=10; % 070120-increase for penalty of unmet demand
% 0425-intentionally make pijt this small to avoid all 0 in Uijt solution
hit=10;
alpha=1;

%%%-----0607----
%%% How to add loop?
%%%----

% 0610-here starts the loop

P_sav = zeros(20,20,40);
U_sav = zeros(20,20,40);
min_k = 0;
P_bench = zeros(20,20,40);

% tik = 2;

satis_P = 0;
feasib_P = 0;
consec_P = 0;

% while (satis_P<i_up*j_up*t_up) || (feasib_P<i_up*j_up*t_up)
% while Iter_count < 2
while (consec_P<3)

Iter_count = Iter_count+1; 


for i=1:i_up 
    for j=1:j_up
        for t=1:t_up
          P_sav(i,j,t)=P(i,j,t); %0608-initial P should be infinity, I use 2 for test
        end
    end    
end  


for i=1:i_up 
    for j=1:j_up
        for t=1:t_up
            U_sav(i,j,t)=U(i,j,t); 
        end
    end    
end


min_k_sav = min_k;


%1.2-add operational level constraints

C = [];

% c15
for i=1:i_up %c15
   for j=1:j_up
     for t=1:t_up
       if (t==1)
         C = [C, U(i,j,1) >= U_t0(i,j,1)+(b*d(i,j,1)+cijt*d(i,j,1)*m-n*T_abs*U(i,j,1))./(a+b+m*(cijt+P(i,j,t)))-X(i,j,1)]; %c15: t=1
       else       
         C = [C, U(i,j,t) >= U(i,j,t-1)+(b*d(i,j,t)+cijt*d(i,j,t)*m-n*T_abs*U(i,j,t))./(a+b+m*(cijt+P(i,j,t)))-X(i,j,t)]; %c15: t from 2
       end
     end
   end    
end

% 0424-c16-I suppose time window is M period for c16 - relax this one, hope get solution
% this also considers Xijt beyond T
for i=1:i_up %c16
    for j=1:j_up
      for t=1:t_up
        if (t==1)
          C = [C, sum(X(i,j,t:(t+M))) >= U_t0(i,j,1)+(b*d(i,j,t)+cijt*d(i,j,t)*m-n*T_abs*U(i,j,t))./(a+b+m*(cijt+P(i,j,t)))]; %c16: t=1
        else
          C = [C, sum(X(i,j,t:(t+M))) >= U(i,j,t-1)+(b*d(i,j,t)+cijt*d(i,j,t)*m-n*T_abs*U(i,j,t))./(a+b+m*(cijt+P(i,j,t)))]; %c16: t from 2
        end
      end  
    end     
end

% 0424 - c17 is in upper level, no need in lower level
% for i=1:i_up %c17
%   for j=1:j_up
%     for t=1:t_up
%       C = [C, ((b-1)*d(i,j,t)+cijt*d(i,j,t)*m)./(n*T_abs) <= U(i,j,t)];
%       C = [C, U(i,j,t) <= (b*d(i,j,t)+cijt*d(i,j,t)*m)./(n*T_abs)];
%     end
%   end
% end

% 0424-c18-here makes <= for Vi,t cuz alpha all set to 1, too relaxed, should be variaty
for i=1:i_up %c18
  for t=1:t_up
    if (t==1)
      C = [C, V(i,t) <= V_t0(i,1) + sum(alpha*X_t0(1:j_up,i,1)) - sum(X(i,1:j_up,t)) ]; %c18: t=1
    else
      C = [C, V(i,t) <= V(i,t-1) + sum(sum(alpha*X(1:j_up,i,1:(t-1)))) - sum(X(i,1:j_up,t)) ]; %c18: t from 2 (remember to use 2 sum for 2nd)
    end
  end
end

% 0508-c19-not added, cuz it looks already relaxed in math function, hidden property of our formula

for i=1:i_up %c20: Xijt
    for j=1:j_up
        for t=1:(t_up+M)
          C = [C, X(i,j,t) >= 0];
        end
    end    
end    

for i=1:i_up %c20: X_t0(i,j,1) is Xijt=0
    for j=1:j_up
        C = [C, X_t0(i,j,1) >= 0];
    end    
end    

for i=1:i_up %c20: Uijt
    for j=1:j_up
        for t=1:t_up
            C = [C, U(i,j,t) >= 0];
        end
    end    
end    

for i=1:i_up %c20: U_t0(i,j,1) is Uijt=0
    for j=1:j_up
        C = [C, U_t0(i,j,1) >= 0];
    end    
end 

for i=1:i_up %c20: Vit
    for t=1:t_up
      C = [C, V(i,t) >= 0];
    end    
end    

for i=1:i_up %c20: V_t0(i,1) is V(i,0)
  C = [C, V_t0(i,1) >= 0];
end

%1.3-solve operational level: settings: use cplex as solver

%0515-following is set optimality for MIP, but kang said make it LNContinuous-https://www.ibm.com/support/knowledgecenter/SSSA5P_12.9.0/ilog.odms.cplex.help/CPLEX/Parameters/topics/OptimalityTarget.html
%0515-'cplex.optimalitytarget'=1: Searches for a globally optimal solution to a convex model.
ops = sdpsettings ('solver','cplex','verbose',2,'cplex.optimalitytarget',1);%0508-verbose 2 shows iterations, 1 means shown normal message, 0 shows silent message
%0515-MIP limit tolerance gap to 1-https://www.ibm.com/support/knowledgecenter/SSSA5P_12.9.0/ilog.odms.cplex.help/CPLEX/Parameters/topics/EpGap.html
%ops = sdpsettings ('solver','cplex','verbose',2,'cplex.mip.tolerances.mipgap',1);%0508-verbose 2 shows iterations, 1 means shown normal message, 0 shows silent message

%0306-limit barrier iteration-https://www.ibm.com/support/knowledgecenter/SSSA5P_12.9.0/ilog.odms.cplex.help/CPLEX/Parameters/topics/BarItLim.html
%Cplex.Param.barrier.limits.iteration = 4;[not work]

%0306-limit network simplex iteration-https://www.ibm.com/support/knowledgecenter/SSSA5P_12.9.0/ilog.odms.cplex.help/CPLEX/Parameters/topics/NetItLim.html
%Cplex.Param.network.iterations = 2;[not work]

%ops.mip.limits.nodes = 4;[not work]

%0304-use the online method for out-of-memory: https://www.ilovematlab.cn/thread-266244-1-1.html
%ops.mip.strategy.file = 3; %[not work: should put in sdpsetting] this has same results 1159110 node, but reduces time for 3 min.

%0305-use online method for out-of-memory: https://www.ibm.com/developerworks/community/forums/html/topic?id=e3ea344c-7249-4edd-a47d-29e1f80fa480
%Cplex.Param.mip.limits.auxrootthreads = 2; %[not work: should put in sdpsetting] -1 or 1

%obj-operational obj z_op
z_op = sum(sum(sum(cijt.*X(1:i_up,1:j_up,1:t_up)))) + sum(sum(cijt.*X_t0(1:i_up,1:j_up,1)))...
  + sum(sum(sum(pijt.*U(1:i_up,1:j_up,1:t_up)))) + sum(sum(pijt.*U_t0(1:i_up,1:j_up,1)))...
  + sum(sum(hit.*V(1:i_up,1:t_up))) + sum(hit.*V_t0(1:i_up,1));

%1.35-express D_ij_SAV
for i=1:i_up
  for j=1:j_up
    for t=1:t_up
      D_ij_SAV(i,j,t) = (b*d(i,j,t)+cijt*d(i,j,t)*m-n*T_abs*U(i,j,t))./(a+b+m*(cijt+P(i,j,t)));
    end
  end
end

% 0604-calculate minimum k ratio covered by SAVs
for i=1:i_up
  for j=1:j_up
    for t=1:t_up
      k(i,j,t) = D_ij_SAV(i,j,t)./d(i,j,t);
    end
  end
end

% k_min = min(k(1,1,1),k(1,2,1),k(1,3,1))

% k_min = 0.99;
% % 
% if ((k_min) >= (k(1,1,1)))
%   k_min = k(1,1,1)                
% end  

% 0604-find minimum covered ratio k
% for i=1:i_up
%   for j=1:j_up
%     for t=1:t_up
%       if (k_min >= k(i,j,t))
%         k_min = k(i,j,t)                
%       end  
%     end
%   end
% end

%1.4-check if solved and output
result = optimize(C,z_op,ops); % 0508-origin
% result = solvesdp(C,z_op,ops); % 0508-this works

%result = optimize(C,z);
% if result.problem == 0 % problem =0 means solve successfully,only print Uijt here

%   fprintf('\nk value: ') %ratio k covered by SAVs
%   value(k)

  
%   fprintf('\nP_Ini_(1,1,1) = ')
%   value(P(1,1,1))

%   fprintf('\nMin cost by Lower = ')
%   value(z_op)

  all_k = [];

  for i=1:i_up
    for j=1:j_up
      for t=1:t_up
        all_k(end+1) = k(i,j,t);
      end
    end
  end
  
  min_k = min(all_k);
  
  k_avg = mean(all_k);
  
%   fprintf('\n min k: ')
%   value(min_k)
  
  max_k = max(all_k);
  
%   fprintf('\n max k: ')
%   value(max_k)
  
  % 0615-save max_k in txt

%   fid = fopen(['IniP=',num2str(2),'max_k.txt'],'a'); %save to last 'a'
%   fprintf(fid,'%.4f\t',double(max_k)); %0516-keep 4 decimals
%   fclose(fid);
  
%   fprintf('\nU(1,1,1) value: ')
%   value(U(1,1,1))
        
%   fprintf('\nX(1,1,1) value: ') 
%   value(X(1,1,1))

fprintf('\n Iteration No. = ')
value(Iter_count)



% fprintf('\n D_ij_SAV(1,1,1) value: ') %0515-as for now hide D_ij_SAV
% value(D_ij_SAV(1,1,1))
  
%   fprintf('\n k(1,1,1) ')
%   value(k(1,1,1))

%   fprintf('\n')
%     fprintf('\nV value: ')
%     value(V)
%     
%     fprintf('\nU_t0 value: ')
%     value(U_t0)
%     
%     fprintf('\nX_t0 value: ') 
%     value(X_t0)
%     
%     fprintf('\nV_t0 value: ')
%     value(V_t0)

% else
%   disp('Not solved!')
% 
% end

% 0519-save total cost z_op in txt
% fid = fopen(['Iter=',num2str(Iter_count),'_OperCost.txt'],'wt'); %creat & write z_op to txt
% fprintf(fid,'%.4f\n',double(z_op)); %0516-keep 4 decimals
% fclose(fid);

% 0605-save min ratio k covered by SAVs in txt
% fid_ratiok = fopen(['Iter=',num2str(Iter_count),'_MinRatiok.txt'],'wt'); %creat & write z_op to txt
% fprintf(fid_ratiok,'%.4f\n',double(min_k)); %0604-keep 4 decimals
% fclose(fid_ratiok);

% 0517-following is save Uijt in txt
%1.6-save Uijt as txt for upper use
% for U_t=1:t_up 
%      fid = fopen(['IniP=',num2str(IniP),'Iter=',num2str(Iter_count),'_Uijt=',num2str(U_t),'.txt'],'wt'); %creat & write Uijt to txt
%      for i=1:i_up
%          for j=1:j_up
%              if j == j_up
%                  fprintf(fid,'%.4f\n',double(U(i,j,U_t))); %0516-keep 4 decimals
%              else
%                  fprintf(fid,'%.4f\t',double(U(i,j,U_t))); %0516-keep 4 decimals
%              end     
%          end    
%      end    
%      fclose(fid);
% end

% 0607-calculate Upper Pijt upper bound
for i=1:i_up
  for j=1:j_up
    for t=1:t_up
      P(i,j,t)=(b*d(i,j,t)+cijt*d(i,j,t)*m-n*T_abs*U(i,j,t)-(a+b)*min_k*d(i,j,t))./(m*min_k*d(i,j,t))-cijt;
    end
  end
end

% fprintf('\nUpper bound P(1,1,1) = %.4f\n',double(P(1,1,1)));

%0608-----
%calculate Upper revenue
%0608-----

revenue=0;

for i=1:i_up
  for j=1:j_up
    for t=1:t_up
      if (t==1)||(t==2)||(t==9)||(t==10)||(t==17)||(t==18)||(t==31)||(t==32)||(t==39)||(t==40)
%         T_abs=2;
        revenue = revenue-(P(i,j,t).*(b.*d(i,j,t)+cijt.*d(i,j,t)*m-n*T_abs*U(i,j,t)))./(a+b+m.*(cijt+P(i,j,t))); 
      else
%         T_abs=2; %here should be 3, but temporarily set to 2, aligh with Oper level
        revenue = revenue-(P(i,j,t).*(b.*d(i,j,t)+cijt.*d(i,j,t)*m-n*T_abs*U(i,j,t)))./(a+b+m.*(cijt+P(i,j,t)));
      end
    end
  end
end

% 070320-avg. all P

 all_P = [];

  for i=1:i_up
    for j=1:j_up
      for t=1:t_up
        all_P = P(i,j,t);
      end
    end
  end

avg_P = all_P./(i_up*j_up*t_up);



% 070320-avg. all D

 all_Dijt = [];

  for i=1:i_up
    for j=1:j_up
      for t=1:t_up
        all_Dijt = D_ij_SAV(i,j,t);
      end
    end
  end

avg_Dijt = all_Dijt./(i_up*j_up*t_up);

% 070320-print results

fprintf('\n Upper bound P(1,1,1) = %.4f\n',double(P(1,1,1))); % particular P each iteration

fprintf('\n');

fprintf('\n D(1,1,1) = %.4f\n',double(D_ij_SAV(1,1,1))); % particular Dijt each iteration

fprintf('\n');

fprintf('\n Avg_P = %.4f\n',double(avg_P)); % particular avg_P each iteration

fprintf('\n');

fprintf('\n Avg_Dijt = %.4f\n',double(avg_Dijt)); % particular avg_Dijt each iteration

fprintf('\n');

fprintf('\n Max Revenue by Upper= %.4f\n',double(-revenue));

fprintf('\n');

fprintf('\n Min Cost by Lower= %.4f\n',double(z_op));

fprintf('\n');

fprintf('\n Profit = %.4f\n',double(-revenue-z_op));

fprintf('\n');

% test how many new P > old P

satis_P = 0;

for i=1:i_up
  for j=1:j_up
    for t=1:t_up
%       if (double(P_sav(i,j,t))<=double(P(i,j,t))) % 0610-check new P greater than old ones
        %0610-if new P no more than 1 cent than old, than add to satis_P terminate condition
      %if ((double(P_sav(i,j,t)))*1.0001 >= double(P(i,j,t))) 
%       if ((double(P_sav(i,j,t)))+0.01 >= double(P(i,j,t))) 
      if ((P_sav(i,j,t))+0.01 > P(i,j,t)) 
        satis_P = (satis_P)+1;
      end  
    end
  end
end

satis_P;

fprintf('\n');

% test how many new P all feasible that new_P <= new P upper bound condition


if (Iter_count==1)
  feasib_P = 0;
else
  
  feasib_P = 0;
  
  for i=1:i_up
    for j=1:j_up
      for t=1:t_up
        P_bench(i,j,t) = (b*d(i,j,t)+cijt*d(i,j,t)*m-n*T_abs*U(i,j,t)-(a+b)*min_k*d(i,j,t))./(m*min_k*d(i,j,t))-cijt;        
      end
    end
  end
  
    
  for i=1:i_up
    for j=1:j_up
      for t=1:t_up
  %       if (double(P_sav(i,j,t))<=double(P(i,j,t))) % 0610-check new P greater than old ones
          %0610-if new P no more than 1 cent than old, than add to satis_P terminate condition
        %if ((double(P_sav(i,j,t)))*1.0001 >= double(P(i,j,t))) 
        
%         if (double(P_sav(i,j,t)) <= double(P_bench(i,j,t)))
        if (P_bench(i,j,t) >= P_sav(i,j,t))
%           i,j,t
%           fprintf("\n");
          feasib_P = (feasib_P)+1;
        end  
      end
    end
  end

end

feasib_P;

fprintf('\n');

if (satis_P==i_up*j_up*t_up) && (feasib_P==i_up*j_up*t_up)
  consec_P = 1+consec_P;
else
  consec_P = 0;
end

% tik = tik-1;

end

fprintf('\n');
fprintf('\n');
fprintf('\n');



%-------case 5.1------

'------Case 5.1------'

IniP=1;

%initially price Pijt all the same
P=zeros(20,20,40);%create 20*20*40 matrix
for i=1:i_up 
    for j=1:j_up
        for t=1:t_up
          P(i,j,t)=IniP; %0608-initial P should be infinity, I use 2 for test
        end
    end    
end

% Pijt=2;

Iter_count = 0;

% 0517-M is c16 time window, set at 1 gives it relaxation
% 0424-M is c16 time window, set at 3 gives it relaxation
M=1; 

%0515-change all var to continours instead of MIP
X = sdpvar(i_up,j_up,t_up+M); %0424-Xi,j,t+M, cuz time window constraint c16 Xijt is beyond t
%X = intvar(i_up,j_up,t_up+M); %0424-Xi,j,t+M, cuz time window constraint c16 Xijt is beyond t
X_t0 = sdpvar(i_up,j_up,1); %0424-Xi,j,0 means moves initially for c18: inventory balance eq
%X_t0 = intvar(i_up,j_up,1); %0424-Xi,j,0 means moves initially for c18: inventory balance eq
U = sdpvar(i_up,j_up,t_up);
%U = intvar(i_up,j_up,t_up); 
U_t0 = sdpvar(i_up,j_up,1); %0424-this is for Ui,j,0 
%U_t0 = intvar(i_up,j_up,1); %0424-this is for Ui,j,0 
V = sdpvar(i_up,t_up); 
%V = intvar(i_up,t_up); 
V_t0 = sdpvar(i_up,1); %0424-this is especially for Vi,0
%V_t0 = intvar(i_up,1); %0424-this is especially for Vi,0
D_ij_SAV = sdpvar(i_up,j_up,t_up); %0425-this is to show D_ij_SAV, for futher plot, then find upper bound of Pijt
%D_ij_SAV = intvar(i_up,j_up,t_up); %0425-this is to show D_ij_SAV, for futher plot, then find upper bound of Pijt
% 0604-k ratio
k = sdpvar(i_up,j_up,t_up);

% 070220-set place for Wijt
W=zeros(20,20,40);%create 20*20*40 matrix

%1.1-parameters
% a=b=m=n=1,T=40
a=1;
b=1;
m=1;
n=1;
T=8; %070120-total 40 time periods, 5 time block, so every block is 8
T_abs=8; %070120-T_abs should be same as T
% cijt=40, pijt=10, hit=10, alphaijt=0.1 (all c,p,h,alpha are subjectively to be same)
cijt=40;
pijt=10; % 070120-increase for penalty of unmet demand
% 0425-intentionally make pijt this small to avoid all 0 in Uijt solution
hit=10;
alpha=1;

%%%-----0607----
%%% How to add loop?
%%%----

% 0610-here starts the loop

P_sav = zeros(20,20,40);
U_sav = zeros(20,20,40);
min_k = 0;
P_bench = zeros(20,20,40);

% tik = 2;

satis_P = 0;
feasib_P = 0;
consec_P = 0;

% while (satis_P<i_up*j_up*t_up) || (feasib_P<i_up*j_up*t_up)
% while Iter_count < 2
while (consec_P<3)

Iter_count = Iter_count+1; 


for i=1:i_up 
    for j=1:j_up
        for t=1:t_up
          P_sav(i,j,t)=P(i,j,t); %0608-initial P should be infinity, I use 2 for test
        end
    end    
end  


for i=1:i_up 
    for j=1:j_up
        for t=1:t_up
            U_sav(i,j,t)=U(i,j,t); 
        end
    end    
end


min_k_sav = min_k;


%1.2-add operational level constraints

C = [];

% c15
for i=1:i_up %c15
   for j=1:j_up
     for t=1:t_up
       if (t==1)
         C = [C, U(i,j,1) >= U_t0(i,j,1)+(b*d(i,j,1)+cijt*d(i,j,1)*m-n*T_abs*U(i,j,1))./(a+b+m*(cijt+P(i,j,t)))-X(i,j,1)]; %c15: t=1
       else       
         C = [C, U(i,j,t) >= U(i,j,t-1)+(b*d(i,j,t)+cijt*d(i,j,t)*m-n*T_abs*U(i,j,t))./(a+b+m*(cijt+P(i,j,t)))-X(i,j,t)]; %c15: t from 2
       end
     end
   end    
end

% 0424-c16-I suppose time window is M period for c16 - relax this one, hope get solution
% this also considers Xijt beyond T
for i=1:i_up %c16
    for j=1:j_up
      for t=1:t_up
        if (t==1)
          C = [C, sum(X(i,j,t:(t+M))) >= U_t0(i,j,1)+(b*d(i,j,t)+cijt*d(i,j,t)*m-n*T_abs*U(i,j,t))./(a+b+m*(cijt+P(i,j,t)))]; %c16: t=1
        else
          C = [C, sum(X(i,j,t:(t+M))) >= U(i,j,t-1)+(b*d(i,j,t)+cijt*d(i,j,t)*m-n*T_abs*U(i,j,t))./(a+b+m*(cijt+P(i,j,t)))]; %c16: t from 2
        end
      end  
    end     
end

% 0424 - c17 is in upper level, no need in lower level
% for i=1:i_up %c17
%   for j=1:j_up
%     for t=1:t_up
%       C = [C, ((b-1)*d(i,j,t)+cijt*d(i,j,t)*m)./(n*T_abs) <= U(i,j,t)];
%       C = [C, U(i,j,t) <= (b*d(i,j,t)+cijt*d(i,j,t)*m)./(n*T_abs)];
%     end
%   end
% end

% 0424-c18-here makes <= for Vi,t cuz alpha all set to 1, too relaxed, should be variaty
for i=1:i_up %c18
  for t=1:t_up
    if (t==1)
      C = [C, V(i,t) <= V_t0(i,1) + sum(alpha*X_t0(1:j_up,i,1)) - sum(X(i,1:j_up,t)) ]; %c18: t=1
    else
      C = [C, V(i,t) <= V(i,t-1) + sum(sum(alpha*X(1:j_up,i,1:(t-1)))) - sum(X(i,1:j_up,t)) ]; %c18: t from 2 (remember to use 2 sum for 2nd)
    end
  end
end

% 0508-c19-not added, cuz it looks already relaxed in math function, hidden property of our formula

for i=1:i_up %c20: Xijt
    for j=1:j_up
        for t=1:(t_up+M)
          C = [C, X(i,j,t) >= 0];
        end
    end    
end    

for i=1:i_up %c20: X_t0(i,j,1) is Xijt=0
    for j=1:j_up
        C = [C, X_t0(i,j,1) >= 0];
    end    
end    

for i=1:i_up %c20: Uijt
    for j=1:j_up
        for t=1:t_up
            C = [C, U(i,j,t) >= 0];
        end
    end    
end    

for i=1:i_up %c20: U_t0(i,j,1) is Uijt=0
    for j=1:j_up
        C = [C, U_t0(i,j,1) >= 0];
    end    
end 

for i=1:i_up %c20: Vit
    for t=1:t_up
      C = [C, V(i,t) >= 0];
    end    
end    

for i=1:i_up %c20: V_t0(i,1) is V(i,0)
  C = [C, V_t0(i,1) >= 0];
end

%1.3-solve operational level: settings: use cplex as solver

%0515-following is set optimality for MIP, but kang said make it LNContinuous-https://www.ibm.com/support/knowledgecenter/SSSA5P_12.9.0/ilog.odms.cplex.help/CPLEX/Parameters/topics/OptimalityTarget.html
%0515-'cplex.optimalitytarget'=1: Searches for a globally optimal solution to a convex model.
ops = sdpsettings ('solver','cplex','verbose',2,'cplex.optimalitytarget',1);%0508-verbose 2 shows iterations, 1 means shown normal message, 0 shows silent message
%0515-MIP limit tolerance gap to 1-https://www.ibm.com/support/knowledgecenter/SSSA5P_12.9.0/ilog.odms.cplex.help/CPLEX/Parameters/topics/EpGap.html
%ops = sdpsettings ('solver','cplex','verbose',2,'cplex.mip.tolerances.mipgap',1);%0508-verbose 2 shows iterations, 1 means shown normal message, 0 shows silent message

%0306-limit barrier iteration-https://www.ibm.com/support/knowledgecenter/SSSA5P_12.9.0/ilog.odms.cplex.help/CPLEX/Parameters/topics/BarItLim.html
%Cplex.Param.barrier.limits.iteration = 4;[not work]

%0306-limit network simplex iteration-https://www.ibm.com/support/knowledgecenter/SSSA5P_12.9.0/ilog.odms.cplex.help/CPLEX/Parameters/topics/NetItLim.html
%Cplex.Param.network.iterations = 2;[not work]

%ops.mip.limits.nodes = 4;[not work]

%0304-use the online method for out-of-memory: https://www.ilovematlab.cn/thread-266244-1-1.html
%ops.mip.strategy.file = 3; %[not work: should put in sdpsetting] this has same results 1159110 node, but reduces time for 3 min.

%0305-use online method for out-of-memory: https://www.ibm.com/developerworks/community/forums/html/topic?id=e3ea344c-7249-4edd-a47d-29e1f80fa480
%Cplex.Param.mip.limits.auxrootthreads = 2; %[not work: should put in sdpsetting] -1 or 1

%obj-operational obj z_op
z_op = sum(sum(sum(cijt.*X(1:i_up,1:j_up,1:t_up)))) + sum(sum(cijt.*X_t0(1:i_up,1:j_up,1)))...
  + sum(sum(sum(pijt.*U(1:i_up,1:j_up,1:t_up)))) + sum(sum(pijt.*U_t0(1:i_up,1:j_up,1)))...
  + sum(sum(hit.*V(1:i_up,1:t_up))) + sum(hit.*V_t0(1:i_up,1));

%1.35-express D_ij_SAV
for i=1:i_up
  for j=1:j_up
    for t=1:t_up
      D_ij_SAV(i,j,t) = (b*d(i,j,t)+cijt*d(i,j,t)*m-n*T_abs*U(i,j,t))./(a+b+m*(cijt+P(i,j,t)));
    end
  end
end

% 0604-calculate minimum k ratio covered by SAVs
for i=1:i_up
  for j=1:j_up
    for t=1:t_up
      k(i,j,t) = D_ij_SAV(i,j,t)./d(i,j,t);
    end
  end
end

% k_min = min(k(1,1,1),k(1,2,1),k(1,3,1))

% k_min = 0.99;
% % 
% if ((k_min) >= (k(1,1,1)))
%   k_min = k(1,1,1)                
% end  

% 0604-find minimum covered ratio k
% for i=1:i_up
%   for j=1:j_up
%     for t=1:t_up
%       if (k_min >= k(i,j,t))
%         k_min = k(i,j,t)                
%       end  
%     end
%   end
% end

%1.4-check if solved and output
result = optimize(C,z_op,ops); % 0508-origin
% result = solvesdp(C,z_op,ops); % 0508-this works

%result = optimize(C,z);
% if result.problem == 0 % problem =0 means solve successfully,only print Uijt here

%   fprintf('\nk value: ') %ratio k covered by SAVs
%   value(k)

  
%   fprintf('\nP_Ini_(1,1,1) = ')
%   value(P(1,1,1))

%   fprintf('\nMin cost by Lower = ')
%   value(z_op)

  all_k = [];

  for i=1:i_up
    for j=1:j_up
      for t=1:t_up
        all_k(end+1) = k(i,j,t);
      end
    end
  end
  
  min_k = min(all_k);
  
  k_avg = mean(all_k);
  
%   fprintf('\n min k: ')
%   value(min_k)
  
  max_k = max(all_k);
  
%   fprintf('\n max k: ')
%   value(max_k)
  
  % 0615-save max_k in txt

%   fid = fopen(['IniP=',num2str(2),'max_k.txt'],'a'); %save to last 'a'
%   fprintf(fid,'%.4f\t',double(max_k)); %0516-keep 4 decimals
%   fclose(fid);
  
%   fprintf('\nU(1,1,1) value: ')
%   value(U(1,1,1))
        
%   fprintf('\nX(1,1,1) value: ') 
%   value(X(1,1,1))

fprintf('\n Iteration No. = ')
value(Iter_count)



% fprintf('\n D_ij_SAV(1,1,1) value: ') %0515-as for now hide D_ij_SAV
% value(D_ij_SAV(1,1,1))
  
%   fprintf('\n k(1,1,1) ')
%   value(k(1,1,1))

%   fprintf('\n')
%     fprintf('\nV value: ')
%     value(V)
%     
%     fprintf('\nU_t0 value: ')
%     value(U_t0)
%     
%     fprintf('\nX_t0 value: ') 
%     value(X_t0)
%     
%     fprintf('\nV_t0 value: ')
%     value(V_t0)

% else
%   disp('Not solved!')
% 
% end

% 0519-save total cost z_op in txt
% fid = fopen(['Iter=',num2str(Iter_count),'_OperCost.txt'],'wt'); %creat & write z_op to txt
% fprintf(fid,'%.4f\n',double(z_op)); %0516-keep 4 decimals
% fclose(fid);

% 0605-save min ratio k covered by SAVs in txt
% fid_ratiok = fopen(['Iter=',num2str(Iter_count),'_MinRatiok.txt'],'wt'); %creat & write z_op to txt
% fprintf(fid_ratiok,'%.4f\n',double(min_k)); %0604-keep 4 decimals
% fclose(fid_ratiok);

% 0517-following is save Uijt in txt
%1.6-save Uijt as txt for upper use
% for U_t=1:t_up 
%      fid = fopen(['IniP=',num2str(IniP),'Iter=',num2str(Iter_count),'_Uijt=',num2str(U_t),'.txt'],'wt'); %creat & write Uijt to txt
%      for i=1:i_up
%          for j=1:j_up
%              if j == j_up
%                  fprintf(fid,'%.4f\n',double(U(i,j,U_t))); %0516-keep 4 decimals
%              else
%                  fprintf(fid,'%.4f\t',double(U(i,j,U_t))); %0516-keep 4 decimals
%              end     
%          end    
%      end    
%      fclose(fid);
% end

% 0607-calculate Upper Pijt upper bound
for i=1:i_up
  for j=1:j_up
    for t=1:t_up
      P(i,j,t)=(b*d(i,j,t)+cijt*d(i,j,t)*m-n*T_abs*U(i,j,t)-(a+b)*min_k*d(i,j,t))./(m*min_k*d(i,j,t))-cijt;
    end
  end
end

% fprintf('\nUpper bound P(1,1,1) = %.4f\n',double(P(1,1,1)));

%0608-----
%calculate Upper revenue
%0608-----

revenue=0;

for i=1:i_up
  for j=1:j_up
    for t=1:t_up
      if (t==1)||(t==2)||(t==9)||(t==10)||(t==17)||(t==18)||(t==31)||(t==32)||(t==39)||(t==40)
%         T_abs=2;
        revenue = revenue-(P(i,j,t).*(b.*d(i,j,t)+cijt.*d(i,j,t)*m-n*T_abs*U(i,j,t)))./(a+b+m.*(cijt+P(i,j,t))); 
      else
%         T_abs=2; %here should be 3, but temporarily set to 2, aligh with Oper level
        revenue = revenue-(P(i,j,t).*(b.*d(i,j,t)+cijt.*d(i,j,t)*m-n*T_abs*U(i,j,t)))./(a+b+m.*(cijt+P(i,j,t)));
      end
    end
  end
end

% 070320-avg. all P

 all_P = [];

  for i=1:i_up
    for j=1:j_up
      for t=1:t_up
        all_P = P(i,j,t);
      end
    end
  end

avg_P = all_P./(i_up*j_up*t_up);



% 070320-avg. all D

 all_Dijt = [];

  for i=1:i_up
    for j=1:j_up
      for t=1:t_up
        all_Dijt = D_ij_SAV(i,j,t);
      end
    end
  end

avg_Dijt = all_Dijt./(i_up*j_up*t_up);

% 070320-print results

fprintf('\n Upper bound P(1,1,1) = %.4f\n',double(P(1,1,1))); % particular P each iteration

fprintf('\n');

fprintf('\n D(1,1,1) = %.4f\n',double(D_ij_SAV(1,1,1))); % particular Dijt each iteration

fprintf('\n');

fprintf('\n Avg_P = %.4f\n',double(avg_P)); % particular avg_P each iteration

fprintf('\n');

fprintf('\n Avg_Dijt = %.4f\n',double(avg_Dijt)); % particular avg_Dijt each iteration

fprintf('\n');

fprintf('\n Max Revenue by Upper= %.4f\n',double(-revenue));

fprintf('\n');

fprintf('\n Min Cost by Lower= %.4f\n',double(z_op));

fprintf('\n');

fprintf('\n Profit = %.4f\n',double(-revenue-z_op));

fprintf('\n');

% test how many new P > old P

satis_P = 0;

for i=1:i_up
  for j=1:j_up
    for t=1:t_up
%       if (double(P_sav(i,j,t))<=double(P(i,j,t))) % 0610-check new P greater than old ones
        %0610-if new P no more than 1 cent than old, than add to satis_P terminate condition
      %if ((double(P_sav(i,j,t)))*1.0001 >= double(P(i,j,t))) 
%       if ((double(P_sav(i,j,t)))+0.01 >= double(P(i,j,t))) 
      if ((P_sav(i,j,t))+0.01 > P(i,j,t)) 
        satis_P = (satis_P)+1;
      end  
    end
  end
end

satis_P;

fprintf('\n');

% test how many new P all feasible that new_P <= new P upper bound condition


if (Iter_count==1)
  feasib_P = 0;
else
  
  feasib_P = 0;
  
  for i=1:i_up
    for j=1:j_up
      for t=1:t_up
        P_bench(i,j,t) = (b*d(i,j,t)+cijt*d(i,j,t)*m-n*T_abs*U(i,j,t)-(a+b)*min_k*d(i,j,t))./(m*min_k*d(i,j,t))-cijt;        
      end
    end
  end
  
    
  for i=1:i_up
    for j=1:j_up
      for t=1:t_up
  %       if (double(P_sav(i,j,t))<=double(P(i,j,t))) % 0610-check new P greater than old ones
          %0610-if new P no more than 1 cent than old, than add to satis_P terminate condition
        %if ((double(P_sav(i,j,t)))*1.0001 >= double(P(i,j,t))) 
        
%         if (double(P_sav(i,j,t)) <= double(P_bench(i,j,t)))
        if (P_bench(i,j,t) >= P_sav(i,j,t))
%           i,j,t
%           fprintf("\n");
          feasib_P = (feasib_P)+1;
        end  
      end
    end
  end

end

feasib_P;

fprintf('\n');

if (satis_P==i_up*j_up*t_up) && (feasib_P==i_up*j_up*t_up)
  consec_P = 1+consec_P;
else
  consec_P = 0;
end

% tik = tik-1;

end

fprintf('\n');
fprintf('\n');
fprintf('\n');


%-------case 5.2------

'-------Case 5.2-------'

IniP=10;

%initially price Pijt all the same
P=zeros(20,20,40);%create 20*20*40 matrix
for i=1:i_up 
    for j=1:j_up
        for t=1:t_up
          P(i,j,t)=IniP; %0608-initial P should be infinity, I use 2 for test
        end
    end    
end

% Pijt=2;

Iter_count = 0;

% 0517-M is c16 time window, set at 1 gives it relaxation
% 0424-M is c16 time window, set at 3 gives it relaxation
M=1; 

%0515-change all var to continours instead of MIP
X = sdpvar(i_up,j_up,t_up+M); %0424-Xi,j,t+M, cuz time window constraint c16 Xijt is beyond t
%X = intvar(i_up,j_up,t_up+M); %0424-Xi,j,t+M, cuz time window constraint c16 Xijt is beyond t
X_t0 = sdpvar(i_up,j_up,1); %0424-Xi,j,0 means moves initially for c18: inventory balance eq
%X_t0 = intvar(i_up,j_up,1); %0424-Xi,j,0 means moves initially for c18: inventory balance eq
U = sdpvar(i_up,j_up,t_up);
%U = intvar(i_up,j_up,t_up); 
U_t0 = sdpvar(i_up,j_up,1); %0424-this is for Ui,j,0 
%U_t0 = intvar(i_up,j_up,1); %0424-this is for Ui,j,0 
V = sdpvar(i_up,t_up); 
%V = intvar(i_up,t_up); 
V_t0 = sdpvar(i_up,1); %0424-this is especially for Vi,0
%V_t0 = intvar(i_up,1); %0424-this is especially for Vi,0
D_ij_SAV = sdpvar(i_up,j_up,t_up); %0425-this is to show D_ij_SAV, for futher plot, then find upper bound of Pijt
%D_ij_SAV = intvar(i_up,j_up,t_up); %0425-this is to show D_ij_SAV, for futher plot, then find upper bound of Pijt
% 0604-k ratio
k = sdpvar(i_up,j_up,t_up);

% 070220-set place for Wijt
W=zeros(20,20,40);%create 20*20*40 matrix

%1.1-parameters
% a=b=m=n=1,T=40
a=1;
b=1;
m=1;
n=1;
T=8; %070120-total 40 time periods, 5 time block, so every block is 8
T_abs=8; %070120-T_abs should be same as T
% cijt=40, pijt=10, hit=10, alphaijt=0.1 (all c,p,h,alpha are subjectively to be same)
cijt=40;
pijt=10; % 070120-increase for penalty of unmet demand
% 0425-intentionally make pijt this small to avoid all 0 in Uijt solution
hit=10;
alpha=1;

%%%-----0607----
%%% How to add loop?
%%%----

% 0610-here starts the loop

P_sav = zeros(20,20,40);
U_sav = zeros(20,20,40);
min_k = 0;
P_bench = zeros(20,20,40);

% tik = 2;

satis_P = 0;
feasib_P = 0;
consec_P = 0;

% while (satis_P<i_up*j_up*t_up) || (feasib_P<i_up*j_up*t_up)
% while Iter_count < 2
while (consec_P<3)

Iter_count = Iter_count+1; 


for i=1:i_up 
    for j=1:j_up
        for t=1:t_up
          P_sav(i,j,t)=P(i,j,t); %0608-initial P should be infinity, I use 2 for test
        end
    end    
end  


for i=1:i_up 
    for j=1:j_up
        for t=1:t_up
            U_sav(i,j,t)=U(i,j,t); 
        end
    end    
end


min_k_sav = min_k;


%1.2-add operational level constraints

C = [];

% c15
for i=1:i_up %c15
   for j=1:j_up
     for t=1:t_up
       if (t==1)
         C = [C, U(i,j,1) >= U_t0(i,j,1)+(b*d(i,j,1)+cijt*d(i,j,1)*m-n*T_abs*U(i,j,1))./(a+b+m*(cijt+P(i,j,t)))-X(i,j,1)]; %c15: t=1
       else       
         C = [C, U(i,j,t) >= U(i,j,t-1)+(b*d(i,j,t)+cijt*d(i,j,t)*m-n*T_abs*U(i,j,t))./(a+b+m*(cijt+P(i,j,t)))-X(i,j,t)]; %c15: t from 2
       end
     end
   end    
end

% 0424-c16-I suppose time window is M period for c16 - relax this one, hope get solution
% this also considers Xijt beyond T
for i=1:i_up %c16
    for j=1:j_up
      for t=1:t_up
        if (t==1)
          C = [C, sum(X(i,j,t:(t+M))) >= U_t0(i,j,1)+(b*d(i,j,t)+cijt*d(i,j,t)*m-n*T_abs*U(i,j,t))./(a+b+m*(cijt+P(i,j,t)))]; %c16: t=1
        else
          C = [C, sum(X(i,j,t:(t+M))) >= U(i,j,t-1)+(b*d(i,j,t)+cijt*d(i,j,t)*m-n*T_abs*U(i,j,t))./(a+b+m*(cijt+P(i,j,t)))]; %c16: t from 2
        end
      end  
    end     
end

% 0424 - c17 is in upper level, no need in lower level
% for i=1:i_up %c17
%   for j=1:j_up
%     for t=1:t_up
%       C = [C, ((b-1)*d(i,j,t)+cijt*d(i,j,t)*m)./(n*T_abs) <= U(i,j,t)];
%       C = [C, U(i,j,t) <= (b*d(i,j,t)+cijt*d(i,j,t)*m)./(n*T_abs)];
%     end
%   end
% end

% 0424-c18-here makes <= for Vi,t cuz alpha all set to 1, too relaxed, should be variaty
for i=1:i_up %c18
  for t=1:t_up
    if (t==1)
      C = [C, V(i,t) <= V_t0(i,1) + sum(alpha*X_t0(1:j_up,i,1)) - sum(X(i,1:j_up,t)) ]; %c18: t=1
    else
      C = [C, V(i,t) <= V(i,t-1) + sum(sum(alpha*X(1:j_up,i,1:(t-1)))) - sum(X(i,1:j_up,t)) ]; %c18: t from 2 (remember to use 2 sum for 2nd)
    end
  end
end

% 0508-c19-not added, cuz it looks already relaxed in math function, hidden property of our formula

for i=1:i_up %c20: Xijt
    for j=1:j_up
        for t=1:(t_up+M)
          C = [C, X(i,j,t) >= 0];
        end
    end    
end    

for i=1:i_up %c20: X_t0(i,j,1) is Xijt=0
    for j=1:j_up
        C = [C, X_t0(i,j,1) >= 0];
    end    
end    

for i=1:i_up %c20: Uijt
    for j=1:j_up
        for t=1:t_up
            C = [C, U(i,j,t) >= 0];
        end
    end    
end    

for i=1:i_up %c20: U_t0(i,j,1) is Uijt=0
    for j=1:j_up
        C = [C, U_t0(i,j,1) >= 0];
    end    
end 

for i=1:i_up %c20: Vit
    for t=1:t_up
      C = [C, V(i,t) >= 0];
    end    
end    

for i=1:i_up %c20: V_t0(i,1) is V(i,0)
  C = [C, V_t0(i,1) >= 0];
end

%1.3-solve operational level: settings: use cplex as solver

%0515-following is set optimality for MIP, but kang said make it LNContinuous-https://www.ibm.com/support/knowledgecenter/SSSA5P_12.9.0/ilog.odms.cplex.help/CPLEX/Parameters/topics/OptimalityTarget.html
%0515-'cplex.optimalitytarget'=1: Searches for a globally optimal solution to a convex model.
ops = sdpsettings ('solver','cplex','verbose',2,'cplex.optimalitytarget',1);%0508-verbose 2 shows iterations, 1 means shown normal message, 0 shows silent message
%0515-MIP limit tolerance gap to 1-https://www.ibm.com/support/knowledgecenter/SSSA5P_12.9.0/ilog.odms.cplex.help/CPLEX/Parameters/topics/EpGap.html
%ops = sdpsettings ('solver','cplex','verbose',2,'cplex.mip.tolerances.mipgap',1);%0508-verbose 2 shows iterations, 1 means shown normal message, 0 shows silent message

%0306-limit barrier iteration-https://www.ibm.com/support/knowledgecenter/SSSA5P_12.9.0/ilog.odms.cplex.help/CPLEX/Parameters/topics/BarItLim.html
%Cplex.Param.barrier.limits.iteration = 4;[not work]

%0306-limit network simplex iteration-https://www.ibm.com/support/knowledgecenter/SSSA5P_12.9.0/ilog.odms.cplex.help/CPLEX/Parameters/topics/NetItLim.html
%Cplex.Param.network.iterations = 2;[not work]

%ops.mip.limits.nodes = 4;[not work]

%0304-use the online method for out-of-memory: https://www.ilovematlab.cn/thread-266244-1-1.html
%ops.mip.strategy.file = 3; %[not work: should put in sdpsetting] this has same results 1159110 node, but reduces time for 3 min.

%0305-use online method for out-of-memory: https://www.ibm.com/developerworks/community/forums/html/topic?id=e3ea344c-7249-4edd-a47d-29e1f80fa480
%Cplex.Param.mip.limits.auxrootthreads = 2; %[not work: should put in sdpsetting] -1 or 1

%obj-operational obj z_op
z_op = sum(sum(sum(cijt.*X(1:i_up,1:j_up,1:t_up)))) + sum(sum(cijt.*X_t0(1:i_up,1:j_up,1)))...
  + sum(sum(sum(pijt.*U(1:i_up,1:j_up,1:t_up)))) + sum(sum(pijt.*U_t0(1:i_up,1:j_up,1)))...
  + sum(sum(hit.*V(1:i_up,1:t_up))) + sum(hit.*V_t0(1:i_up,1));

%1.35-express D_ij_SAV
for i=1:i_up
  for j=1:j_up
    for t=1:t_up
      D_ij_SAV(i,j,t) = (b*d(i,j,t)+cijt*d(i,j,t)*m-n*T_abs*U(i,j,t))./(a+b+m*(cijt+P(i,j,t)));
    end
  end
end

% 0604-calculate minimum k ratio covered by SAVs
for i=1:i_up
  for j=1:j_up
    for t=1:t_up
      k(i,j,t) = D_ij_SAV(i,j,t)./d(i,j,t);
    end
  end
end

% k_min = min(k(1,1,1),k(1,2,1),k(1,3,1))

% k_min = 0.99;
% % 
% if ((k_min) >= (k(1,1,1)))
%   k_min = k(1,1,1)                
% end  

% 0604-find minimum covered ratio k
% for i=1:i_up
%   for j=1:j_up
%     for t=1:t_up
%       if (k_min >= k(i,j,t))
%         k_min = k(i,j,t)                
%       end  
%     end
%   end
% end

%1.4-check if solved and output
result = optimize(C,z_op,ops); % 0508-origin
% result = solvesdp(C,z_op,ops); % 0508-this works

%result = optimize(C,z);
% if result.problem == 0 % problem =0 means solve successfully,only print Uijt here

%   fprintf('\nk value: ') %ratio k covered by SAVs
%   value(k)

  
%   fprintf('\nP_Ini_(1,1,1) = ')
%   value(P(1,1,1))

%   fprintf('\nMin cost by Lower = ')
%   value(z_op)

  all_k = [];

  for i=1:i_up
    for j=1:j_up
      for t=1:t_up
        all_k(end+1) = k(i,j,t);
      end
    end
  end
  
  min_k = min(all_k);
  
  k_avg = mean(all_k);
  
%   fprintf('\n min k: ')
%   value(min_k)
  
  max_k = max(all_k);
  
%   fprintf('\n max k: ')
%   value(max_k)
  
  % 0615-save max_k in txt

%   fid = fopen(['IniP=',num2str(2),'max_k.txt'],'a'); %save to last 'a'
%   fprintf(fid,'%.4f\t',double(max_k)); %0516-keep 4 decimals
%   fclose(fid);
  
%   fprintf('\nU(1,1,1) value: ')
%   value(U(1,1,1))
        
%   fprintf('\nX(1,1,1) value: ') 
%   value(X(1,1,1))

fprintf('\n Iteration No. = ')
value(Iter_count)



% fprintf('\n D_ij_SAV(1,1,1) value: ') %0515-as for now hide D_ij_SAV
% value(D_ij_SAV(1,1,1))
  
%   fprintf('\n k(1,1,1) ')
%   value(k(1,1,1))

%   fprintf('\n')
%     fprintf('\nV value: ')
%     value(V)
%     
%     fprintf('\nU_t0 value: ')
%     value(U_t0)
%     
%     fprintf('\nX_t0 value: ') 
%     value(X_t0)
%     
%     fprintf('\nV_t0 value: ')
%     value(V_t0)

% else
%   disp('Not solved!')
% 
% end

% 0519-save total cost z_op in txt
% fid = fopen(['Iter=',num2str(Iter_count),'_OperCost.txt'],'wt'); %creat & write z_op to txt
% fprintf(fid,'%.4f\n',double(z_op)); %0516-keep 4 decimals
% fclose(fid);

% 0605-save min ratio k covered by SAVs in txt
% fid_ratiok = fopen(['Iter=',num2str(Iter_count),'_MinRatiok.txt'],'wt'); %creat & write z_op to txt
% fprintf(fid_ratiok,'%.4f\n',double(min_k)); %0604-keep 4 decimals
% fclose(fid_ratiok);

% 0517-following is save Uijt in txt
%1.6-save Uijt as txt for upper use
% for U_t=1:t_up 
%      fid = fopen(['IniP=',num2str(IniP),'Iter=',num2str(Iter_count),'_Uijt=',num2str(U_t),'.txt'],'wt'); %creat & write Uijt to txt
%      for i=1:i_up
%          for j=1:j_up
%              if j == j_up
%                  fprintf(fid,'%.4f\n',double(U(i,j,U_t))); %0516-keep 4 decimals
%              else
%                  fprintf(fid,'%.4f\t',double(U(i,j,U_t))); %0516-keep 4 decimals
%              end     
%          end    
%      end    
%      fclose(fid);
% end

% 0607-calculate Upper Pijt upper bound
for i=1:i_up
  for j=1:j_up
    for t=1:t_up
      P(i,j,t)=(b*d(i,j,t)+cijt*d(i,j,t)*m-n*T_abs*U(i,j,t)-(a+b)*min_k*d(i,j,t))./(m*min_k*d(i,j,t))-cijt;
    end
  end
end

% fprintf('\nUpper bound P(1,1,1) = %.4f\n',double(P(1,1,1)));

%0608-----
%calculate Upper revenue
%0608-----

revenue=0;

for i=1:i_up
  for j=1:j_up
    for t=1:t_up
      if (t==1)||(t==2)||(t==9)||(t==10)||(t==17)||(t==18)||(t==31)||(t==32)||(t==39)||(t==40)
%         T_abs=2;
        revenue = revenue-(P(i,j,t).*(b.*d(i,j,t)+cijt.*d(i,j,t)*m-n*T_abs*U(i,j,t)))./(a+b+m.*(cijt+P(i,j,t))); 
      else
%         T_abs=2; %here should be 3, but temporarily set to 2, aligh with Oper level
        revenue = revenue-(P(i,j,t).*(b.*d(i,j,t)+cijt.*d(i,j,t)*m-n*T_abs*U(i,j,t)))./(a+b+m.*(cijt+P(i,j,t)));
      end
    end
  end
end

% 070320-avg. all P

 all_P = [];

  for i=1:i_up
    for j=1:j_up
      for t=1:t_up
        all_P = P(i,j,t);
      end
    end
  end

avg_P = all_P./(i_up*j_up*t_up);



% 070320-avg. all D

 all_Dijt = [];

  for i=1:i_up
    for j=1:j_up
      for t=1:t_up
        all_Dijt = D_ij_SAV(i,j,t);
      end
    end
  end

avg_Dijt = all_Dijt./(i_up*j_up*t_up);

% 070320-print results

fprintf('\n Upper bound P(1,1,1) = %.4f\n',double(P(1,1,1))); % particular P each iteration

fprintf('\n');

fprintf('\n D(1,1,1) = %.4f\n',double(D_ij_SAV(1,1,1))); % particular Dijt each iteration

fprintf('\n');

fprintf('\n Avg_P = %.4f\n',double(avg_P)); % particular avg_P each iteration

fprintf('\n');

fprintf('\n Avg_Dijt = %.4f\n',double(avg_Dijt)); % particular avg_Dijt each iteration

fprintf('\n');

fprintf('\n Max Revenue by Upper= %.4f\n',double(-revenue));

fprintf('\n');

fprintf('\n Min Cost by Lower= %.4f\n',double(z_op));

fprintf('\n');

fprintf('\n Profit = %.4f\n',double(-revenue-z_op));

fprintf('\n');

% test how many new P > old P

satis_P = 0;

for i=1:i_up
  for j=1:j_up
    for t=1:t_up
%       if (double(P_sav(i,j,t))<=double(P(i,j,t))) % 0610-check new P greater than old ones
        %0610-if new P no more than 1 cent than old, than add to satis_P terminate condition
      %if ((double(P_sav(i,j,t)))*1.0001 >= double(P(i,j,t))) 
%       if ((double(P_sav(i,j,t)))+0.01 >= double(P(i,j,t))) 
      if ((P_sav(i,j,t))+0.01 > P(i,j,t)) 
        satis_P = (satis_P)+1;
      end  
    end
  end
end

satis_P;

fprintf('\n');

% test how many new P all feasible that new_P <= new P upper bound condition


if (Iter_count==1)
  feasib_P = 0;
else
  
  feasib_P = 0;
  
  for i=1:i_up
    for j=1:j_up
      for t=1:t_up
        P_bench(i,j,t) = (b*d(i,j,t)+cijt*d(i,j,t)*m-n*T_abs*U(i,j,t)-(a+b)*min_k*d(i,j,t))./(m*min_k*d(i,j,t))-cijt;        
      end
    end
  end
  
    
  for i=1:i_up
    for j=1:j_up
      for t=1:t_up
  %       if (double(P_sav(i,j,t))<=double(P(i,j,t))) % 0610-check new P greater than old ones
          %0610-if new P no more than 1 cent than old, than add to satis_P terminate condition
        %if ((double(P_sav(i,j,t)))*1.0001 >= double(P(i,j,t))) 
        
%         if (double(P_sav(i,j,t)) <= double(P_bench(i,j,t)))
        if (P_bench(i,j,t) >= P_sav(i,j,t))
%           i,j,t
%           fprintf("\n");
          feasib_P = (feasib_P)+1;
        end  
      end
    end
  end

end

feasib_P;

fprintf('\n');

if (satis_P==i_up*j_up*t_up) && (feasib_P==i_up*j_up*t_up)
  consec_P = 1+consec_P;
else
  consec_P = 0;
end

% tik = tik-1;

end

fprintf('\n');
fprintf('\n');
fprintf('\n');


%-------case 5.3------

'-------Case 5.3-----'

IniP=100;

%initially price Pijt all the same
P=zeros(20,20,40);%create 20*20*40 matrix
for i=1:i_up 
    for j=1:j_up
        for t=1:t_up
          P(i,j,t)=IniP; %0608-initial P should be infinity, I use 2 for test
        end
    end    
end

% Pijt=2;

Iter_count = 0;

% 0517-M is c16 time window, set at 1 gives it relaxation
% 0424-M is c16 time window, set at 3 gives it relaxation
M=1; 

%0515-change all var to continours instead of MIP
X = sdpvar(i_up,j_up,t_up+M); %0424-Xi,j,t+M, cuz time window constraint c16 Xijt is beyond t
%X = intvar(i_up,j_up,t_up+M); %0424-Xi,j,t+M, cuz time window constraint c16 Xijt is beyond t
X_t0 = sdpvar(i_up,j_up,1); %0424-Xi,j,0 means moves initially for c18: inventory balance eq
%X_t0 = intvar(i_up,j_up,1); %0424-Xi,j,0 means moves initially for c18: inventory balance eq
U = sdpvar(i_up,j_up,t_up);
%U = intvar(i_up,j_up,t_up); 
U_t0 = sdpvar(i_up,j_up,1); %0424-this is for Ui,j,0 
%U_t0 = intvar(i_up,j_up,1); %0424-this is for Ui,j,0 
V = sdpvar(i_up,t_up); 
%V = intvar(i_up,t_up); 
V_t0 = sdpvar(i_up,1); %0424-this is especially for Vi,0
%V_t0 = intvar(i_up,1); %0424-this is especially for Vi,0
D_ij_SAV = sdpvar(i_up,j_up,t_up); %0425-this is to show D_ij_SAV, for futher plot, then find upper bound of Pijt
%D_ij_SAV = intvar(i_up,j_up,t_up); %0425-this is to show D_ij_SAV, for futher plot, then find upper bound of Pijt
% 0604-k ratio
k = sdpvar(i_up,j_up,t_up);

% 070220-set place for Wijt
W=zeros(20,20,40);%create 20*20*40 matrix

%1.1-parameters
% a=b=m=n=1,T=40
a=1;
b=1;
m=1;
n=1;
T=8; %070120-total 40 time periods, 5 time block, so every block is 8
T_abs=8; %070120-T_abs should be same as T
% cijt=40, pijt=10, hit=10, alphaijt=0.1 (all c,p,h,alpha are subjectively to be same)
cijt=40;
pijt=10; % 070120-increase for penalty of unmet demand
% 0425-intentionally make pijt this small to avoid all 0 in Uijt solution
hit=10;
alpha=1;

%%%-----0607----
%%% How to add loop?
%%%----

% 0610-here starts the loop

P_sav = zeros(20,20,40);
U_sav = zeros(20,20,40);
min_k = 0;
P_bench = zeros(20,20,40);

% tik = 2;

satis_P = 0;
feasib_P = 0;
consec_P = 0;

% while (satis_P<i_up*j_up*t_up) || (feasib_P<i_up*j_up*t_up)
% while Iter_count < 2
while (consec_P<3)

Iter_count = Iter_count+1; 


for i=1:i_up 
    for j=1:j_up
        for t=1:t_up
          P_sav(i,j,t)=P(i,j,t); %0608-initial P should be infinity, I use 2 for test
        end
    end    
end  


for i=1:i_up 
    for j=1:j_up
        for t=1:t_up
            U_sav(i,j,t)=U(i,j,t); 
        end
    end    
end


min_k_sav = min_k;


%1.2-add operational level constraints

C = [];

% c15
for i=1:i_up %c15
   for j=1:j_up
     for t=1:t_up
       if (t==1)
         C = [C, U(i,j,1) >= U_t0(i,j,1)+(b*d(i,j,1)+cijt*d(i,j,1)*m-n*T_abs*U(i,j,1))./(a+b+m*(cijt+P(i,j,t)))-X(i,j,1)]; %c15: t=1
       else       
         C = [C, U(i,j,t) >= U(i,j,t-1)+(b*d(i,j,t)+cijt*d(i,j,t)*m-n*T_abs*U(i,j,t))./(a+b+m*(cijt+P(i,j,t)))-X(i,j,t)]; %c15: t from 2
       end
     end
   end    
end

% 0424-c16-I suppose time window is M period for c16 - relax this one, hope get solution
% this also considers Xijt beyond T
for i=1:i_up %c16
    for j=1:j_up
      for t=1:t_up
        if (t==1)
          C = [C, sum(X(i,j,t:(t+M))) >= U_t0(i,j,1)+(b*d(i,j,t)+cijt*d(i,j,t)*m-n*T_abs*U(i,j,t))./(a+b+m*(cijt+P(i,j,t)))]; %c16: t=1
        else
          C = [C, sum(X(i,j,t:(t+M))) >= U(i,j,t-1)+(b*d(i,j,t)+cijt*d(i,j,t)*m-n*T_abs*U(i,j,t))./(a+b+m*(cijt+P(i,j,t)))]; %c16: t from 2
        end
      end  
    end     
end

% 0424 - c17 is in upper level, no need in lower level
% for i=1:i_up %c17
%   for j=1:j_up
%     for t=1:t_up
%       C = [C, ((b-1)*d(i,j,t)+cijt*d(i,j,t)*m)./(n*T_abs) <= U(i,j,t)];
%       C = [C, U(i,j,t) <= (b*d(i,j,t)+cijt*d(i,j,t)*m)./(n*T_abs)];
%     end
%   end
% end

% 0424-c18-here makes <= for Vi,t cuz alpha all set to 1, too relaxed, should be variaty
for i=1:i_up %c18
  for t=1:t_up
    if (t==1)
      C = [C, V(i,t) <= V_t0(i,1) + sum(alpha*X_t0(1:j_up,i,1)) - sum(X(i,1:j_up,t)) ]; %c18: t=1
    else
      C = [C, V(i,t) <= V(i,t-1) + sum(sum(alpha*X(1:j_up,i,1:(t-1)))) - sum(X(i,1:j_up,t)) ]; %c18: t from 2 (remember to use 2 sum for 2nd)
    end
  end
end

% 0508-c19-not added, cuz it looks already relaxed in math function, hidden property of our formula

for i=1:i_up %c20: Xijt
    for j=1:j_up
        for t=1:(t_up+M)
          C = [C, X(i,j,t) >= 0];
        end
    end    
end    

for i=1:i_up %c20: X_t0(i,j,1) is Xijt=0
    for j=1:j_up
        C = [C, X_t0(i,j,1) >= 0];
    end    
end    

for i=1:i_up %c20: Uijt
    for j=1:j_up
        for t=1:t_up
            C = [C, U(i,j,t) >= 0];
        end
    end    
end    

for i=1:i_up %c20: U_t0(i,j,1) is Uijt=0
    for j=1:j_up
        C = [C, U_t0(i,j,1) >= 0];
    end    
end 

for i=1:i_up %c20: Vit
    for t=1:t_up
      C = [C, V(i,t) >= 0];
    end    
end    

for i=1:i_up %c20: V_t0(i,1) is V(i,0)
  C = [C, V_t0(i,1) >= 0];
end

%1.3-solve operational level: settings: use cplex as solver

%0515-following is set optimality for MIP, but kang said make it LNContinuous-https://www.ibm.com/support/knowledgecenter/SSSA5P_12.9.0/ilog.odms.cplex.help/CPLEX/Parameters/topics/OptimalityTarget.html
%0515-'cplex.optimalitytarget'=1: Searches for a globally optimal solution to a convex model.
ops = sdpsettings ('solver','cplex','verbose',2,'cplex.optimalitytarget',1);%0508-verbose 2 shows iterations, 1 means shown normal message, 0 shows silent message
%0515-MIP limit tolerance gap to 1-https://www.ibm.com/support/knowledgecenter/SSSA5P_12.9.0/ilog.odms.cplex.help/CPLEX/Parameters/topics/EpGap.html
%ops = sdpsettings ('solver','cplex','verbose',2,'cplex.mip.tolerances.mipgap',1);%0508-verbose 2 shows iterations, 1 means shown normal message, 0 shows silent message

%0306-limit barrier iteration-https://www.ibm.com/support/knowledgecenter/SSSA5P_12.9.0/ilog.odms.cplex.help/CPLEX/Parameters/topics/BarItLim.html
%Cplex.Param.barrier.limits.iteration = 4;[not work]

%0306-limit network simplex iteration-https://www.ibm.com/support/knowledgecenter/SSSA5P_12.9.0/ilog.odms.cplex.help/CPLEX/Parameters/topics/NetItLim.html
%Cplex.Param.network.iterations = 2;[not work]

%ops.mip.limits.nodes = 4;[not work]

%0304-use the online method for out-of-memory: https://www.ilovematlab.cn/thread-266244-1-1.html
%ops.mip.strategy.file = 3; %[not work: should put in sdpsetting] this has same results 1159110 node, but reduces time for 3 min.

%0305-use online method for out-of-memory: https://www.ibm.com/developerworks/community/forums/html/topic?id=e3ea344c-7249-4edd-a47d-29e1f80fa480
%Cplex.Param.mip.limits.auxrootthreads = 2; %[not work: should put in sdpsetting] -1 or 1

%obj-operational obj z_op
z_op = sum(sum(sum(cijt.*X(1:i_up,1:j_up,1:t_up)))) + sum(sum(cijt.*X_t0(1:i_up,1:j_up,1)))...
  + sum(sum(sum(pijt.*U(1:i_up,1:j_up,1:t_up)))) + sum(sum(pijt.*U_t0(1:i_up,1:j_up,1)))...
  + sum(sum(hit.*V(1:i_up,1:t_up))) + sum(hit.*V_t0(1:i_up,1));

%1.35-express D_ij_SAV
for i=1:i_up
  for j=1:j_up
    for t=1:t_up
      D_ij_SAV(i,j,t) = (b*d(i,j,t)+cijt*d(i,j,t)*m-n*T_abs*U(i,j,t))./(a+b+m*(cijt+P(i,j,t)));
    end
  end
end

% 0604-calculate minimum k ratio covered by SAVs
for i=1:i_up
  for j=1:j_up
    for t=1:t_up
      k(i,j,t) = D_ij_SAV(i,j,t)./d(i,j,t);
    end
  end
end

% k_min = min(k(1,1,1),k(1,2,1),k(1,3,1))

% k_min = 0.99;
% % 
% if ((k_min) >= (k(1,1,1)))
%   k_min = k(1,1,1)                
% end  

% 0604-find minimum covered ratio k
% for i=1:i_up
%   for j=1:j_up
%     for t=1:t_up
%       if (k_min >= k(i,j,t))
%         k_min = k(i,j,t)                
%       end  
%     end
%   end
% end

%1.4-check if solved and output
result = optimize(C,z_op,ops); % 0508-origin
% result = solvesdp(C,z_op,ops); % 0508-this works

%result = optimize(C,z);
% if result.problem == 0 % problem =0 means solve successfully,only print Uijt here

%   fprintf('\nk value: ') %ratio k covered by SAVs
%   value(k)

  
%   fprintf('\nP_Ini_(1,1,1) = ')
%   value(P(1,1,1))

%   fprintf('\nMin cost by Lower = ')
%   value(z_op)

  all_k = [];

  for i=1:i_up
    for j=1:j_up
      for t=1:t_up
        all_k(end+1) = k(i,j,t);
      end
    end
  end
  
  min_k = min(all_k);
  
  k_avg = mean(all_k);
  
%   fprintf('\n min k: ')
%   value(min_k)
  
  max_k = max(all_k);
  
%   fprintf('\n max k: ')
%   value(max_k)
  
  % 0615-save max_k in txt

%   fid = fopen(['IniP=',num2str(2),'max_k.txt'],'a'); %save to last 'a'
%   fprintf(fid,'%.4f\t',double(max_k)); %0516-keep 4 decimals
%   fclose(fid);
  
%   fprintf('\nU(1,1,1) value: ')
%   value(U(1,1,1))
        
%   fprintf('\nX(1,1,1) value: ') 
%   value(X(1,1,1))

fprintf('\n Iteration No. = ')
value(Iter_count)



% fprintf('\n D_ij_SAV(1,1,1) value: ') %0515-as for now hide D_ij_SAV
% value(D_ij_SAV(1,1,1))
  
%   fprintf('\n k(1,1,1) ')
%   value(k(1,1,1))

%   fprintf('\n')
%     fprintf('\nV value: ')
%     value(V)
%     
%     fprintf('\nU_t0 value: ')
%     value(U_t0)
%     
%     fprintf('\nX_t0 value: ') 
%     value(X_t0)
%     
%     fprintf('\nV_t0 value: ')
%     value(V_t0)

% else
%   disp('Not solved!')
% 
% end

% 0519-save total cost z_op in txt
% fid = fopen(['Iter=',num2str(Iter_count),'_OperCost.txt'],'wt'); %creat & write z_op to txt
% fprintf(fid,'%.4f\n',double(z_op)); %0516-keep 4 decimals
% fclose(fid);

% 0605-save min ratio k covered by SAVs in txt
% fid_ratiok = fopen(['Iter=',num2str(Iter_count),'_MinRatiok.txt'],'wt'); %creat & write z_op to txt
% fprintf(fid_ratiok,'%.4f\n',double(min_k)); %0604-keep 4 decimals
% fclose(fid_ratiok);

% 0517-following is save Uijt in txt
%1.6-save Uijt as txt for upper use
% for U_t=1:t_up 
%      fid = fopen(['IniP=',num2str(IniP),'Iter=',num2str(Iter_count),'_Uijt=',num2str(U_t),'.txt'],'wt'); %creat & write Uijt to txt
%      for i=1:i_up
%          for j=1:j_up
%              if j == j_up
%                  fprintf(fid,'%.4f\n',double(U(i,j,U_t))); %0516-keep 4 decimals
%              else
%                  fprintf(fid,'%.4f\t',double(U(i,j,U_t))); %0516-keep 4 decimals
%              end     
%          end    
%      end    
%      fclose(fid);
% end

% 0607-calculate Upper Pijt upper bound
for i=1:i_up
  for j=1:j_up
    for t=1:t_up
      P(i,j,t)=(b*d(i,j,t)+cijt*d(i,j,t)*m-n*T_abs*U(i,j,t)-(a+b)*min_k*d(i,j,t))./(m*min_k*d(i,j,t))-cijt;
    end
  end
end

% fprintf('\nUpper bound P(1,1,1) = %.4f\n',double(P(1,1,1)));

%0608-----
%calculate Upper revenue
%0608-----

revenue=0;

for i=1:i_up
  for j=1:j_up
    for t=1:t_up
      if (t==1)||(t==2)||(t==9)||(t==10)||(t==17)||(t==18)||(t==31)||(t==32)||(t==39)||(t==40)
%         T_abs=2;
        revenue = revenue-(P(i,j,t).*(b.*d(i,j,t)+cijt.*d(i,j,t)*m-n*T_abs*U(i,j,t)))./(a+b+m.*(cijt+P(i,j,t))); 
      else
%         T_abs=2; %here should be 3, but temporarily set to 2, aligh with Oper level
        revenue = revenue-(P(i,j,t).*(b.*d(i,j,t)+cijt.*d(i,j,t)*m-n*T_abs*U(i,j,t)))./(a+b+m.*(cijt+P(i,j,t)));
      end
    end
  end
end

% 070320-avg. all P

 all_P = [];

  for i=1:i_up
    for j=1:j_up
      for t=1:t_up
        all_P = P(i,j,t);
      end
    end
  end

avg_P = all_P./(i_up*j_up*t_up);



% 070320-avg. all D

 all_Dijt = [];

  for i=1:i_up
    for j=1:j_up
      for t=1:t_up
        all_Dijt = D_ij_SAV(i,j,t);
      end
    end
  end

avg_Dijt = all_Dijt./(i_up*j_up*t_up);

% 070320-print results

fprintf('\n Upper bound P(1,1,1) = %.4f\n',double(P(1,1,1))); % particular P each iteration

fprintf('\n');

fprintf('\n D(1,1,1) = %.4f\n',double(D_ij_SAV(1,1,1))); % particular Dijt each iteration

fprintf('\n');

fprintf('\n Avg_P = %.4f\n',double(avg_P)); % particular avg_P each iteration

fprintf('\n');

fprintf('\n Avg_Dijt = %.4f\n',double(avg_Dijt)); % particular avg_Dijt each iteration

fprintf('\n');

fprintf('\n Max Revenue by Upper= %.4f\n',double(-revenue));

fprintf('\n');

fprintf('\n Min Cost by Lower= %.4f\n',double(z_op));

fprintf('\n');

fprintf('\n Profit = %.4f\n',double(-revenue-z_op));

fprintf('\n');

% test how many new P > old P

satis_P = 0;

for i=1:i_up
  for j=1:j_up
    for t=1:t_up
%       if (double(P_sav(i,j,t))<=double(P(i,j,t))) % 0610-check new P greater than old ones
        %0610-if new P no more than 1 cent than old, than add to satis_P terminate condition
      %if ((double(P_sav(i,j,t)))*1.0001 >= double(P(i,j,t))) 
%       if ((double(P_sav(i,j,t)))+0.01 >= double(P(i,j,t))) 
      if ((P_sav(i,j,t))+0.01 > P(i,j,t)) 
        satis_P = (satis_P)+1;
      end  
    end
  end
end

satis_P;

fprintf('\n');

% test how many new P all feasible that new_P <= new P upper bound condition


if (Iter_count==1)
  feasib_P = 0;
else
  
  feasib_P = 0;
  
  for i=1:i_up
    for j=1:j_up
      for t=1:t_up
        P_bench(i,j,t) = (b*d(i,j,t)+cijt*d(i,j,t)*m-n*T_abs*U(i,j,t)-(a+b)*min_k*d(i,j,t))./(m*min_k*d(i,j,t))-cijt;        
      end
    end
  end
  
    
  for i=1:i_up
    for j=1:j_up
      for t=1:t_up
  %       if (double(P_sav(i,j,t))<=double(P(i,j,t))) % 0610-check new P greater than old ones
          %0610-if new P no more than 1 cent than old, than add to satis_P terminate condition
        %if ((double(P_sav(i,j,t)))*1.0001 >= double(P(i,j,t))) 
        
%         if (double(P_sav(i,j,t)) <= double(P_bench(i,j,t)))
        if (P_bench(i,j,t) >= P_sav(i,j,t))
%           i,j,t
%           fprintf("\n");
          feasib_P = (feasib_P)+1;
        end  
      end
    end
  end

end

feasib_P;

fprintf('\n');

if (satis_P==i_up*j_up*t_up) && (feasib_P==i_up*j_up*t_up)
  consec_P = 1+consec_P;
else
  consec_P = 0;
end

% tik = tik-1;

end

fprintf('\n');
fprintf('\n');
fprintf('\n');



%-------case 6.1------

'------Case 6.1------'

IniP=1;

%initially price Pijt all the same
P=zeros(20,20,40);%create 20*20*40 matrix
for i=1:i_up 
    for j=1:j_up
        for t=1:t_up
          P(i,j,t)=IniP; %0608-initial P should be infinity, I use 2 for test
        end
    end    
end

% Pijt=2;

Iter_count = 0;

% 0517-M is c16 time window, set at 1 gives it relaxation
% 0424-M is c16 time window, set at 3 gives it relaxation
M=1; 

%0515-change all var to continours instead of MIP
X = sdpvar(i_up,j_up,t_up+M); %0424-Xi,j,t+M, cuz time window constraint c16 Xijt is beyond t
%X = intvar(i_up,j_up,t_up+M); %0424-Xi,j,t+M, cuz time window constraint c16 Xijt is beyond t
X_t0 = sdpvar(i_up,j_up,1); %0424-Xi,j,0 means moves initially for c18: inventory balance eq
%X_t0 = intvar(i_up,j_up,1); %0424-Xi,j,0 means moves initially for c18: inventory balance eq
U = sdpvar(i_up,j_up,t_up);
%U = intvar(i_up,j_up,t_up); 
U_t0 = sdpvar(i_up,j_up,1); %0424-this is for Ui,j,0 
%U_t0 = intvar(i_up,j_up,1); %0424-this is for Ui,j,0 
V = sdpvar(i_up,t_up); 
%V = intvar(i_up,t_up); 
V_t0 = sdpvar(i_up,1); %0424-this is especially for Vi,0
%V_t0 = intvar(i_up,1); %0424-this is especially for Vi,0
D_ij_SAV = sdpvar(i_up,j_up,t_up); %0425-this is to show D_ij_SAV, for futher plot, then find upper bound of Pijt
%D_ij_SAV = intvar(i_up,j_up,t_up); %0425-this is to show D_ij_SAV, for futher plot, then find upper bound of Pijt
% 0604-k ratio
k = sdpvar(i_up,j_up,t_up);

% 070220-set place for Wijt
W=zeros(20,20,40);%create 20*20*40 matrix

%1.1-parameters
% a=b=m=n=1,T=40
a=1;
b=1;
m=1;
n=1;
T=8; %070120-total 40 time periods, 5 time block, so every block is 8
T_abs=8; %070120-T_abs should be same as T
% cijt=40, pijt=10, hit=10, alphaijt=0.1 (all c,p,h,alpha are subjectively to be same)
cijt=20;
pijt=30; % 070120-increase for penalty of unmet demand
% 0425-intentionally make pijt this small to avoid all 0 in Uijt solution
hit=10;
alpha=1;

%%%-----0607----
%%% How to add loop?
%%%----

% 0610-here starts the loop

P_sav = zeros(20,20,40);
U_sav = zeros(20,20,40);
min_k = 0;
P_bench = zeros(20,20,40);

% tik = 2;

satis_P = 0;
feasib_P = 0;
consec_P = 0;

% while (satis_P<i_up*j_up*t_up) || (feasib_P<i_up*j_up*t_up)
% while Iter_count < 2
while (consec_P<3)

Iter_count = Iter_count+1; 


for i=1:i_up 
    for j=1:j_up
        for t=1:t_up
          P_sav(i,j,t)=P(i,j,t); %0608-initial P should be infinity, I use 2 for test
        end
    end    
end  


for i=1:i_up 
    for j=1:j_up
        for t=1:t_up
            U_sav(i,j,t)=U(i,j,t); 
        end
    end    
end


min_k_sav = min_k;


%1.2-add operational level constraints

C = [];

% c15
for i=1:i_up %c15
   for j=1:j_up
     for t=1:t_up
       if (t==1)
         C = [C, U(i,j,1) >= U_t0(i,j,1)+(b*d(i,j,1)+cijt*d(i,j,1)*m-n*T_abs*U(i,j,1))./(a+b+m*(cijt+P(i,j,t)))-X(i,j,1)]; %c15: t=1
       else       
         C = [C, U(i,j,t) >= U(i,j,t-1)+(b*d(i,j,t)+cijt*d(i,j,t)*m-n*T_abs*U(i,j,t))./(a+b+m*(cijt+P(i,j,t)))-X(i,j,t)]; %c15: t from 2
       end
     end
   end    
end

% 0424-c16-I suppose time window is M period for c16 - relax this one, hope get solution
% this also considers Xijt beyond T
for i=1:i_up %c16
    for j=1:j_up
      for t=1:t_up
        if (t==1)
          C = [C, sum(X(i,j,t:(t+M))) >= U_t0(i,j,1)+(b*d(i,j,t)+cijt*d(i,j,t)*m-n*T_abs*U(i,j,t))./(a+b+m*(cijt+P(i,j,t)))]; %c16: t=1
        else
          C = [C, sum(X(i,j,t:(t+M))) >= U(i,j,t-1)+(b*d(i,j,t)+cijt*d(i,j,t)*m-n*T_abs*U(i,j,t))./(a+b+m*(cijt+P(i,j,t)))]; %c16: t from 2
        end
      end  
    end     
end

% 0424 - c17 is in upper level, no need in lower level
% for i=1:i_up %c17
%   for j=1:j_up
%     for t=1:t_up
%       C = [C, ((b-1)*d(i,j,t)+cijt*d(i,j,t)*m)./(n*T_abs) <= U(i,j,t)];
%       C = [C, U(i,j,t) <= (b*d(i,j,t)+cijt*d(i,j,t)*m)./(n*T_abs)];
%     end
%   end
% end

% 0424-c18-here makes <= for Vi,t cuz alpha all set to 1, too relaxed, should be variaty
for i=1:i_up %c18
  for t=1:t_up
    if (t==1)
      C = [C, V(i,t) <= V_t0(i,1) + sum(alpha*X_t0(1:j_up,i,1)) - sum(X(i,1:j_up,t)) ]; %c18: t=1
    else
      C = [C, V(i,t) <= V(i,t-1) + sum(sum(alpha*X(1:j_up,i,1:(t-1)))) - sum(X(i,1:j_up,t)) ]; %c18: t from 2 (remember to use 2 sum for 2nd)
    end
  end
end

% 0508-c19-not added, cuz it looks already relaxed in math function, hidden property of our formula

for i=1:i_up %c20: Xijt
    for j=1:j_up
        for t=1:(t_up+M)
          C = [C, X(i,j,t) >= 0];
        end
    end    
end    

for i=1:i_up %c20: X_t0(i,j,1) is Xijt=0
    for j=1:j_up
        C = [C, X_t0(i,j,1) >= 0];
    end    
end    

for i=1:i_up %c20: Uijt
    for j=1:j_up
        for t=1:t_up
            C = [C, U(i,j,t) >= 0];
        end
    end    
end    

for i=1:i_up %c20: U_t0(i,j,1) is Uijt=0
    for j=1:j_up
        C = [C, U_t0(i,j,1) >= 0];
    end    
end 

for i=1:i_up %c20: Vit
    for t=1:t_up
      C = [C, V(i,t) >= 0];
    end    
end    

for i=1:i_up %c20: V_t0(i,1) is V(i,0)
  C = [C, V_t0(i,1) >= 0];
end

%1.3-solve operational level: settings: use cplex as solver

%0515-following is set optimality for MIP, but kang said make it LNContinuous-https://www.ibm.com/support/knowledgecenter/SSSA5P_12.9.0/ilog.odms.cplex.help/CPLEX/Parameters/topics/OptimalityTarget.html
%0515-'cplex.optimalitytarget'=1: Searches for a globally optimal solution to a convex model.
ops = sdpsettings ('solver','cplex','verbose',2,'cplex.optimalitytarget',1);%0508-verbose 2 shows iterations, 1 means shown normal message, 0 shows silent message
%0515-MIP limit tolerance gap to 1-https://www.ibm.com/support/knowledgecenter/SSSA5P_12.9.0/ilog.odms.cplex.help/CPLEX/Parameters/topics/EpGap.html
%ops = sdpsettings ('solver','cplex','verbose',2,'cplex.mip.tolerances.mipgap',1);%0508-verbose 2 shows iterations, 1 means shown normal message, 0 shows silent message

%0306-limit barrier iteration-https://www.ibm.com/support/knowledgecenter/SSSA5P_12.9.0/ilog.odms.cplex.help/CPLEX/Parameters/topics/BarItLim.html
%Cplex.Param.barrier.limits.iteration = 4;[not work]

%0306-limit network simplex iteration-https://www.ibm.com/support/knowledgecenter/SSSA5P_12.9.0/ilog.odms.cplex.help/CPLEX/Parameters/topics/NetItLim.html
%Cplex.Param.network.iterations = 2;[not work]

%ops.mip.limits.nodes = 4;[not work]

%0304-use the online method for out-of-memory: https://www.ilovematlab.cn/thread-266244-1-1.html
%ops.mip.strategy.file = 3; %[not work: should put in sdpsetting] this has same results 1159110 node, but reduces time for 3 min.

%0305-use online method for out-of-memory: https://www.ibm.com/developerworks/community/forums/html/topic?id=e3ea344c-7249-4edd-a47d-29e1f80fa480
%Cplex.Param.mip.limits.auxrootthreads = 2; %[not work: should put in sdpsetting] -1 or 1

%obj-operational obj z_op
z_op = sum(sum(sum(cijt.*X(1:i_up,1:j_up,1:t_up)))) + sum(sum(cijt.*X_t0(1:i_up,1:j_up,1)))...
  + sum(sum(sum(pijt.*U(1:i_up,1:j_up,1:t_up)))) + sum(sum(pijt.*U_t0(1:i_up,1:j_up,1)))...
  + sum(sum(hit.*V(1:i_up,1:t_up))) + sum(hit.*V_t0(1:i_up,1));

%1.35-express D_ij_SAV
for i=1:i_up
  for j=1:j_up
    for t=1:t_up
      D_ij_SAV(i,j,t) = (b*d(i,j,t)+cijt*d(i,j,t)*m-n*T_abs*U(i,j,t))./(a+b+m*(cijt+P(i,j,t)));
    end
  end
end

% 0604-calculate minimum k ratio covered by SAVs
for i=1:i_up
  for j=1:j_up
    for t=1:t_up
      k(i,j,t) = D_ij_SAV(i,j,t)./d(i,j,t);
    end
  end
end

% k_min = min(k(1,1,1),k(1,2,1),k(1,3,1))

% k_min = 0.99;
% % 
% if ((k_min) >= (k(1,1,1)))
%   k_min = k(1,1,1)                
% end  

% 0604-find minimum covered ratio k
% for i=1:i_up
%   for j=1:j_up
%     for t=1:t_up
%       if (k_min >= k(i,j,t))
%         k_min = k(i,j,t)                
%       end  
%     end
%   end
% end

%1.4-check if solved and output
result = optimize(C,z_op,ops); % 0508-origin
% result = solvesdp(C,z_op,ops); % 0508-this works

%result = optimize(C,z);
% if result.problem == 0 % problem =0 means solve successfully,only print Uijt here

%   fprintf('\nk value: ') %ratio k covered by SAVs
%   value(k)

  
%   fprintf('\nP_Ini_(1,1,1) = ')
%   value(P(1,1,1))

%   fprintf('\nMin cost by Lower = ')
%   value(z_op)

  all_k = [];

  for i=1:i_up
    for j=1:j_up
      for t=1:t_up
        all_k(end+1) = k(i,j,t);
      end
    end
  end
  
  min_k = min(all_k);
  
  k_avg = mean(all_k);
  
%   fprintf('\n min k: ')
%   value(min_k)
  
  max_k = max(all_k);
  
%   fprintf('\n max k: ')
%   value(max_k)
  
  % 0615-save max_k in txt

%   fid = fopen(['IniP=',num2str(2),'max_k.txt'],'a'); %save to last 'a'
%   fprintf(fid,'%.4f\t',double(max_k)); %0516-keep 4 decimals
%   fclose(fid);
  
%   fprintf('\nU(1,1,1) value: ')
%   value(U(1,1,1))
        
%   fprintf('\nX(1,1,1) value: ') 
%   value(X(1,1,1))

fprintf('\n Iteration No. = ')
value(Iter_count)



% fprintf('\n D_ij_SAV(1,1,1) value: ') %0515-as for now hide D_ij_SAV
% value(D_ij_SAV(1,1,1))
  
%   fprintf('\n k(1,1,1) ')
%   value(k(1,1,1))

%   fprintf('\n')
%     fprintf('\nV value: ')
%     value(V)
%     
%     fprintf('\nU_t0 value: ')
%     value(U_t0)
%     
%     fprintf('\nX_t0 value: ') 
%     value(X_t0)
%     
%     fprintf('\nV_t0 value: ')
%     value(V_t0)

% else
%   disp('Not solved!')
% 
% end

% 0519-save total cost z_op in txt
% fid = fopen(['Iter=',num2str(Iter_count),'_OperCost.txt'],'wt'); %creat & write z_op to txt
% fprintf(fid,'%.4f\n',double(z_op)); %0516-keep 4 decimals
% fclose(fid);

% 0605-save min ratio k covered by SAVs in txt
% fid_ratiok = fopen(['Iter=',num2str(Iter_count),'_MinRatiok.txt'],'wt'); %creat & write z_op to txt
% fprintf(fid_ratiok,'%.4f\n',double(min_k)); %0604-keep 4 decimals
% fclose(fid_ratiok);

% 0517-following is save Uijt in txt
%1.6-save Uijt as txt for upper use
% for U_t=1:t_up 
%      fid = fopen(['IniP=',num2str(IniP),'Iter=',num2str(Iter_count),'_Uijt=',num2str(U_t),'.txt'],'wt'); %creat & write Uijt to txt
%      for i=1:i_up
%          for j=1:j_up
%              if j == j_up
%                  fprintf(fid,'%.4f\n',double(U(i,j,U_t))); %0516-keep 4 decimals
%              else
%                  fprintf(fid,'%.4f\t',double(U(i,j,U_t))); %0516-keep 4 decimals
%              end     
%          end    
%      end    
%      fclose(fid);
% end

% 0607-calculate Upper Pijt upper bound
for i=1:i_up
  for j=1:j_up
    for t=1:t_up
      P(i,j,t)=(b*d(i,j,t)+cijt*d(i,j,t)*m-n*T_abs*U(i,j,t)-(a+b)*min_k*d(i,j,t))./(m*min_k*d(i,j,t))-cijt;
    end
  end
end

% fprintf('\nUpper bound P(1,1,1) = %.4f\n',double(P(1,1,1)));

%0608-----
%calculate Upper revenue
%0608-----

revenue=0;

for i=1:i_up
  for j=1:j_up
    for t=1:t_up
      if (t==1)||(t==2)||(t==9)||(t==10)||(t==17)||(t==18)||(t==31)||(t==32)||(t==39)||(t==40)
%         T_abs=2;
        revenue = revenue-(P(i,j,t).*(b.*d(i,j,t)+cijt.*d(i,j,t)*m-n*T_abs*U(i,j,t)))./(a+b+m.*(cijt+P(i,j,t))); 
      else
%         T_abs=2; %here should be 3, but temporarily set to 2, aligh with Oper level
        revenue = revenue-(P(i,j,t).*(b.*d(i,j,t)+cijt.*d(i,j,t)*m-n*T_abs*U(i,j,t)))./(a+b+m.*(cijt+P(i,j,t)));
      end
    end
  end
end

% 070320-avg. all P

 all_P = [];

  for i=1:i_up
    for j=1:j_up
      for t=1:t_up
        all_P = P(i,j,t);
      end
    end
  end

avg_P = all_P./(i_up*j_up*t_up);



% 070320-avg. all D

 all_Dijt = [];

  for i=1:i_up
    for j=1:j_up
      for t=1:t_up
        all_Dijt = D_ij_SAV(i,j,t);
      end
    end
  end

avg_Dijt = all_Dijt./(i_up*j_up*t_up);

% 070320-print results

fprintf('\n Upper bound P(1,1,1) = %.4f\n',double(P(1,1,1))); % particular P each iteration

fprintf('\n');

fprintf('\n D(1,1,1) = %.4f\n',double(D_ij_SAV(1,1,1))); % particular Dijt each iteration

fprintf('\n');

fprintf('\n Avg_P = %.4f\n',double(avg_P)); % particular avg_P each iteration

fprintf('\n');

fprintf('\n Avg_Dijt = %.4f\n',double(avg_Dijt)); % particular avg_Dijt each iteration

fprintf('\n');

fprintf('\n Max Revenue by Upper= %.4f\n',double(-revenue));

fprintf('\n');

fprintf('\n Min Cost by Lower= %.4f\n',double(z_op));

fprintf('\n');

fprintf('\n Profit = %.4f\n',double(-revenue-z_op));

fprintf('\n');

% test how many new P > old P

satis_P = 0;

for i=1:i_up
  for j=1:j_up
    for t=1:t_up
%       if (double(P_sav(i,j,t))<=double(P(i,j,t))) % 0610-check new P greater than old ones
        %0610-if new P no more than 1 cent than old, than add to satis_P terminate condition
      %if ((double(P_sav(i,j,t)))*1.0001 >= double(P(i,j,t))) 
%       if ((double(P_sav(i,j,t)))+0.01 >= double(P(i,j,t))) 
      if ((P_sav(i,j,t))+0.01 > P(i,j,t)) 
        satis_P = (satis_P)+1;
      end  
    end
  end
end

satis_P;

fprintf('\n');

% test how many new P all feasible that new_P <= new P upper bound condition


if (Iter_count==1)
  feasib_P = 0;
else
  
  feasib_P = 0;
  
  for i=1:i_up
    for j=1:j_up
      for t=1:t_up
        P_bench(i,j,t) = (b*d(i,j,t)+cijt*d(i,j,t)*m-n*T_abs*U(i,j,t)-(a+b)*min_k*d(i,j,t))./(m*min_k*d(i,j,t))-cijt;        
      end
    end
  end
  
    
  for i=1:i_up
    for j=1:j_up
      for t=1:t_up
  %       if (double(P_sav(i,j,t))<=double(P(i,j,t))) % 0610-check new P greater than old ones
          %0610-if new P no more than 1 cent than old, than add to satis_P terminate condition
        %if ((double(P_sav(i,j,t)))*1.0001 >= double(P(i,j,t))) 
        
%         if (double(P_sav(i,j,t)) <= double(P_bench(i,j,t)))
        if (P_bench(i,j,t) >= P_sav(i,j,t))
%           i,j,t
%           fprintf("\n");
          feasib_P = (feasib_P)+1;
        end  
      end
    end
  end

end

feasib_P;

fprintf('\n');

if (satis_P==i_up*j_up*t_up) && (feasib_P==i_up*j_up*t_up)
  consec_P = 1+consec_P;
else
  consec_P = 0;
end

% tik = tik-1;

end

fprintf('\n');
fprintf('\n');
fprintf('\n');


%-------case 6.2------

'-------Case 6.2-------'

IniP=10;

%initially price Pijt all the same
P=zeros(20,20,40);%create 20*20*40 matrix
for i=1:i_up 
    for j=1:j_up
        for t=1:t_up
          P(i,j,t)=IniP; %0608-initial P should be infinity, I use 2 for test
        end
    end    
end

% Pijt=2;

Iter_count = 0;

% 0517-M is c16 time window, set at 1 gives it relaxation
% 0424-M is c16 time window, set at 3 gives it relaxation
M=1; 

%0515-change all var to continours instead of MIP
X = sdpvar(i_up,j_up,t_up+M); %0424-Xi,j,t+M, cuz time window constraint c16 Xijt is beyond t
%X = intvar(i_up,j_up,t_up+M); %0424-Xi,j,t+M, cuz time window constraint c16 Xijt is beyond t
X_t0 = sdpvar(i_up,j_up,1); %0424-Xi,j,0 means moves initially for c18: inventory balance eq
%X_t0 = intvar(i_up,j_up,1); %0424-Xi,j,0 means moves initially for c18: inventory balance eq
U = sdpvar(i_up,j_up,t_up);
%U = intvar(i_up,j_up,t_up); 
U_t0 = sdpvar(i_up,j_up,1); %0424-this is for Ui,j,0 
%U_t0 = intvar(i_up,j_up,1); %0424-this is for Ui,j,0 
V = sdpvar(i_up,t_up); 
%V = intvar(i_up,t_up); 
V_t0 = sdpvar(i_up,1); %0424-this is especially for Vi,0
%V_t0 = intvar(i_up,1); %0424-this is especially for Vi,0
D_ij_SAV = sdpvar(i_up,j_up,t_up); %0425-this is to show D_ij_SAV, for futher plot, then find upper bound of Pijt
%D_ij_SAV = intvar(i_up,j_up,t_up); %0425-this is to show D_ij_SAV, for futher plot, then find upper bound of Pijt
% 0604-k ratio
k = sdpvar(i_up,j_up,t_up);

% 070220-set place for Wijt
W=zeros(20,20,40);%create 20*20*40 matrix

%1.1-parameters
% a=b=m=n=1,T=40
a=1;
b=1;
m=1;
n=1;
T=8; %070120-total 40 time periods, 5 time block, so every block is 8
T_abs=8; %070120-T_abs should be same as T
% cijt=40, pijt=10, hit=10, alphaijt=0.1 (all c,p,h,alpha are subjectively to be same)
cijt=20;
pijt=30; % 070120-increase for penalty of unmet demand
% 0425-intentionally make pijt this small to avoid all 0 in Uijt solution
hit=10;
alpha=1;

%%%-----0607----
%%% How to add loop?
%%%----

% 0610-here starts the loop

P_sav = zeros(20,20,40);
U_sav = zeros(20,20,40);
min_k = 0;
P_bench = zeros(20,20,40);

% tik = 2;

satis_P = 0;
feasib_P = 0;
consec_P = 0;

% while (satis_P<i_up*j_up*t_up) || (feasib_P<i_up*j_up*t_up)
% while Iter_count < 2
while (consec_P<3)

Iter_count = Iter_count+1; 


for i=1:i_up 
    for j=1:j_up
        for t=1:t_up
          P_sav(i,j,t)=P(i,j,t); %0608-initial P should be infinity, I use 2 for test
        end
    end    
end  


for i=1:i_up 
    for j=1:j_up
        for t=1:t_up
            U_sav(i,j,t)=U(i,j,t); 
        end
    end    
end


min_k_sav = min_k;


%1.2-add operational level constraints

C = [];

% c15
for i=1:i_up %c15
   for j=1:j_up
     for t=1:t_up
       if (t==1)
         C = [C, U(i,j,1) >= U_t0(i,j,1)+(b*d(i,j,1)+cijt*d(i,j,1)*m-n*T_abs*U(i,j,1))./(a+b+m*(cijt+P(i,j,t)))-X(i,j,1)]; %c15: t=1
       else       
         C = [C, U(i,j,t) >= U(i,j,t-1)+(b*d(i,j,t)+cijt*d(i,j,t)*m-n*T_abs*U(i,j,t))./(a+b+m*(cijt+P(i,j,t)))-X(i,j,t)]; %c15: t from 2
       end
     end
   end    
end

% 0424-c16-I suppose time window is M period for c16 - relax this one, hope get solution
% this also considers Xijt beyond T
for i=1:i_up %c16
    for j=1:j_up
      for t=1:t_up
        if (t==1)
          C = [C, sum(X(i,j,t:(t+M))) >= U_t0(i,j,1)+(b*d(i,j,t)+cijt*d(i,j,t)*m-n*T_abs*U(i,j,t))./(a+b+m*(cijt+P(i,j,t)))]; %c16: t=1
        else
          C = [C, sum(X(i,j,t:(t+M))) >= U(i,j,t-1)+(b*d(i,j,t)+cijt*d(i,j,t)*m-n*T_abs*U(i,j,t))./(a+b+m*(cijt+P(i,j,t)))]; %c16: t from 2
        end
      end  
    end     
end

% 0424 - c17 is in upper level, no need in lower level
% for i=1:i_up %c17
%   for j=1:j_up
%     for t=1:t_up
%       C = [C, ((b-1)*d(i,j,t)+cijt*d(i,j,t)*m)./(n*T_abs) <= U(i,j,t)];
%       C = [C, U(i,j,t) <= (b*d(i,j,t)+cijt*d(i,j,t)*m)./(n*T_abs)];
%     end
%   end
% end

% 0424-c18-here makes <= for Vi,t cuz alpha all set to 1, too relaxed, should be variaty
for i=1:i_up %c18
  for t=1:t_up
    if (t==1)
      C = [C, V(i,t) <= V_t0(i,1) + sum(alpha*X_t0(1:j_up,i,1)) - sum(X(i,1:j_up,t)) ]; %c18: t=1
    else
      C = [C, V(i,t) <= V(i,t-1) + sum(sum(alpha*X(1:j_up,i,1:(t-1)))) - sum(X(i,1:j_up,t)) ]; %c18: t from 2 (remember to use 2 sum for 2nd)
    end
  end
end

% 0508-c19-not added, cuz it looks already relaxed in math function, hidden property of our formula

for i=1:i_up %c20: Xijt
    for j=1:j_up
        for t=1:(t_up+M)
          C = [C, X(i,j,t) >= 0];
        end
    end    
end    

for i=1:i_up %c20: X_t0(i,j,1) is Xijt=0
    for j=1:j_up
        C = [C, X_t0(i,j,1) >= 0];
    end    
end    

for i=1:i_up %c20: Uijt
    for j=1:j_up
        for t=1:t_up
            C = [C, U(i,j,t) >= 0];
        end
    end    
end    

for i=1:i_up %c20: U_t0(i,j,1) is Uijt=0
    for j=1:j_up
        C = [C, U_t0(i,j,1) >= 0];
    end    
end 

for i=1:i_up %c20: Vit
    for t=1:t_up
      C = [C, V(i,t) >= 0];
    end    
end    

for i=1:i_up %c20: V_t0(i,1) is V(i,0)
  C = [C, V_t0(i,1) >= 0];
end

%1.3-solve operational level: settings: use cplex as solver

%0515-following is set optimality for MIP, but kang said make it LNContinuous-https://www.ibm.com/support/knowledgecenter/SSSA5P_12.9.0/ilog.odms.cplex.help/CPLEX/Parameters/topics/OptimalityTarget.html
%0515-'cplex.optimalitytarget'=1: Searches for a globally optimal solution to a convex model.
ops = sdpsettings ('solver','cplex','verbose',2,'cplex.optimalitytarget',1);%0508-verbose 2 shows iterations, 1 means shown normal message, 0 shows silent message
%0515-MIP limit tolerance gap to 1-https://www.ibm.com/support/knowledgecenter/SSSA5P_12.9.0/ilog.odms.cplex.help/CPLEX/Parameters/topics/EpGap.html
%ops = sdpsettings ('solver','cplex','verbose',2,'cplex.mip.tolerances.mipgap',1);%0508-verbose 2 shows iterations, 1 means shown normal message, 0 shows silent message

%0306-limit barrier iteration-https://www.ibm.com/support/knowledgecenter/SSSA5P_12.9.0/ilog.odms.cplex.help/CPLEX/Parameters/topics/BarItLim.html
%Cplex.Param.barrier.limits.iteration = 4;[not work]

%0306-limit network simplex iteration-https://www.ibm.com/support/knowledgecenter/SSSA5P_12.9.0/ilog.odms.cplex.help/CPLEX/Parameters/topics/NetItLim.html
%Cplex.Param.network.iterations = 2;[not work]

%ops.mip.limits.nodes = 4;[not work]

%0304-use the online method for out-of-memory: https://www.ilovematlab.cn/thread-266244-1-1.html
%ops.mip.strategy.file = 3; %[not work: should put in sdpsetting] this has same results 1159110 node, but reduces time for 3 min.

%0305-use online method for out-of-memory: https://www.ibm.com/developerworks/community/forums/html/topic?id=e3ea344c-7249-4edd-a47d-29e1f80fa480
%Cplex.Param.mip.limits.auxrootthreads = 2; %[not work: should put in sdpsetting] -1 or 1

%obj-operational obj z_op
z_op = sum(sum(sum(cijt.*X(1:i_up,1:j_up,1:t_up)))) + sum(sum(cijt.*X_t0(1:i_up,1:j_up,1)))...
  + sum(sum(sum(pijt.*U(1:i_up,1:j_up,1:t_up)))) + sum(sum(pijt.*U_t0(1:i_up,1:j_up,1)))...
  + sum(sum(hit.*V(1:i_up,1:t_up))) + sum(hit.*V_t0(1:i_up,1));

%1.35-express D_ij_SAV
for i=1:i_up
  for j=1:j_up
    for t=1:t_up
      D_ij_SAV(i,j,t) = (b*d(i,j,t)+cijt*d(i,j,t)*m-n*T_abs*U(i,j,t))./(a+b+m*(cijt+P(i,j,t)));
    end
  end
end

% 0604-calculate minimum k ratio covered by SAVs
for i=1:i_up
  for j=1:j_up
    for t=1:t_up
      k(i,j,t) = D_ij_SAV(i,j,t)./d(i,j,t);
    end
  end
end

% k_min = min(k(1,1,1),k(1,2,1),k(1,3,1))

% k_min = 0.99;
% % 
% if ((k_min) >= (k(1,1,1)))
%   k_min = k(1,1,1)                
% end  

% 0604-find minimum covered ratio k
% for i=1:i_up
%   for j=1:j_up
%     for t=1:t_up
%       if (k_min >= k(i,j,t))
%         k_min = k(i,j,t)                
%       end  
%     end
%   end
% end

%1.4-check if solved and output
result = optimize(C,z_op,ops); % 0508-origin
% result = solvesdp(C,z_op,ops); % 0508-this works

%result = optimize(C,z);
% if result.problem == 0 % problem =0 means solve successfully,only print Uijt here

%   fprintf('\nk value: ') %ratio k covered by SAVs
%   value(k)

  
%   fprintf('\nP_Ini_(1,1,1) = ')
%   value(P(1,1,1))

%   fprintf('\nMin cost by Lower = ')
%   value(z_op)

  all_k = [];

  for i=1:i_up
    for j=1:j_up
      for t=1:t_up
        all_k(end+1) = k(i,j,t);
      end
    end
  end
  
  min_k = min(all_k);
  
  k_avg = mean(all_k);
  
%   fprintf('\n min k: ')
%   value(min_k)
  
  max_k = max(all_k);
  
%   fprintf('\n max k: ')
%   value(max_k)
  
  % 0615-save max_k in txt

%   fid = fopen(['IniP=',num2str(2),'max_k.txt'],'a'); %save to last 'a'
%   fprintf(fid,'%.4f\t',double(max_k)); %0516-keep 4 decimals
%   fclose(fid);
  
%   fprintf('\nU(1,1,1) value: ')
%   value(U(1,1,1))
        
%   fprintf('\nX(1,1,1) value: ') 
%   value(X(1,1,1))

fprintf('\n Iteration No. = ')
value(Iter_count)



% fprintf('\n D_ij_SAV(1,1,1) value: ') %0515-as for now hide D_ij_SAV
% value(D_ij_SAV(1,1,1))
  
%   fprintf('\n k(1,1,1) ')
%   value(k(1,1,1))

%   fprintf('\n')
%     fprintf('\nV value: ')
%     value(V)
%     
%     fprintf('\nU_t0 value: ')
%     value(U_t0)
%     
%     fprintf('\nX_t0 value: ') 
%     value(X_t0)
%     
%     fprintf('\nV_t0 value: ')
%     value(V_t0)

% else
%   disp('Not solved!')
% 
% end

% 0519-save total cost z_op in txt
% fid = fopen(['Iter=',num2str(Iter_count),'_OperCost.txt'],'wt'); %creat & write z_op to txt
% fprintf(fid,'%.4f\n',double(z_op)); %0516-keep 4 decimals
% fclose(fid);

% 0605-save min ratio k covered by SAVs in txt
% fid_ratiok = fopen(['Iter=',num2str(Iter_count),'_MinRatiok.txt'],'wt'); %creat & write z_op to txt
% fprintf(fid_ratiok,'%.4f\n',double(min_k)); %0604-keep 4 decimals
% fclose(fid_ratiok);

% 0517-following is save Uijt in txt
%1.6-save Uijt as txt for upper use
% for U_t=1:t_up 
%      fid = fopen(['IniP=',num2str(IniP),'Iter=',num2str(Iter_count),'_Uijt=',num2str(U_t),'.txt'],'wt'); %creat & write Uijt to txt
%      for i=1:i_up
%          for j=1:j_up
%              if j == j_up
%                  fprintf(fid,'%.4f\n',double(U(i,j,U_t))); %0516-keep 4 decimals
%              else
%                  fprintf(fid,'%.4f\t',double(U(i,j,U_t))); %0516-keep 4 decimals
%              end     
%          end    
%      end    
%      fclose(fid);
% end

% 0607-calculate Upper Pijt upper bound
for i=1:i_up
  for j=1:j_up
    for t=1:t_up
      P(i,j,t)=(b*d(i,j,t)+cijt*d(i,j,t)*m-n*T_abs*U(i,j,t)-(a+b)*min_k*d(i,j,t))./(m*min_k*d(i,j,t))-cijt;
    end
  end
end

% fprintf('\nUpper bound P(1,1,1) = %.4f\n',double(P(1,1,1)));

%0608-----
%calculate Upper revenue
%0608-----

revenue=0;

for i=1:i_up
  for j=1:j_up
    for t=1:t_up
      if (t==1)||(t==2)||(t==9)||(t==10)||(t==17)||(t==18)||(t==31)||(t==32)||(t==39)||(t==40)
%         T_abs=2;
        revenue = revenue-(P(i,j,t).*(b.*d(i,j,t)+cijt.*d(i,j,t)*m-n*T_abs*U(i,j,t)))./(a+b+m.*(cijt+P(i,j,t))); 
      else
%         T_abs=2; %here should be 3, but temporarily set to 2, aligh with Oper level
        revenue = revenue-(P(i,j,t).*(b.*d(i,j,t)+cijt.*d(i,j,t)*m-n*T_abs*U(i,j,t)))./(a+b+m.*(cijt+P(i,j,t)));
      end
    end
  end
end

% 070320-avg. all P

 all_P = [];

  for i=1:i_up
    for j=1:j_up
      for t=1:t_up
        all_P = P(i,j,t);
      end
    end
  end

avg_P = all_P./(i_up*j_up*t_up);



% 070320-avg. all D

 all_Dijt = [];

  for i=1:i_up
    for j=1:j_up
      for t=1:t_up
        all_Dijt = D_ij_SAV(i,j,t);
      end
    end
  end

avg_Dijt = all_Dijt./(i_up*j_up*t_up);

% 070320-print results

fprintf('\n Upper bound P(1,1,1) = %.4f\n',double(P(1,1,1))); % particular P each iteration

fprintf('\n');

fprintf('\n D(1,1,1) = %.4f\n',double(D_ij_SAV(1,1,1))); % particular Dijt each iteration

fprintf('\n');

fprintf('\n Avg_P = %.4f\n',double(avg_P)); % particular avg_P each iteration

fprintf('\n');

fprintf('\n Avg_Dijt = %.4f\n',double(avg_Dijt)); % particular avg_Dijt each iteration

fprintf('\n');

fprintf('\n Max Revenue by Upper= %.4f\n',double(-revenue));

fprintf('\n');

fprintf('\n Min Cost by Lower= %.4f\n',double(z_op));

fprintf('\n');

fprintf('\n Profit = %.4f\n',double(-revenue-z_op));

fprintf('\n');

% test how many new P > old P

satis_P = 0;

for i=1:i_up
  for j=1:j_up
    for t=1:t_up
%       if (double(P_sav(i,j,t))<=double(P(i,j,t))) % 0610-check new P greater than old ones
        %0610-if new P no more than 1 cent than old, than add to satis_P terminate condition
      %if ((double(P_sav(i,j,t)))*1.0001 >= double(P(i,j,t))) 
%       if ((double(P_sav(i,j,t)))+0.01 >= double(P(i,j,t))) 
      if ((P_sav(i,j,t))+0.01 > P(i,j,t)) 
        satis_P = (satis_P)+1;
      end  
    end
  end
end

satis_P;

fprintf('\n');

% test how many new P all feasible that new_P <= new P upper bound condition


if (Iter_count==1)
  feasib_P = 0;
else
  
  feasib_P = 0;
  
  for i=1:i_up
    for j=1:j_up
      for t=1:t_up
        P_bench(i,j,t) = (b*d(i,j,t)+cijt*d(i,j,t)*m-n*T_abs*U(i,j,t)-(a+b)*min_k*d(i,j,t))./(m*min_k*d(i,j,t))-cijt;        
      end
    end
  end
  
    
  for i=1:i_up
    for j=1:j_up
      for t=1:t_up
  %       if (double(P_sav(i,j,t))<=double(P(i,j,t))) % 0610-check new P greater than old ones
          %0610-if new P no more than 1 cent than old, than add to satis_P terminate condition
        %if ((double(P_sav(i,j,t)))*1.0001 >= double(P(i,j,t))) 
        
%         if (double(P_sav(i,j,t)) <= double(P_bench(i,j,t)))
        if (P_bench(i,j,t) >= P_sav(i,j,t))
%           i,j,t
%           fprintf("\n");
          feasib_P = (feasib_P)+1;
        end  
      end
    end
  end

end

feasib_P;

fprintf('\n');

if (satis_P==i_up*j_up*t_up) && (feasib_P==i_up*j_up*t_up)
  consec_P = 1+consec_P;
else
  consec_P = 0;
end

% tik = tik-1;

end

fprintf('\n');
fprintf('\n');
fprintf('\n');


%-------case 6.3------

'-------Case 6.3-----'

IniP=100;

%initially price Pijt all the same
P=zeros(20,20,40);%create 20*20*40 matrix
for i=1:i_up 
    for j=1:j_up
        for t=1:t_up
          P(i,j,t)=IniP; %0608-initial P should be infinity, I use 2 for test
        end
    end    
end

% Pijt=2;

Iter_count = 0;

% 0517-M is c16 time window, set at 1 gives it relaxation
% 0424-M is c16 time window, set at 3 gives it relaxation
M=1; 

%0515-change all var to continours instead of MIP
X = sdpvar(i_up,j_up,t_up+M); %0424-Xi,j,t+M, cuz time window constraint c16 Xijt is beyond t
%X = intvar(i_up,j_up,t_up+M); %0424-Xi,j,t+M, cuz time window constraint c16 Xijt is beyond t
X_t0 = sdpvar(i_up,j_up,1); %0424-Xi,j,0 means moves initially for c18: inventory balance eq
%X_t0 = intvar(i_up,j_up,1); %0424-Xi,j,0 means moves initially for c18: inventory balance eq
U = sdpvar(i_up,j_up,t_up);
%U = intvar(i_up,j_up,t_up); 
U_t0 = sdpvar(i_up,j_up,1); %0424-this is for Ui,j,0 
%U_t0 = intvar(i_up,j_up,1); %0424-this is for Ui,j,0 
V = sdpvar(i_up,t_up); 
%V = intvar(i_up,t_up); 
V_t0 = sdpvar(i_up,1); %0424-this is especially for Vi,0
%V_t0 = intvar(i_up,1); %0424-this is especially for Vi,0
D_ij_SAV = sdpvar(i_up,j_up,t_up); %0425-this is to show D_ij_SAV, for futher plot, then find upper bound of Pijt
%D_ij_SAV = intvar(i_up,j_up,t_up); %0425-this is to show D_ij_SAV, for futher plot, then find upper bound of Pijt
% 0604-k ratio
k = sdpvar(i_up,j_up,t_up);

% 070220-set place for Wijt
W=zeros(20,20,40);%create 20*20*40 matrix

%1.1-parameters
% a=b=m=n=1,T=40
a=1;
b=1;
m=1;
n=1;
T=8; %070120-total 40 time periods, 5 time block, so every block is 8
T_abs=8; %070120-T_abs should be same as T
% cijt=40, pijt=10, hit=10, alphaijt=0.1 (all c,p,h,alpha are subjectively to be same)
cijt=20;
pijt=30; % 070120-increase for penalty of unmet demand
% 0425-intentionally make pijt this small to avoid all 0 in Uijt solution
hit=10;
alpha=1;

%%%-----0607----
%%% How to add loop?
%%%----

% 0610-here starts the loop

P_sav = zeros(20,20,40);
U_sav = zeros(20,20,40);
min_k = 0;
P_bench = zeros(20,20,40);

% tik = 2;

satis_P = 0;
feasib_P = 0;
consec_P = 0;

% while (satis_P<i_up*j_up*t_up) || (feasib_P<i_up*j_up*t_up)
% while Iter_count < 2
while (consec_P<3)

Iter_count = Iter_count+1; 


for i=1:i_up 
    for j=1:j_up
        for t=1:t_up
          P_sav(i,j,t)=P(i,j,t); %0608-initial P should be infinity, I use 2 for test
        end
    end    
end  


for i=1:i_up 
    for j=1:j_up
        for t=1:t_up
            U_sav(i,j,t)=U(i,j,t); 
        end
    end    
end


min_k_sav = min_k;


%1.2-add operational level constraints

C = [];

% c15
for i=1:i_up %c15
   for j=1:j_up
     for t=1:t_up
       if (t==1)
         C = [C, U(i,j,1) >= U_t0(i,j,1)+(b*d(i,j,1)+cijt*d(i,j,1)*m-n*T_abs*U(i,j,1))./(a+b+m*(cijt+P(i,j,t)))-X(i,j,1)]; %c15: t=1
       else       
         C = [C, U(i,j,t) >= U(i,j,t-1)+(b*d(i,j,t)+cijt*d(i,j,t)*m-n*T_abs*U(i,j,t))./(a+b+m*(cijt+P(i,j,t)))-X(i,j,t)]; %c15: t from 2
       end
     end
   end    
end

% 0424-c16-I suppose time window is M period for c16 - relax this one, hope get solution
% this also considers Xijt beyond T
for i=1:i_up %c16
    for j=1:j_up
      for t=1:t_up
        if (t==1)
          C = [C, sum(X(i,j,t:(t+M))) >= U_t0(i,j,1)+(b*d(i,j,t)+cijt*d(i,j,t)*m-n*T_abs*U(i,j,t))./(a+b+m*(cijt+P(i,j,t)))]; %c16: t=1
        else
          C = [C, sum(X(i,j,t:(t+M))) >= U(i,j,t-1)+(b*d(i,j,t)+cijt*d(i,j,t)*m-n*T_abs*U(i,j,t))./(a+b+m*(cijt+P(i,j,t)))]; %c16: t from 2
        end
      end  
    end     
end

% 0424 - c17 is in upper level, no need in lower level
% for i=1:i_up %c17
%   for j=1:j_up
%     for t=1:t_up
%       C = [C, ((b-1)*d(i,j,t)+cijt*d(i,j,t)*m)./(n*T_abs) <= U(i,j,t)];
%       C = [C, U(i,j,t) <= (b*d(i,j,t)+cijt*d(i,j,t)*m)./(n*T_abs)];
%     end
%   end
% end

% 0424-c18-here makes <= for Vi,t cuz alpha all set to 1, too relaxed, should be variaty
for i=1:i_up %c18
  for t=1:t_up
    if (t==1)
      C = [C, V(i,t) <= V_t0(i,1) + sum(alpha*X_t0(1:j_up,i,1)) - sum(X(i,1:j_up,t)) ]; %c18: t=1
    else
      C = [C, V(i,t) <= V(i,t-1) + sum(sum(alpha*X(1:j_up,i,1:(t-1)))) - sum(X(i,1:j_up,t)) ]; %c18: t from 2 (remember to use 2 sum for 2nd)
    end
  end
end

% 0508-c19-not added, cuz it looks already relaxed in math function, hidden property of our formula

for i=1:i_up %c20: Xijt
    for j=1:j_up
        for t=1:(t_up+M)
          C = [C, X(i,j,t) >= 0];
        end
    end    
end    

for i=1:i_up %c20: X_t0(i,j,1) is Xijt=0
    for j=1:j_up
        C = [C, X_t0(i,j,1) >= 0];
    end    
end    

for i=1:i_up %c20: Uijt
    for j=1:j_up
        for t=1:t_up
            C = [C, U(i,j,t) >= 0];
        end
    end    
end    

for i=1:i_up %c20: U_t0(i,j,1) is Uijt=0
    for j=1:j_up
        C = [C, U_t0(i,j,1) >= 0];
    end    
end 

for i=1:i_up %c20: Vit
    for t=1:t_up
      C = [C, V(i,t) >= 0];
    end    
end    

for i=1:i_up %c20: V_t0(i,1) is V(i,0)
  C = [C, V_t0(i,1) >= 0];
end

%1.3-solve operational level: settings: use cplex as solver

%0515-following is set optimality for MIP, but kang said make it LNContinuous-https://www.ibm.com/support/knowledgecenter/SSSA5P_12.9.0/ilog.odms.cplex.help/CPLEX/Parameters/topics/OptimalityTarget.html
%0515-'cplex.optimalitytarget'=1: Searches for a globally optimal solution to a convex model.
ops = sdpsettings ('solver','cplex','verbose',2,'cplex.optimalitytarget',1);%0508-verbose 2 shows iterations, 1 means shown normal message, 0 shows silent message
%0515-MIP limit tolerance gap to 1-https://www.ibm.com/support/knowledgecenter/SSSA5P_12.9.0/ilog.odms.cplex.help/CPLEX/Parameters/topics/EpGap.html
%ops = sdpsettings ('solver','cplex','verbose',2,'cplex.mip.tolerances.mipgap',1);%0508-verbose 2 shows iterations, 1 means shown normal message, 0 shows silent message

%0306-limit barrier iteration-https://www.ibm.com/support/knowledgecenter/SSSA5P_12.9.0/ilog.odms.cplex.help/CPLEX/Parameters/topics/BarItLim.html
%Cplex.Param.barrier.limits.iteration = 4;[not work]

%0306-limit network simplex iteration-https://www.ibm.com/support/knowledgecenter/SSSA5P_12.9.0/ilog.odms.cplex.help/CPLEX/Parameters/topics/NetItLim.html
%Cplex.Param.network.iterations = 2;[not work]

%ops.mip.limits.nodes = 4;[not work]

%0304-use the online method for out-of-memory: https://www.ilovematlab.cn/thread-266244-1-1.html
%ops.mip.strategy.file = 3; %[not work: should put in sdpsetting] this has same results 1159110 node, but reduces time for 3 min.

%0305-use online method for out-of-memory: https://www.ibm.com/developerworks/community/forums/html/topic?id=e3ea344c-7249-4edd-a47d-29e1f80fa480
%Cplex.Param.mip.limits.auxrootthreads = 2; %[not work: should put in sdpsetting] -1 or 1

%obj-operational obj z_op
z_op = sum(sum(sum(cijt.*X(1:i_up,1:j_up,1:t_up)))) + sum(sum(cijt.*X_t0(1:i_up,1:j_up,1)))...
  + sum(sum(sum(pijt.*U(1:i_up,1:j_up,1:t_up)))) + sum(sum(pijt.*U_t0(1:i_up,1:j_up,1)))...
  + sum(sum(hit.*V(1:i_up,1:t_up))) + sum(hit.*V_t0(1:i_up,1));

%1.35-express D_ij_SAV
for i=1:i_up
  for j=1:j_up
    for t=1:t_up
      D_ij_SAV(i,j,t) = (b*d(i,j,t)+cijt*d(i,j,t)*m-n*T_abs*U(i,j,t))./(a+b+m*(cijt+P(i,j,t)));
    end
  end
end

% 0604-calculate minimum k ratio covered by SAVs
for i=1:i_up
  for j=1:j_up
    for t=1:t_up
      k(i,j,t) = D_ij_SAV(i,j,t)./d(i,j,t);
    end
  end
end

% k_min = min(k(1,1,1),k(1,2,1),k(1,3,1))

% k_min = 0.99;
% % 
% if ((k_min) >= (k(1,1,1)))
%   k_min = k(1,1,1)                
% end  

% 0604-find minimum covered ratio k
% for i=1:i_up
%   for j=1:j_up
%     for t=1:t_up
%       if (k_min >= k(i,j,t))
%         k_min = k(i,j,t)                
%       end  
%     end
%   end
% end

%1.4-check if solved and output
result = optimize(C,z_op,ops); % 0508-origin
% result = solvesdp(C,z_op,ops); % 0508-this works

%result = optimize(C,z);
% if result.problem == 0 % problem =0 means solve successfully,only print Uijt here

%   fprintf('\nk value: ') %ratio k covered by SAVs
%   value(k)

  
%   fprintf('\nP_Ini_(1,1,1) = ')
%   value(P(1,1,1))

%   fprintf('\nMin cost by Lower = ')
%   value(z_op)

  all_k = [];

  for i=1:i_up
    for j=1:j_up
      for t=1:t_up
        all_k(end+1) = k(i,j,t);
      end
    end
  end
  
  min_k = min(all_k);
  
  k_avg = mean(all_k);
  
%   fprintf('\n min k: ')
%   value(min_k)
  
  max_k = max(all_k);
  
%   fprintf('\n max k: ')
%   value(max_k)
  
  % 0615-save max_k in txt

%   fid = fopen(['IniP=',num2str(2),'max_k.txt'],'a'); %save to last 'a'
%   fprintf(fid,'%.4f\t',double(max_k)); %0516-keep 4 decimals
%   fclose(fid);
  
%   fprintf('\nU(1,1,1) value: ')
%   value(U(1,1,1))
        
%   fprintf('\nX(1,1,1) value: ') 
%   value(X(1,1,1))

fprintf('\n Iteration No. = ')
value(Iter_count)



% fprintf('\n D_ij_SAV(1,1,1) value: ') %0515-as for now hide D_ij_SAV
% value(D_ij_SAV(1,1,1))
  
%   fprintf('\n k(1,1,1) ')
%   value(k(1,1,1))

%   fprintf('\n')
%     fprintf('\nV value: ')
%     value(V)
%     
%     fprintf('\nU_t0 value: ')
%     value(U_t0)
%     
%     fprintf('\nX_t0 value: ') 
%     value(X_t0)
%     
%     fprintf('\nV_t0 value: ')
%     value(V_t0)

% else
%   disp('Not solved!')
% 
% end

% 0519-save total cost z_op in txt
% fid = fopen(['Iter=',num2str(Iter_count),'_OperCost.txt'],'wt'); %creat & write z_op to txt
% fprintf(fid,'%.4f\n',double(z_op)); %0516-keep 4 decimals
% fclose(fid);

% 0605-save min ratio k covered by SAVs in txt
% fid_ratiok = fopen(['Iter=',num2str(Iter_count),'_MinRatiok.txt'],'wt'); %creat & write z_op to txt
% fprintf(fid_ratiok,'%.4f\n',double(min_k)); %0604-keep 4 decimals
% fclose(fid_ratiok);

% 0517-following is save Uijt in txt
%1.6-save Uijt as txt for upper use
% for U_t=1:t_up 
%      fid = fopen(['IniP=',num2str(IniP),'Iter=',num2str(Iter_count),'_Uijt=',num2str(U_t),'.txt'],'wt'); %creat & write Uijt to txt
%      for i=1:i_up
%          for j=1:j_up
%              if j == j_up
%                  fprintf(fid,'%.4f\n',double(U(i,j,U_t))); %0516-keep 4 decimals
%              else
%                  fprintf(fid,'%.4f\t',double(U(i,j,U_t))); %0516-keep 4 decimals
%              end     
%          end    
%      end    
%      fclose(fid);
% end

% 0607-calculate Upper Pijt upper bound
for i=1:i_up
  for j=1:j_up
    for t=1:t_up
      P(i,j,t)=(b*d(i,j,t)+cijt*d(i,j,t)*m-n*T_abs*U(i,j,t)-(a+b)*min_k*d(i,j,t))./(m*min_k*d(i,j,t))-cijt;
    end
  end
end

% fprintf('\nUpper bound P(1,1,1) = %.4f\n',double(P(1,1,1)));

%0608-----
%calculate Upper revenue
%0608-----

revenue=0;

for i=1:i_up
  for j=1:j_up
    for t=1:t_up
      if (t==1)||(t==2)||(t==9)||(t==10)||(t==17)||(t==18)||(t==31)||(t==32)||(t==39)||(t==40)
%         T_abs=2;
        revenue = revenue-(P(i,j,t).*(b.*d(i,j,t)+cijt.*d(i,j,t)*m-n*T_abs*U(i,j,t)))./(a+b+m.*(cijt+P(i,j,t))); 
      else
%         T_abs=2; %here should be 3, but temporarily set to 2, aligh with Oper level
        revenue = revenue-(P(i,j,t).*(b.*d(i,j,t)+cijt.*d(i,j,t)*m-n*T_abs*U(i,j,t)))./(a+b+m.*(cijt+P(i,j,t)));
      end
    end
  end
end

% 070320-avg. all P

 all_P = [];

  for i=1:i_up
    for j=1:j_up
      for t=1:t_up
        all_P = P(i,j,t);
      end
    end
  end

avg_P = all_P./(i_up*j_up*t_up);



% 070320-avg. all D

 all_Dijt = [];

  for i=1:i_up
    for j=1:j_up
      for t=1:t_up
        all_Dijt = D_ij_SAV(i,j,t);
      end
    end
  end

avg_Dijt = all_Dijt./(i_up*j_up*t_up);

% 070320-print results

fprintf('\n Upper bound P(1,1,1) = %.4f\n',double(P(1,1,1))); % particular P each iteration

fprintf('\n');

fprintf('\n D(1,1,1) = %.4f\n',double(D_ij_SAV(1,1,1))); % particular Dijt each iteration

fprintf('\n');

fprintf('\n Avg_P = %.4f\n',double(avg_P)); % particular avg_P each iteration

fprintf('\n');

fprintf('\n Avg_Dijt = %.4f\n',double(avg_Dijt)); % particular avg_Dijt each iteration

fprintf('\n');

fprintf('\n Max Revenue by Upper= %.4f\n',double(-revenue));

fprintf('\n');

fprintf('\n Min Cost by Lower= %.4f\n',double(z_op));

fprintf('\n');

fprintf('\n Profit = %.4f\n',double(-revenue-z_op));

fprintf('\n');

% test how many new P > old P

satis_P = 0;

for i=1:i_up
  for j=1:j_up
    for t=1:t_up
%       if (double(P_sav(i,j,t))<=double(P(i,j,t))) % 0610-check new P greater than old ones
        %0610-if new P no more than 1 cent than old, than add to satis_P terminate condition
      %if ((double(P_sav(i,j,t)))*1.0001 >= double(P(i,j,t))) 
%       if ((double(P_sav(i,j,t)))+0.01 >= double(P(i,j,t))) 
      if ((P_sav(i,j,t))+0.01 > P(i,j,t)) 
        satis_P = (satis_P)+1;
      end  
    end
  end
end

satis_P;

fprintf('\n');

% test how many new P all feasible that new_P <= new P upper bound condition


if (Iter_count==1)
  feasib_P = 0;
else
  
  feasib_P = 0;
  
  for i=1:i_up
    for j=1:j_up
      for t=1:t_up
        P_bench(i,j,t) = (b*d(i,j,t)+cijt*d(i,j,t)*m-n*T_abs*U(i,j,t)-(a+b)*min_k*d(i,j,t))./(m*min_k*d(i,j,t))-cijt;        
      end
    end
  end
  
    
  for i=1:i_up
    for j=1:j_up
      for t=1:t_up
  %       if (double(P_sav(i,j,t))<=double(P(i,j,t))) % 0610-check new P greater than old ones
          %0610-if new P no more than 1 cent than old, than add to satis_P terminate condition
        %if ((double(P_sav(i,j,t)))*1.0001 >= double(P(i,j,t))) 
        
%         if (double(P_sav(i,j,t)) <= double(P_bench(i,j,t)))
        if (P_bench(i,j,t) >= P_sav(i,j,t))
%           i,j,t
%           fprintf("\n");
          feasib_P = (feasib_P)+1;
        end  
      end
    end
  end

end

feasib_P;

fprintf('\n');

if (satis_P==i_up*j_up*t_up) && (feasib_P==i_up*j_up*t_up)
  consec_P = 1+consec_P;
else
  consec_P = 0;
end

% tik = tik-1;

end

fprintf('\n');
fprintf('\n');
fprintf('\n');




%-------case 7.1------

'------Case 7.1------'

IniP=1;

%initially price Pijt all the same
P=zeros(20,20,40);%create 20*20*40 matrix
for i=1:i_up 
    for j=1:j_up
        for t=1:t_up
          P(i,j,t)=IniP; %0608-initial P should be infinity, I use 2 for test
        end
    end    
end

% Pijt=2;

Iter_count = 0;

% 0517-M is c16 time window, set at 1 gives it relaxation
% 0424-M is c16 time window, set at 3 gives it relaxation
M=1; 

%0515-change all var to continours instead of MIP
X = sdpvar(i_up,j_up,t_up+M); %0424-Xi,j,t+M, cuz time window constraint c16 Xijt is beyond t
%X = intvar(i_up,j_up,t_up+M); %0424-Xi,j,t+M, cuz time window constraint c16 Xijt is beyond t
X_t0 = sdpvar(i_up,j_up,1); %0424-Xi,j,0 means moves initially for c18: inventory balance eq
%X_t0 = intvar(i_up,j_up,1); %0424-Xi,j,0 means moves initially for c18: inventory balance eq
U = sdpvar(i_up,j_up,t_up);
%U = intvar(i_up,j_up,t_up); 
U_t0 = sdpvar(i_up,j_up,1); %0424-this is for Ui,j,0 
%U_t0 = intvar(i_up,j_up,1); %0424-this is for Ui,j,0 
V = sdpvar(i_up,t_up); 
%V = intvar(i_up,t_up); 
V_t0 = sdpvar(i_up,1); %0424-this is especially for Vi,0
%V_t0 = intvar(i_up,1); %0424-this is especially for Vi,0
D_ij_SAV = sdpvar(i_up,j_up,t_up); %0425-this is to show D_ij_SAV, for futher plot, then find upper bound of Pijt
%D_ij_SAV = intvar(i_up,j_up,t_up); %0425-this is to show D_ij_SAV, for futher plot, then find upper bound of Pijt
% 0604-k ratio
k = sdpvar(i_up,j_up,t_up);

% 070220-set place for Wijt
W=zeros(20,20,40);%create 20*20*40 matrix

%1.1-parameters
% a=b=m=n=1,T=40
a=1;
b=1;
m=1;
n=1;
T=8; %070120-total 40 time periods, 5 time block, so every block is 8
T_abs=8; %070120-T_abs should be same as T
% cijt=40, pijt=10, hit=10, alphaijt=0.1 (all c,p,h,alpha are subjectively to be same)
cijt=20;
pijt=10; % 070120-increase for penalty of unmet demand
% 0425-intentionally make pijt this small to avoid all 0 in Uijt solution
hit=30;
alpha=1;

%%%-----0607----
%%% How to add loop?
%%%----

% 0610-here starts the loop

P_sav = zeros(20,20,40);
U_sav = zeros(20,20,40);
min_k = 0;
P_bench = zeros(20,20,40);

% tik = 2;

satis_P = 0;
feasib_P = 0;
consec_P = 0;

% while (satis_P<i_up*j_up*t_up) || (feasib_P<i_up*j_up*t_up)
% while Iter_count < 2
while (consec_P<3)

Iter_count = Iter_count+1; 


for i=1:i_up 
    for j=1:j_up
        for t=1:t_up
          P_sav(i,j,t)=P(i,j,t); %0608-initial P should be infinity, I use 2 for test
        end
    end    
end  


for i=1:i_up 
    for j=1:j_up
        for t=1:t_up
            U_sav(i,j,t)=U(i,j,t); 
        end
    end    
end


min_k_sav = min_k;


%1.2-add operational level constraints

C = [];

% c15
for i=1:i_up %c15
   for j=1:j_up
     for t=1:t_up
       if (t==1)
         C = [C, U(i,j,1) >= U_t0(i,j,1)+(b*d(i,j,1)+cijt*d(i,j,1)*m-n*T_abs*U(i,j,1))./(a+b+m*(cijt+P(i,j,t)))-X(i,j,1)]; %c15: t=1
       else       
         C = [C, U(i,j,t) >= U(i,j,t-1)+(b*d(i,j,t)+cijt*d(i,j,t)*m-n*T_abs*U(i,j,t))./(a+b+m*(cijt+P(i,j,t)))-X(i,j,t)]; %c15: t from 2
       end
     end
   end    
end

% 0424-c16-I suppose time window is M period for c16 - relax this one, hope get solution
% this also considers Xijt beyond T
for i=1:i_up %c16
    for j=1:j_up
      for t=1:t_up
        if (t==1)
          C = [C, sum(X(i,j,t:(t+M))) >= U_t0(i,j,1)+(b*d(i,j,t)+cijt*d(i,j,t)*m-n*T_abs*U(i,j,t))./(a+b+m*(cijt+P(i,j,t)))]; %c16: t=1
        else
          C = [C, sum(X(i,j,t:(t+M))) >= U(i,j,t-1)+(b*d(i,j,t)+cijt*d(i,j,t)*m-n*T_abs*U(i,j,t))./(a+b+m*(cijt+P(i,j,t)))]; %c16: t from 2
        end
      end  
    end     
end

% 0424 - c17 is in upper level, no need in lower level
% for i=1:i_up %c17
%   for j=1:j_up
%     for t=1:t_up
%       C = [C, ((b-1)*d(i,j,t)+cijt*d(i,j,t)*m)./(n*T_abs) <= U(i,j,t)];
%       C = [C, U(i,j,t) <= (b*d(i,j,t)+cijt*d(i,j,t)*m)./(n*T_abs)];
%     end
%   end
% end

% 0424-c18-here makes <= for Vi,t cuz alpha all set to 1, too relaxed, should be variaty
for i=1:i_up %c18
  for t=1:t_up
    if (t==1)
      C = [C, V(i,t) <= V_t0(i,1) + sum(alpha*X_t0(1:j_up,i,1)) - sum(X(i,1:j_up,t)) ]; %c18: t=1
    else
      C = [C, V(i,t) <= V(i,t-1) + sum(sum(alpha*X(1:j_up,i,1:(t-1)))) - sum(X(i,1:j_up,t)) ]; %c18: t from 2 (remember to use 2 sum for 2nd)
    end
  end
end

% 0508-c19-not added, cuz it looks already relaxed in math function, hidden property of our formula

for i=1:i_up %c20: Xijt
    for j=1:j_up
        for t=1:(t_up+M)
          C = [C, X(i,j,t) >= 0];
        end
    end    
end    

for i=1:i_up %c20: X_t0(i,j,1) is Xijt=0
    for j=1:j_up
        C = [C, X_t0(i,j,1) >= 0];
    end    
end    

for i=1:i_up %c20: Uijt
    for j=1:j_up
        for t=1:t_up
            C = [C, U(i,j,t) >= 0];
        end
    end    
end    

for i=1:i_up %c20: U_t0(i,j,1) is Uijt=0
    for j=1:j_up
        C = [C, U_t0(i,j,1) >= 0];
    end    
end 

for i=1:i_up %c20: Vit
    for t=1:t_up
      C = [C, V(i,t) >= 0];
    end    
end    

for i=1:i_up %c20: V_t0(i,1) is V(i,0)
  C = [C, V_t0(i,1) >= 0];
end

%1.3-solve operational level: settings: use cplex as solver

%0515-following is set optimality for MIP, but kang said make it LNContinuous-https://www.ibm.com/support/knowledgecenter/SSSA5P_12.9.0/ilog.odms.cplex.help/CPLEX/Parameters/topics/OptimalityTarget.html
%0515-'cplex.optimalitytarget'=1: Searches for a globally optimal solution to a convex model.
ops = sdpsettings ('solver','cplex','verbose',2,'cplex.optimalitytarget',1);%0508-verbose 2 shows iterations, 1 means shown normal message, 0 shows silent message
%0515-MIP limit tolerance gap to 1-https://www.ibm.com/support/knowledgecenter/SSSA5P_12.9.0/ilog.odms.cplex.help/CPLEX/Parameters/topics/EpGap.html
%ops = sdpsettings ('solver','cplex','verbose',2,'cplex.mip.tolerances.mipgap',1);%0508-verbose 2 shows iterations, 1 means shown normal message, 0 shows silent message

%0306-limit barrier iteration-https://www.ibm.com/support/knowledgecenter/SSSA5P_12.9.0/ilog.odms.cplex.help/CPLEX/Parameters/topics/BarItLim.html
%Cplex.Param.barrier.limits.iteration = 4;[not work]

%0306-limit network simplex iteration-https://www.ibm.com/support/knowledgecenter/SSSA5P_12.9.0/ilog.odms.cplex.help/CPLEX/Parameters/topics/NetItLim.html
%Cplex.Param.network.iterations = 2;[not work]

%ops.mip.limits.nodes = 4;[not work]

%0304-use the online method for out-of-memory: https://www.ilovematlab.cn/thread-266244-1-1.html
%ops.mip.strategy.file = 3; %[not work: should put in sdpsetting] this has same results 1159110 node, but reduces time for 3 min.

%0305-use online method for out-of-memory: https://www.ibm.com/developerworks/community/forums/html/topic?id=e3ea344c-7249-4edd-a47d-29e1f80fa480
%Cplex.Param.mip.limits.auxrootthreads = 2; %[not work: should put in sdpsetting] -1 or 1

%obj-operational obj z_op
z_op = sum(sum(sum(cijt.*X(1:i_up,1:j_up,1:t_up)))) + sum(sum(cijt.*X_t0(1:i_up,1:j_up,1)))...
  + sum(sum(sum(pijt.*U(1:i_up,1:j_up,1:t_up)))) + sum(sum(pijt.*U_t0(1:i_up,1:j_up,1)))...
  + sum(sum(hit.*V(1:i_up,1:t_up))) + sum(hit.*V_t0(1:i_up,1));

%1.35-express D_ij_SAV
for i=1:i_up
  for j=1:j_up
    for t=1:t_up
      D_ij_SAV(i,j,t) = (b*d(i,j,t)+cijt*d(i,j,t)*m-n*T_abs*U(i,j,t))./(a+b+m*(cijt+P(i,j,t)));
    end
  end
end

% 0604-calculate minimum k ratio covered by SAVs
for i=1:i_up
  for j=1:j_up
    for t=1:t_up
      k(i,j,t) = D_ij_SAV(i,j,t)./d(i,j,t);
    end
  end
end

% k_min = min(k(1,1,1),k(1,2,1),k(1,3,1))

% k_min = 0.99;
% % 
% if ((k_min) >= (k(1,1,1)))
%   k_min = k(1,1,1)                
% end  

% 0604-find minimum covered ratio k
% for i=1:i_up
%   for j=1:j_up
%     for t=1:t_up
%       if (k_min >= k(i,j,t))
%         k_min = k(i,j,t)                
%       end  
%     end
%   end
% end

%1.4-check if solved and output
result = optimize(C,z_op,ops); % 0508-origin
% result = solvesdp(C,z_op,ops); % 0508-this works

%result = optimize(C,z);
% if result.problem == 0 % problem =0 means solve successfully,only print Uijt here

%   fprintf('\nk value: ') %ratio k covered by SAVs
%   value(k)

  
%   fprintf('\nP_Ini_(1,1,1) = ')
%   value(P(1,1,1))

%   fprintf('\nMin cost by Lower = ')
%   value(z_op)

  all_k = [];

  for i=1:i_up
    for j=1:j_up
      for t=1:t_up
        all_k(end+1) = k(i,j,t);
      end
    end
  end
  
  min_k = min(all_k);
  
  k_avg = mean(all_k);
  
%   fprintf('\n min k: ')
%   value(min_k)
  
  max_k = max(all_k);
  
%   fprintf('\n max k: ')
%   value(max_k)
  
  % 0615-save max_k in txt

%   fid = fopen(['IniP=',num2str(2),'max_k.txt'],'a'); %save to last 'a'
%   fprintf(fid,'%.4f\t',double(max_k)); %0516-keep 4 decimals
%   fclose(fid);
  
%   fprintf('\nU(1,1,1) value: ')
%   value(U(1,1,1))
        
%   fprintf('\nX(1,1,1) value: ') 
%   value(X(1,1,1))

fprintf('\n Iteration No. = ')
value(Iter_count)



% fprintf('\n D_ij_SAV(1,1,1) value: ') %0515-as for now hide D_ij_SAV
% value(D_ij_SAV(1,1,1))
  
%   fprintf('\n k(1,1,1) ')
%   value(k(1,1,1))

%   fprintf('\n')
%     fprintf('\nV value: ')
%     value(V)
%     
%     fprintf('\nU_t0 value: ')
%     value(U_t0)
%     
%     fprintf('\nX_t0 value: ') 
%     value(X_t0)
%     
%     fprintf('\nV_t0 value: ')
%     value(V_t0)

% else
%   disp('Not solved!')
% 
% end

% 0519-save total cost z_op in txt
% fid = fopen(['Iter=',num2str(Iter_count),'_OperCost.txt'],'wt'); %creat & write z_op to txt
% fprintf(fid,'%.4f\n',double(z_op)); %0516-keep 4 decimals
% fclose(fid);

% 0605-save min ratio k covered by SAVs in txt
% fid_ratiok = fopen(['Iter=',num2str(Iter_count),'_MinRatiok.txt'],'wt'); %creat & write z_op to txt
% fprintf(fid_ratiok,'%.4f\n',double(min_k)); %0604-keep 4 decimals
% fclose(fid_ratiok);

% 0517-following is save Uijt in txt
%1.6-save Uijt as txt for upper use
% for U_t=1:t_up 
%      fid = fopen(['IniP=',num2str(IniP),'Iter=',num2str(Iter_count),'_Uijt=',num2str(U_t),'.txt'],'wt'); %creat & write Uijt to txt
%      for i=1:i_up
%          for j=1:j_up
%              if j == j_up
%                  fprintf(fid,'%.4f\n',double(U(i,j,U_t))); %0516-keep 4 decimals
%              else
%                  fprintf(fid,'%.4f\t',double(U(i,j,U_t))); %0516-keep 4 decimals
%              end     
%          end    
%      end    
%      fclose(fid);
% end

% 0607-calculate Upper Pijt upper bound
for i=1:i_up
  for j=1:j_up
    for t=1:t_up
      P(i,j,t)=(b*d(i,j,t)+cijt*d(i,j,t)*m-n*T_abs*U(i,j,t)-(a+b)*min_k*d(i,j,t))./(m*min_k*d(i,j,t))-cijt;
    end
  end
end

% fprintf('\nUpper bound P(1,1,1) = %.4f\n',double(P(1,1,1)));

%0608-----
%calculate Upper revenue
%0608-----

revenue=0;

for i=1:i_up
  for j=1:j_up
    for t=1:t_up
      if (t==1)||(t==2)||(t==9)||(t==10)||(t==17)||(t==18)||(t==31)||(t==32)||(t==39)||(t==40)
%         T_abs=2;
        revenue = revenue-(P(i,j,t).*(b.*d(i,j,t)+cijt.*d(i,j,t)*m-n*T_abs*U(i,j,t)))./(a+b+m.*(cijt+P(i,j,t))); 
      else
%         T_abs=2; %here should be 3, but temporarily set to 2, aligh with Oper level
        revenue = revenue-(P(i,j,t).*(b.*d(i,j,t)+cijt.*d(i,j,t)*m-n*T_abs*U(i,j,t)))./(a+b+m.*(cijt+P(i,j,t)));
      end
    end
  end
end

% 070320-avg. all P

 all_P = [];

  for i=1:i_up
    for j=1:j_up
      for t=1:t_up
        all_P = P(i,j,t);
      end
    end
  end

avg_P = all_P./(i_up*j_up*t_up);



% 070320-avg. all D

 all_Dijt = [];

  for i=1:i_up
    for j=1:j_up
      for t=1:t_up
        all_Dijt = D_ij_SAV(i,j,t);
      end
    end
  end

avg_Dijt = all_Dijt./(i_up*j_up*t_up);

% 070320-print results

fprintf('\n Upper bound P(1,1,1) = %.4f\n',double(P(1,1,1))); % particular P each iteration

fprintf('\n');

fprintf('\n D(1,1,1) = %.4f\n',double(D_ij_SAV(1,1,1))); % particular Dijt each iteration

fprintf('\n');

fprintf('\n Avg_P = %.4f\n',double(avg_P)); % particular avg_P each iteration

fprintf('\n');

fprintf('\n Avg_Dijt = %.4f\n',double(avg_Dijt)); % particular avg_Dijt each iteration

fprintf('\n');

fprintf('\n Max Revenue by Upper= %.4f\n',double(-revenue));

fprintf('\n');

fprintf('\n Min Cost by Lower= %.4f\n',double(z_op));

fprintf('\n');

fprintf('\n Profit = %.4f\n',double(-revenue-z_op));

fprintf('\n');

% test how many new P > old P

satis_P = 0;

for i=1:i_up
  for j=1:j_up
    for t=1:t_up
%       if (double(P_sav(i,j,t))<=double(P(i,j,t))) % 0610-check new P greater than old ones
        %0610-if new P no more than 1 cent than old, than add to satis_P terminate condition
      %if ((double(P_sav(i,j,t)))*1.0001 >= double(P(i,j,t))) 
%       if ((double(P_sav(i,j,t)))+0.01 >= double(P(i,j,t))) 
      if ((P_sav(i,j,t))+0.01 > P(i,j,t)) 
        satis_P = (satis_P)+1;
      end  
    end
  end
end

satis_P;

fprintf('\n');

% test how many new P all feasible that new_P <= new P upper bound condition


if (Iter_count==1)
  feasib_P = 0;
else
  
  feasib_P = 0;
  
  for i=1:i_up
    for j=1:j_up
      for t=1:t_up
        P_bench(i,j,t) = (b*d(i,j,t)+cijt*d(i,j,t)*m-n*T_abs*U(i,j,t)-(a+b)*min_k*d(i,j,t))./(m*min_k*d(i,j,t))-cijt;        
      end
    end
  end
  
    
  for i=1:i_up
    for j=1:j_up
      for t=1:t_up
  %       if (double(P_sav(i,j,t))<=double(P(i,j,t))) % 0610-check new P greater than old ones
          %0610-if new P no more than 1 cent than old, than add to satis_P terminate condition
        %if ((double(P_sav(i,j,t)))*1.0001 >= double(P(i,j,t))) 
        
%         if (double(P_sav(i,j,t)) <= double(P_bench(i,j,t)))
        if (P_bench(i,j,t) >= P_sav(i,j,t))
%           i,j,t
%           fprintf("\n");
          feasib_P = (feasib_P)+1;
        end  
      end
    end
  end

end

feasib_P;

fprintf('\n');

if (satis_P==i_up*j_up*t_up) && (feasib_P==i_up*j_up*t_up)
  consec_P = 1+consec_P;
else
  consec_P = 0;
end

% tik = tik-1;

end

fprintf('\n');
fprintf('\n');
fprintf('\n');


%-------case 7.2------

'-------Case 7.2-------'

IniP=10;

%initially price Pijt all the same
P=zeros(20,20,40);%create 20*20*40 matrix
for i=1:i_up 
    for j=1:j_up
        for t=1:t_up
          P(i,j,t)=IniP; %0608-initial P should be infinity, I use 2 for test
        end
    end    
end

% Pijt=2;

Iter_count = 0;

% 0517-M is c16 time window, set at 1 gives it relaxation
% 0424-M is c16 time window, set at 3 gives it relaxation
M=1; 

%0515-change all var to continours instead of MIP
X = sdpvar(i_up,j_up,t_up+M); %0424-Xi,j,t+M, cuz time window constraint c16 Xijt is beyond t
%X = intvar(i_up,j_up,t_up+M); %0424-Xi,j,t+M, cuz time window constraint c16 Xijt is beyond t
X_t0 = sdpvar(i_up,j_up,1); %0424-Xi,j,0 means moves initially for c18: inventory balance eq
%X_t0 = intvar(i_up,j_up,1); %0424-Xi,j,0 means moves initially for c18: inventory balance eq
U = sdpvar(i_up,j_up,t_up);
%U = intvar(i_up,j_up,t_up); 
U_t0 = sdpvar(i_up,j_up,1); %0424-this is for Ui,j,0 
%U_t0 = intvar(i_up,j_up,1); %0424-this is for Ui,j,0 
V = sdpvar(i_up,t_up); 
%V = intvar(i_up,t_up); 
V_t0 = sdpvar(i_up,1); %0424-this is especially for Vi,0
%V_t0 = intvar(i_up,1); %0424-this is especially for Vi,0
D_ij_SAV = sdpvar(i_up,j_up,t_up); %0425-this is to show D_ij_SAV, for futher plot, then find upper bound of Pijt
%D_ij_SAV = intvar(i_up,j_up,t_up); %0425-this is to show D_ij_SAV, for futher plot, then find upper bound of Pijt
% 0604-k ratio
k = sdpvar(i_up,j_up,t_up);

% 070220-set place for Wijt
W=zeros(20,20,40);%create 20*20*40 matrix

%1.1-parameters
% a=b=m=n=1,T=40
a=1;
b=1;
m=1;
n=1;
T=8; %070120-total 40 time periods, 5 time block, so every block is 8
T_abs=8; %070120-T_abs should be same as T
% cijt=40, pijt=10, hit=10, alphaijt=0.1 (all c,p,h,alpha are subjectively to be same)
cijt=20;
pijt=10; % 070120-increase for penalty of unmet demand
% 0425-intentionally make pijt this small to avoid all 0 in Uijt solution
hit=30;
alpha=1;

%%%-----0607----
%%% How to add loop?
%%%----

% 0610-here starts the loop

P_sav = zeros(20,20,40);
U_sav = zeros(20,20,40);
min_k = 0;
P_bench = zeros(20,20,40);

% tik = 2;

satis_P = 0;
feasib_P = 0;
consec_P = 0;

% while (satis_P<i_up*j_up*t_up) || (feasib_P<i_up*j_up*t_up)
% while Iter_count < 2
while (consec_P<3)

Iter_count = Iter_count+1; 


for i=1:i_up 
    for j=1:j_up
        for t=1:t_up
          P_sav(i,j,t)=P(i,j,t); %0608-initial P should be infinity, I use 2 for test
        end
    end    
end  


for i=1:i_up 
    for j=1:j_up
        for t=1:t_up
            U_sav(i,j,t)=U(i,j,t); 
        end
    end    
end


min_k_sav = min_k;


%1.2-add operational level constraints

C = [];

% c15
for i=1:i_up %c15
   for j=1:j_up
     for t=1:t_up
       if (t==1)
         C = [C, U(i,j,1) >= U_t0(i,j,1)+(b*d(i,j,1)+cijt*d(i,j,1)*m-n*T_abs*U(i,j,1))./(a+b+m*(cijt+P(i,j,t)))-X(i,j,1)]; %c15: t=1
       else       
         C = [C, U(i,j,t) >= U(i,j,t-1)+(b*d(i,j,t)+cijt*d(i,j,t)*m-n*T_abs*U(i,j,t))./(a+b+m*(cijt+P(i,j,t)))-X(i,j,t)]; %c15: t from 2
       end
     end
   end    
end

% 0424-c16-I suppose time window is M period for c16 - relax this one, hope get solution
% this also considers Xijt beyond T
for i=1:i_up %c16
    for j=1:j_up
      for t=1:t_up
        if (t==1)
          C = [C, sum(X(i,j,t:(t+M))) >= U_t0(i,j,1)+(b*d(i,j,t)+cijt*d(i,j,t)*m-n*T_abs*U(i,j,t))./(a+b+m*(cijt+P(i,j,t)))]; %c16: t=1
        else
          C = [C, sum(X(i,j,t:(t+M))) >= U(i,j,t-1)+(b*d(i,j,t)+cijt*d(i,j,t)*m-n*T_abs*U(i,j,t))./(a+b+m*(cijt+P(i,j,t)))]; %c16: t from 2
        end
      end  
    end     
end

% 0424 - c17 is in upper level, no need in lower level
% for i=1:i_up %c17
%   for j=1:j_up
%     for t=1:t_up
%       C = [C, ((b-1)*d(i,j,t)+cijt*d(i,j,t)*m)./(n*T_abs) <= U(i,j,t)];
%       C = [C, U(i,j,t) <= (b*d(i,j,t)+cijt*d(i,j,t)*m)./(n*T_abs)];
%     end
%   end
% end

% 0424-c18-here makes <= for Vi,t cuz alpha all set to 1, too relaxed, should be variaty
for i=1:i_up %c18
  for t=1:t_up
    if (t==1)
      C = [C, V(i,t) <= V_t0(i,1) + sum(alpha*X_t0(1:j_up,i,1)) - sum(X(i,1:j_up,t)) ]; %c18: t=1
    else
      C = [C, V(i,t) <= V(i,t-1) + sum(sum(alpha*X(1:j_up,i,1:(t-1)))) - sum(X(i,1:j_up,t)) ]; %c18: t from 2 (remember to use 2 sum for 2nd)
    end
  end
end

% 0508-c19-not added, cuz it looks already relaxed in math function, hidden property of our formula

for i=1:i_up %c20: Xijt
    for j=1:j_up
        for t=1:(t_up+M)
          C = [C, X(i,j,t) >= 0];
        end
    end    
end    

for i=1:i_up %c20: X_t0(i,j,1) is Xijt=0
    for j=1:j_up
        C = [C, X_t0(i,j,1) >= 0];
    end    
end    

for i=1:i_up %c20: Uijt
    for j=1:j_up
        for t=1:t_up
            C = [C, U(i,j,t) >= 0];
        end
    end    
end    

for i=1:i_up %c20: U_t0(i,j,1) is Uijt=0
    for j=1:j_up
        C = [C, U_t0(i,j,1) >= 0];
    end    
end 

for i=1:i_up %c20: Vit
    for t=1:t_up
      C = [C, V(i,t) >= 0];
    end    
end    

for i=1:i_up %c20: V_t0(i,1) is V(i,0)
  C = [C, V_t0(i,1) >= 0];
end

%1.3-solve operational level: settings: use cplex as solver

%0515-following is set optimality for MIP, but kang said make it LNContinuous-https://www.ibm.com/support/knowledgecenter/SSSA5P_12.9.0/ilog.odms.cplex.help/CPLEX/Parameters/topics/OptimalityTarget.html
%0515-'cplex.optimalitytarget'=1: Searches for a globally optimal solution to a convex model.
ops = sdpsettings ('solver','cplex','verbose',2,'cplex.optimalitytarget',1);%0508-verbose 2 shows iterations, 1 means shown normal message, 0 shows silent message
%0515-MIP limit tolerance gap to 1-https://www.ibm.com/support/knowledgecenter/SSSA5P_12.9.0/ilog.odms.cplex.help/CPLEX/Parameters/topics/EpGap.html
%ops = sdpsettings ('solver','cplex','verbose',2,'cplex.mip.tolerances.mipgap',1);%0508-verbose 2 shows iterations, 1 means shown normal message, 0 shows silent message

%0306-limit barrier iteration-https://www.ibm.com/support/knowledgecenter/SSSA5P_12.9.0/ilog.odms.cplex.help/CPLEX/Parameters/topics/BarItLim.html
%Cplex.Param.barrier.limits.iteration = 4;[not work]

%0306-limit network simplex iteration-https://www.ibm.com/support/knowledgecenter/SSSA5P_12.9.0/ilog.odms.cplex.help/CPLEX/Parameters/topics/NetItLim.html
%Cplex.Param.network.iterations = 2;[not work]

%ops.mip.limits.nodes = 4;[not work]

%0304-use the online method for out-of-memory: https://www.ilovematlab.cn/thread-266244-1-1.html
%ops.mip.strategy.file = 3; %[not work: should put in sdpsetting] this has same results 1159110 node, but reduces time for 3 min.

%0305-use online method for out-of-memory: https://www.ibm.com/developerworks/community/forums/html/topic?id=e3ea344c-7249-4edd-a47d-29e1f80fa480
%Cplex.Param.mip.limits.auxrootthreads = 2; %[not work: should put in sdpsetting] -1 or 1

%obj-operational obj z_op
z_op = sum(sum(sum(cijt.*X(1:i_up,1:j_up,1:t_up)))) + sum(sum(cijt.*X_t0(1:i_up,1:j_up,1)))...
  + sum(sum(sum(pijt.*U(1:i_up,1:j_up,1:t_up)))) + sum(sum(pijt.*U_t0(1:i_up,1:j_up,1)))...
  + sum(sum(hit.*V(1:i_up,1:t_up))) + sum(hit.*V_t0(1:i_up,1));

%1.35-express D_ij_SAV
for i=1:i_up
  for j=1:j_up
    for t=1:t_up
      D_ij_SAV(i,j,t) = (b*d(i,j,t)+cijt*d(i,j,t)*m-n*T_abs*U(i,j,t))./(a+b+m*(cijt+P(i,j,t)));
    end
  end
end

% 0604-calculate minimum k ratio covered by SAVs
for i=1:i_up
  for j=1:j_up
    for t=1:t_up
      k(i,j,t) = D_ij_SAV(i,j,t)./d(i,j,t);
    end
  end
end

% k_min = min(k(1,1,1),k(1,2,1),k(1,3,1))

% k_min = 0.99;
% % 
% if ((k_min) >= (k(1,1,1)))
%   k_min = k(1,1,1)                
% end  

% 0604-find minimum covered ratio k
% for i=1:i_up
%   for j=1:j_up
%     for t=1:t_up
%       if (k_min >= k(i,j,t))
%         k_min = k(i,j,t)                
%       end  
%     end
%   end
% end

%1.4-check if solved and output
result = optimize(C,z_op,ops); % 0508-origin
% result = solvesdp(C,z_op,ops); % 0508-this works

%result = optimize(C,z);
% if result.problem == 0 % problem =0 means solve successfully,only print Uijt here

%   fprintf('\nk value: ') %ratio k covered by SAVs
%   value(k)

  
%   fprintf('\nP_Ini_(1,1,1) = ')
%   value(P(1,1,1))

%   fprintf('\nMin cost by Lower = ')
%   value(z_op)

  all_k = [];

  for i=1:i_up
    for j=1:j_up
      for t=1:t_up
        all_k(end+1) = k(i,j,t);
      end
    end
  end
  
  min_k = min(all_k);
  
  k_avg = mean(all_k);
  
%   fprintf('\n min k: ')
%   value(min_k)
  
  max_k = max(all_k);
  
%   fprintf('\n max k: ')
%   value(max_k)
  
  % 0615-save max_k in txt

%   fid = fopen(['IniP=',num2str(2),'max_k.txt'],'a'); %save to last 'a'
%   fprintf(fid,'%.4f\t',double(max_k)); %0516-keep 4 decimals
%   fclose(fid);
  
%   fprintf('\nU(1,1,1) value: ')
%   value(U(1,1,1))
        
%   fprintf('\nX(1,1,1) value: ') 
%   value(X(1,1,1))

fprintf('\n Iteration No. = ')
value(Iter_count)



% fprintf('\n D_ij_SAV(1,1,1) value: ') %0515-as for now hide D_ij_SAV
% value(D_ij_SAV(1,1,1))
  
%   fprintf('\n k(1,1,1) ')
%   value(k(1,1,1))

%   fprintf('\n')
%     fprintf('\nV value: ')
%     value(V)
%     
%     fprintf('\nU_t0 value: ')
%     value(U_t0)
%     
%     fprintf('\nX_t0 value: ') 
%     value(X_t0)
%     
%     fprintf('\nV_t0 value: ')
%     value(V_t0)

% else
%   disp('Not solved!')
% 
% end

% 0519-save total cost z_op in txt
% fid = fopen(['Iter=',num2str(Iter_count),'_OperCost.txt'],'wt'); %creat & write z_op to txt
% fprintf(fid,'%.4f\n',double(z_op)); %0516-keep 4 decimals
% fclose(fid);

% 0605-save min ratio k covered by SAVs in txt
% fid_ratiok = fopen(['Iter=',num2str(Iter_count),'_MinRatiok.txt'],'wt'); %creat & write z_op to txt
% fprintf(fid_ratiok,'%.4f\n',double(min_k)); %0604-keep 4 decimals
% fclose(fid_ratiok);

% 0517-following is save Uijt in txt
%1.6-save Uijt as txt for upper use
% for U_t=1:t_up 
%      fid = fopen(['IniP=',num2str(IniP),'Iter=',num2str(Iter_count),'_Uijt=',num2str(U_t),'.txt'],'wt'); %creat & write Uijt to txt
%      for i=1:i_up
%          for j=1:j_up
%              if j == j_up
%                  fprintf(fid,'%.4f\n',double(U(i,j,U_t))); %0516-keep 4 decimals
%              else
%                  fprintf(fid,'%.4f\t',double(U(i,j,U_t))); %0516-keep 4 decimals
%              end     
%          end    
%      end    
%      fclose(fid);
% end

% 0607-calculate Upper Pijt upper bound
for i=1:i_up
  for j=1:j_up
    for t=1:t_up
      P(i,j,t)=(b*d(i,j,t)+cijt*d(i,j,t)*m-n*T_abs*U(i,j,t)-(a+b)*min_k*d(i,j,t))./(m*min_k*d(i,j,t))-cijt;
    end
  end
end

% fprintf('\nUpper bound P(1,1,1) = %.4f\n',double(P(1,1,1)));

%0608-----
%calculate Upper revenue
%0608-----

revenue=0;

for i=1:i_up
  for j=1:j_up
    for t=1:t_up
      if (t==1)||(t==2)||(t==9)||(t==10)||(t==17)||(t==18)||(t==31)||(t==32)||(t==39)||(t==40)
%         T_abs=2;
        revenue = revenue-(P(i,j,t).*(b.*d(i,j,t)+cijt.*d(i,j,t)*m-n*T_abs*U(i,j,t)))./(a+b+m.*(cijt+P(i,j,t))); 
      else
%         T_abs=2; %here should be 3, but temporarily set to 2, aligh with Oper level
        revenue = revenue-(P(i,j,t).*(b.*d(i,j,t)+cijt.*d(i,j,t)*m-n*T_abs*U(i,j,t)))./(a+b+m.*(cijt+P(i,j,t)));
      end
    end
  end
end

% 070320-avg. all P

 all_P = [];

  for i=1:i_up
    for j=1:j_up
      for t=1:t_up
        all_P = P(i,j,t);
      end
    end
  end

avg_P = all_P./(i_up*j_up*t_up);



% 070320-avg. all D

 all_Dijt = [];

  for i=1:i_up
    for j=1:j_up
      for t=1:t_up
        all_Dijt = D_ij_SAV(i,j,t);
      end
    end
  end

avg_Dijt = all_Dijt./(i_up*j_up*t_up);

% 070320-print results

fprintf('\n Upper bound P(1,1,1) = %.4f\n',double(P(1,1,1))); % particular P each iteration

fprintf('\n');

fprintf('\n D(1,1,1) = %.4f\n',double(D_ij_SAV(1,1,1))); % particular Dijt each iteration

fprintf('\n');

fprintf('\n Avg_P = %.4f\n',double(avg_P)); % particular avg_P each iteration

fprintf('\n');

fprintf('\n Avg_Dijt = %.4f\n',double(avg_Dijt)); % particular avg_Dijt each iteration

fprintf('\n');

fprintf('\n Max Revenue by Upper= %.4f\n',double(-revenue));

fprintf('\n');

fprintf('\n Min Cost by Lower= %.4f\n',double(z_op));

fprintf('\n');

fprintf('\n Profit = %.4f\n',double(-revenue-z_op));

fprintf('\n');

% test how many new P > old P

satis_P = 0;

for i=1:i_up
  for j=1:j_up
    for t=1:t_up
%       if (double(P_sav(i,j,t))<=double(P(i,j,t))) % 0610-check new P greater than old ones
        %0610-if new P no more than 1 cent than old, than add to satis_P terminate condition
      %if ((double(P_sav(i,j,t)))*1.0001 >= double(P(i,j,t))) 
%       if ((double(P_sav(i,j,t)))+0.01 >= double(P(i,j,t))) 
      if ((P_sav(i,j,t))+0.01 > P(i,j,t)) 
        satis_P = (satis_P)+1;
      end  
    end
  end
end

satis_P;

fprintf('\n');

% test how many new P all feasible that new_P <= new P upper bound condition


if (Iter_count==1)
  feasib_P = 0;
else
  
  feasib_P = 0;
  
  for i=1:i_up
    for j=1:j_up
      for t=1:t_up
        P_bench(i,j,t) = (b*d(i,j,t)+cijt*d(i,j,t)*m-n*T_abs*U(i,j,t)-(a+b)*min_k*d(i,j,t))./(m*min_k*d(i,j,t))-cijt;        
      end
    end
  end
  
    
  for i=1:i_up
    for j=1:j_up
      for t=1:t_up
  %       if (double(P_sav(i,j,t))<=double(P(i,j,t))) % 0610-check new P greater than old ones
          %0610-if new P no more than 1 cent than old, than add to satis_P terminate condition
        %if ((double(P_sav(i,j,t)))*1.0001 >= double(P(i,j,t))) 
        
%         if (double(P_sav(i,j,t)) <= double(P_bench(i,j,t)))
        if (P_bench(i,j,t) >= P_sav(i,j,t))
%           i,j,t
%           fprintf("\n");
          feasib_P = (feasib_P)+1;
        end  
      end
    end
  end

end

feasib_P;

fprintf('\n');

if (satis_P==i_up*j_up*t_up) && (feasib_P==i_up*j_up*t_up)
  consec_P = 1+consec_P;
else
  consec_P = 0;
end

% tik = tik-1;

end

fprintf('\n');
fprintf('\n');
fprintf('\n');


%-------case 7.3------

'-------Case 7.3-----'

IniP=100;

%initially price Pijt all the same
P=zeros(20,20,40);%create 20*20*40 matrix
for i=1:i_up 
    for j=1:j_up
        for t=1:t_up
          P(i,j,t)=IniP; %0608-initial P should be infinity, I use 2 for test
        end
    end    
end

% Pijt=2;

Iter_count = 0;

% 0517-M is c16 time window, set at 1 gives it relaxation
% 0424-M is c16 time window, set at 3 gives it relaxation
M=1; 

%0515-change all var to continours instead of MIP
X = sdpvar(i_up,j_up,t_up+M); %0424-Xi,j,t+M, cuz time window constraint c16 Xijt is beyond t
%X = intvar(i_up,j_up,t_up+M); %0424-Xi,j,t+M, cuz time window constraint c16 Xijt is beyond t
X_t0 = sdpvar(i_up,j_up,1); %0424-Xi,j,0 means moves initially for c18: inventory balance eq
%X_t0 = intvar(i_up,j_up,1); %0424-Xi,j,0 means moves initially for c18: inventory balance eq
U = sdpvar(i_up,j_up,t_up);
%U = intvar(i_up,j_up,t_up); 
U_t0 = sdpvar(i_up,j_up,1); %0424-this is for Ui,j,0 
%U_t0 = intvar(i_up,j_up,1); %0424-this is for Ui,j,0 
V = sdpvar(i_up,t_up); 
%V = intvar(i_up,t_up); 
V_t0 = sdpvar(i_up,1); %0424-this is especially for Vi,0
%V_t0 = intvar(i_up,1); %0424-this is especially for Vi,0
D_ij_SAV = sdpvar(i_up,j_up,t_up); %0425-this is to show D_ij_SAV, for futher plot, then find upper bound of Pijt
%D_ij_SAV = intvar(i_up,j_up,t_up); %0425-this is to show D_ij_SAV, for futher plot, then find upper bound of Pijt
% 0604-k ratio
k = sdpvar(i_up,j_up,t_up);

% 070220-set place for Wijt
W=zeros(20,20,40);%create 20*20*40 matrix

%1.1-parameters
% a=b=m=n=1,T=40
a=1;
b=1;
m=1;
n=1;
T=8; %070120-total 40 time periods, 5 time block, so every block is 8
T_abs=8; %070120-T_abs should be same as T
% cijt=40, pijt=10, hit=10, alphaijt=0.1 (all c,p,h,alpha are subjectively to be same)
cijt=20;
pijt=10; % 070120-increase for penalty of unmet demand
% 0425-intentionally make pijt this small to avoid all 0 in Uijt solution
hit=30;
alpha=1;

%%%-----0607----
%%% How to add loop?
%%%----

% 0610-here starts the loop

P_sav = zeros(20,20,40);
U_sav = zeros(20,20,40);
min_k = 0;
P_bench = zeros(20,20,40);

% tik = 2;

satis_P = 0;
feasib_P = 0;
consec_P = 0;

% while (satis_P<i_up*j_up*t_up) || (feasib_P<i_up*j_up*t_up)
% while Iter_count < 2
while (consec_P<3)

Iter_count = Iter_count+1; 


for i=1:i_up 
    for j=1:j_up
        for t=1:t_up
          P_sav(i,j,t)=P(i,j,t); %0608-initial P should be infinity, I use 2 for test
        end
    end    
end  


for i=1:i_up 
    for j=1:j_up
        for t=1:t_up
            U_sav(i,j,t)=U(i,j,t); 
        end
    end    
end


min_k_sav = min_k;


%1.2-add operational level constraints

C = [];

% c15
for i=1:i_up %c15
   for j=1:j_up
     for t=1:t_up
       if (t==1)
         C = [C, U(i,j,1) >= U_t0(i,j,1)+(b*d(i,j,1)+cijt*d(i,j,1)*m-n*T_abs*U(i,j,1))./(a+b+m*(cijt+P(i,j,t)))-X(i,j,1)]; %c15: t=1
       else       
         C = [C, U(i,j,t) >= U(i,j,t-1)+(b*d(i,j,t)+cijt*d(i,j,t)*m-n*T_abs*U(i,j,t))./(a+b+m*(cijt+P(i,j,t)))-X(i,j,t)]; %c15: t from 2
       end
     end
   end    
end

% 0424-c16-I suppose time window is M period for c16 - relax this one, hope get solution
% this also considers Xijt beyond T
for i=1:i_up %c16
    for j=1:j_up
      for t=1:t_up
        if (t==1)
          C = [C, sum(X(i,j,t:(t+M))) >= U_t0(i,j,1)+(b*d(i,j,t)+cijt*d(i,j,t)*m-n*T_abs*U(i,j,t))./(a+b+m*(cijt+P(i,j,t)))]; %c16: t=1
        else
          C = [C, sum(X(i,j,t:(t+M))) >= U(i,j,t-1)+(b*d(i,j,t)+cijt*d(i,j,t)*m-n*T_abs*U(i,j,t))./(a+b+m*(cijt+P(i,j,t)))]; %c16: t from 2
        end
      end  
    end     
end

% 0424 - c17 is in upper level, no need in lower level
% for i=1:i_up %c17
%   for j=1:j_up
%     for t=1:t_up
%       C = [C, ((b-1)*d(i,j,t)+cijt*d(i,j,t)*m)./(n*T_abs) <= U(i,j,t)];
%       C = [C, U(i,j,t) <= (b*d(i,j,t)+cijt*d(i,j,t)*m)./(n*T_abs)];
%     end
%   end
% end

% 0424-c18-here makes <= for Vi,t cuz alpha all set to 1, too relaxed, should be variaty
for i=1:i_up %c18
  for t=1:t_up
    if (t==1)
      C = [C, V(i,t) <= V_t0(i,1) + sum(alpha*X_t0(1:j_up,i,1)) - sum(X(i,1:j_up,t)) ]; %c18: t=1
    else
      C = [C, V(i,t) <= V(i,t-1) + sum(sum(alpha*X(1:j_up,i,1:(t-1)))) - sum(X(i,1:j_up,t)) ]; %c18: t from 2 (remember to use 2 sum for 2nd)
    end
  end
end

% 0508-c19-not added, cuz it looks already relaxed in math function, hidden property of our formula

for i=1:i_up %c20: Xijt
    for j=1:j_up
        for t=1:(t_up+M)
          C = [C, X(i,j,t) >= 0];
        end
    end    
end    

for i=1:i_up %c20: X_t0(i,j,1) is Xijt=0
    for j=1:j_up
        C = [C, X_t0(i,j,1) >= 0];
    end    
end    

for i=1:i_up %c20: Uijt
    for j=1:j_up
        for t=1:t_up
            C = [C, U(i,j,t) >= 0];
        end
    end    
end    

for i=1:i_up %c20: U_t0(i,j,1) is Uijt=0
    for j=1:j_up
        C = [C, U_t0(i,j,1) >= 0];
    end    
end 

for i=1:i_up %c20: Vit
    for t=1:t_up
      C = [C, V(i,t) >= 0];
    end    
end    

for i=1:i_up %c20: V_t0(i,1) is V(i,0)
  C = [C, V_t0(i,1) >= 0];
end

%1.3-solve operational level: settings: use cplex as solver

%0515-following is set optimality for MIP, but kang said make it LNContinuous-https://www.ibm.com/support/knowledgecenter/SSSA5P_12.9.0/ilog.odms.cplex.help/CPLEX/Parameters/topics/OptimalityTarget.html
%0515-'cplex.optimalitytarget'=1: Searches for a globally optimal solution to a convex model.
ops = sdpsettings ('solver','cplex','verbose',2,'cplex.optimalitytarget',1);%0508-verbose 2 shows iterations, 1 means shown normal message, 0 shows silent message
%0515-MIP limit tolerance gap to 1-https://www.ibm.com/support/knowledgecenter/SSSA5P_12.9.0/ilog.odms.cplex.help/CPLEX/Parameters/topics/EpGap.html
%ops = sdpsettings ('solver','cplex','verbose',2,'cplex.mip.tolerances.mipgap',1);%0508-verbose 2 shows iterations, 1 means shown normal message, 0 shows silent message

%0306-limit barrier iteration-https://www.ibm.com/support/knowledgecenter/SSSA5P_12.9.0/ilog.odms.cplex.help/CPLEX/Parameters/topics/BarItLim.html
%Cplex.Param.barrier.limits.iteration = 4;[not work]

%0306-limit network simplex iteration-https://www.ibm.com/support/knowledgecenter/SSSA5P_12.9.0/ilog.odms.cplex.help/CPLEX/Parameters/topics/NetItLim.html
%Cplex.Param.network.iterations = 2;[not work]

%ops.mip.limits.nodes = 4;[not work]

%0304-use the online method for out-of-memory: https://www.ilovematlab.cn/thread-266244-1-1.html
%ops.mip.strategy.file = 3; %[not work: should put in sdpsetting] this has same results 1159110 node, but reduces time for 3 min.

%0305-use online method for out-of-memory: https://www.ibm.com/developerworks/community/forums/html/topic?id=e3ea344c-7249-4edd-a47d-29e1f80fa480
%Cplex.Param.mip.limits.auxrootthreads = 2; %[not work: should put in sdpsetting] -1 or 1

%obj-operational obj z_op
z_op = sum(sum(sum(cijt.*X(1:i_up,1:j_up,1:t_up)))) + sum(sum(cijt.*X_t0(1:i_up,1:j_up,1)))...
  + sum(sum(sum(pijt.*U(1:i_up,1:j_up,1:t_up)))) + sum(sum(pijt.*U_t0(1:i_up,1:j_up,1)))...
  + sum(sum(hit.*V(1:i_up,1:t_up))) + sum(hit.*V_t0(1:i_up,1));

%1.35-express D_ij_SAV
for i=1:i_up
  for j=1:j_up
    for t=1:t_up
      D_ij_SAV(i,j,t) = (b*d(i,j,t)+cijt*d(i,j,t)*m-n*T_abs*U(i,j,t))./(a+b+m*(cijt+P(i,j,t)));
    end
  end
end

% 0604-calculate minimum k ratio covered by SAVs
for i=1:i_up
  for j=1:j_up
    for t=1:t_up
      k(i,j,t) = D_ij_SAV(i,j,t)./d(i,j,t);
    end
  end
end

% k_min = min(k(1,1,1),k(1,2,1),k(1,3,1))

% k_min = 0.99;
% % 
% if ((k_min) >= (k(1,1,1)))
%   k_min = k(1,1,1)                
% end  

% 0604-find minimum covered ratio k
% for i=1:i_up
%   for j=1:j_up
%     for t=1:t_up
%       if (k_min >= k(i,j,t))
%         k_min = k(i,j,t)                
%       end  
%     end
%   end
% end

%1.4-check if solved and output
result = optimize(C,z_op,ops); % 0508-origin
% result = solvesdp(C,z_op,ops); % 0508-this works

%result = optimize(C,z);
% if result.problem == 0 % problem =0 means solve successfully,only print Uijt here

%   fprintf('\nk value: ') %ratio k covered by SAVs
%   value(k)

  
%   fprintf('\nP_Ini_(1,1,1) = ')
%   value(P(1,1,1))

%   fprintf('\nMin cost by Lower = ')
%   value(z_op)

  all_k = [];

  for i=1:i_up
    for j=1:j_up
      for t=1:t_up
        all_k(end+1) = k(i,j,t);
      end
    end
  end
  
  min_k = min(all_k);
  
  k_avg = mean(all_k);
  
%   fprintf('\n min k: ')
%   value(min_k)
  
  max_k = max(all_k);
  
%   fprintf('\n max k: ')
%   value(max_k)
  
  % 0615-save max_k in txt

%   fid = fopen(['IniP=',num2str(2),'max_k.txt'],'a'); %save to last 'a'
%   fprintf(fid,'%.4f\t',double(max_k)); %0516-keep 4 decimals
%   fclose(fid);
  
%   fprintf('\nU(1,1,1) value: ')
%   value(U(1,1,1))
        
%   fprintf('\nX(1,1,1) value: ') 
%   value(X(1,1,1))

fprintf('\n Iteration No. = ')
value(Iter_count)



% fprintf('\n D_ij_SAV(1,1,1) value: ') %0515-as for now hide D_ij_SAV
% value(D_ij_SAV(1,1,1))
  
%   fprintf('\n k(1,1,1) ')
%   value(k(1,1,1))

%   fprintf('\n')
%     fprintf('\nV value: ')
%     value(V)
%     
%     fprintf('\nU_t0 value: ')
%     value(U_t0)
%     
%     fprintf('\nX_t0 value: ') 
%     value(X_t0)
%     
%     fprintf('\nV_t0 value: ')
%     value(V_t0)

% else
%   disp('Not solved!')
% 
% end

% 0519-save total cost z_op in txt
% fid = fopen(['Iter=',num2str(Iter_count),'_OperCost.txt'],'wt'); %creat & write z_op to txt
% fprintf(fid,'%.4f\n',double(z_op)); %0516-keep 4 decimals
% fclose(fid);

% 0605-save min ratio k covered by SAVs in txt
% fid_ratiok = fopen(['Iter=',num2str(Iter_count),'_MinRatiok.txt'],'wt'); %creat & write z_op to txt
% fprintf(fid_ratiok,'%.4f\n',double(min_k)); %0604-keep 4 decimals
% fclose(fid_ratiok);

% 0517-following is save Uijt in txt
%1.6-save Uijt as txt for upper use
% for U_t=1:t_up 
%      fid = fopen(['IniP=',num2str(IniP),'Iter=',num2str(Iter_count),'_Uijt=',num2str(U_t),'.txt'],'wt'); %creat & write Uijt to txt
%      for i=1:i_up
%          for j=1:j_up
%              if j == j_up
%                  fprintf(fid,'%.4f\n',double(U(i,j,U_t))); %0516-keep 4 decimals
%              else
%                  fprintf(fid,'%.4f\t',double(U(i,j,U_t))); %0516-keep 4 decimals
%              end     
%          end    
%      end    
%      fclose(fid);
% end

% 0607-calculate Upper Pijt upper bound
for i=1:i_up
  for j=1:j_up
    for t=1:t_up
      P(i,j,t)=(b*d(i,j,t)+cijt*d(i,j,t)*m-n*T_abs*U(i,j,t)-(a+b)*min_k*d(i,j,t))./(m*min_k*d(i,j,t))-cijt;
    end
  end
end

% fprintf('\nUpper bound P(1,1,1) = %.4f\n',double(P(1,1,1)));

%0608-----
%calculate Upper revenue
%0608-----

revenue=0;

for i=1:i_up
  for j=1:j_up
    for t=1:t_up
      if (t==1)||(t==2)||(t==9)||(t==10)||(t==17)||(t==18)||(t==31)||(t==32)||(t==39)||(t==40)
%         T_abs=2;
        revenue = revenue-(P(i,j,t).*(b.*d(i,j,t)+cijt.*d(i,j,t)*m-n*T_abs*U(i,j,t)))./(a+b+m.*(cijt+P(i,j,t))); 
      else
%         T_abs=2; %here should be 3, but temporarily set to 2, aligh with Oper level
        revenue = revenue-(P(i,j,t).*(b.*d(i,j,t)+cijt.*d(i,j,t)*m-n*T_abs*U(i,j,t)))./(a+b+m.*(cijt+P(i,j,t)));
      end
    end
  end
end

% 070320-avg. all P

 all_P = [];

  for i=1:i_up
    for j=1:j_up
      for t=1:t_up
        all_P = P(i,j,t);
      end
    end
  end

avg_P = all_P./(i_up*j_up*t_up);



% 070320-avg. all D

 all_Dijt = [];

  for i=1:i_up
    for j=1:j_up
      for t=1:t_up
        all_Dijt = D_ij_SAV(i,j,t);
      end
    end
  end

avg_Dijt = all_Dijt./(i_up*j_up*t_up);

% 070320-print results

fprintf('\n Upper bound P(1,1,1) = %.4f\n',double(P(1,1,1))); % particular P each iteration

fprintf('\n');

fprintf('\n D(1,1,1) = %.4f\n',double(D_ij_SAV(1,1,1))); % particular Dijt each iteration

fprintf('\n');

fprintf('\n Avg_P = %.4f\n',double(avg_P)); % particular avg_P each iteration

fprintf('\n');

fprintf('\n Avg_Dijt = %.4f\n',double(avg_Dijt)); % particular avg_Dijt each iteration

fprintf('\n');

fprintf('\n Max Revenue by Upper= %.4f\n',double(-revenue));

fprintf('\n');

fprintf('\n Min Cost by Lower= %.4f\n',double(z_op));

fprintf('\n');

fprintf('\n Profit = %.4f\n',double(-revenue-z_op));

fprintf('\n');

% test how many new P > old P

satis_P = 0;

for i=1:i_up
  for j=1:j_up
    for t=1:t_up
%       if (double(P_sav(i,j,t))<=double(P(i,j,t))) % 0610-check new P greater than old ones
        %0610-if new P no more than 1 cent than old, than add to satis_P terminate condition
      %if ((double(P_sav(i,j,t)))*1.0001 >= double(P(i,j,t))) 
%       if ((double(P_sav(i,j,t)))+0.01 >= double(P(i,j,t))) 
      if ((P_sav(i,j,t))+0.01 > P(i,j,t)) 
        satis_P = (satis_P)+1;
      end  
    end
  end
end

satis_P;

fprintf('\n');

% test how many new P all feasible that new_P <= new P upper bound condition


if (Iter_count==1)
  feasib_P = 0;
else
  
  feasib_P = 0;
  
  for i=1:i_up
    for j=1:j_up
      for t=1:t_up
        P_bench(i,j,t) = (b*d(i,j,t)+cijt*d(i,j,t)*m-n*T_abs*U(i,j,t)-(a+b)*min_k*d(i,j,t))./(m*min_k*d(i,j,t))-cijt;        
      end
    end
  end
  
    
  for i=1:i_up
    for j=1:j_up
      for t=1:t_up
  %       if (double(P_sav(i,j,t))<=double(P(i,j,t))) % 0610-check new P greater than old ones
          %0610-if new P no more than 1 cent than old, than add to satis_P terminate condition
        %if ((double(P_sav(i,j,t)))*1.0001 >= double(P(i,j,t))) 
        
%         if (double(P_sav(i,j,t)) <= double(P_bench(i,j,t)))
        if (P_bench(i,j,t) >= P_sav(i,j,t))
%           i,j,t
%           fprintf("\n");
          feasib_P = (feasib_P)+1;
        end  
      end
    end
  end

end

feasib_P;

fprintf('\n');

if (satis_P==i_up*j_up*t_up) && (feasib_P==i_up*j_up*t_up)
  consec_P = 1+consec_P;
else
  consec_P = 0;
end

% tik = tik-1;

end

fprintf('\n');
fprintf('\n');
fprintf('\n');



%-------case 8.1------

'------Case 8.1------'

IniP=1;

%initially price Pijt all the same
P=zeros(20,20,40);%create 20*20*40 matrix
for i=1:i_up 
    for j=1:j_up
        for t=1:t_up
          P(i,j,t)=IniP; %0608-initial P should be infinity, I use 2 for test
        end
    end    
end

% Pijt=2;

Iter_count = 0;

% 0517-M is c16 time window, set at 1 gives it relaxation
% 0424-M is c16 time window, set at 3 gives it relaxation
M=1; 

%0515-change all var to continours instead of MIP
X = sdpvar(i_up,j_up,t_up+M); %0424-Xi,j,t+M, cuz time window constraint c16 Xijt is beyond t
%X = intvar(i_up,j_up,t_up+M); %0424-Xi,j,t+M, cuz time window constraint c16 Xijt is beyond t
X_t0 = sdpvar(i_up,j_up,1); %0424-Xi,j,0 means moves initially for c18: inventory balance eq
%X_t0 = intvar(i_up,j_up,1); %0424-Xi,j,0 means moves initially for c18: inventory balance eq
U = sdpvar(i_up,j_up,t_up);
%U = intvar(i_up,j_up,t_up); 
U_t0 = sdpvar(i_up,j_up,1); %0424-this is for Ui,j,0 
%U_t0 = intvar(i_up,j_up,1); %0424-this is for Ui,j,0 
V = sdpvar(i_up,t_up); 
%V = intvar(i_up,t_up); 
V_t0 = sdpvar(i_up,1); %0424-this is especially for Vi,0
%V_t0 = intvar(i_up,1); %0424-this is especially for Vi,0
D_ij_SAV = sdpvar(i_up,j_up,t_up); %0425-this is to show D_ij_SAV, for futher plot, then find upper bound of Pijt
%D_ij_SAV = intvar(i_up,j_up,t_up); %0425-this is to show D_ij_SAV, for futher plot, then find upper bound of Pijt
% 0604-k ratio
k = sdpvar(i_up,j_up,t_up);

% 070220-set place for Wijt
W=zeros(20,20,40);%create 20*20*40 matrix

%1.1-parameters
% a=b=m=n=1,T=40
a=1;
b=1;
m=1;
n=1;
T=8; %070120-total 40 time periods, 5 time block, so every block is 8
T_abs=8; %070120-T_abs should be same as T
% cijt=40, pijt=10, hit=10, alphaijt=0.1 (all c,p,h,alpha are subjectively to be same)
cijt=20;
pijt=10; % 070120-increase for penalty of unmet demand
% 0425-intentionally make pijt this small to avoid all 0 in Uijt solution
hit=10;
alpha=0.5;

%%%-----0607----
%%% How to add loop?
%%%----

% 0610-here starts the loop

P_sav = zeros(20,20,40);
U_sav = zeros(20,20,40);
min_k = 0;
P_bench = zeros(20,20,40);

% tik = 2;

satis_P = 0;
feasib_P = 0;
consec_P = 0;

% while (satis_P<i_up*j_up*t_up) || (feasib_P<i_up*j_up*t_up)
% while Iter_count < 2
while (consec_P<3)

Iter_count = Iter_count+1; 


for i=1:i_up 
    for j=1:j_up
        for t=1:t_up
          P_sav(i,j,t)=P(i,j,t); %0608-initial P should be infinity, I use 2 for test
        end
    end    
end  


for i=1:i_up 
    for j=1:j_up
        for t=1:t_up
            U_sav(i,j,t)=U(i,j,t); 
        end
    end    
end


min_k_sav = min_k;


%1.2-add operational level constraints

C = [];

% c15
for i=1:i_up %c15
   for j=1:j_up
     for t=1:t_up
       if (t==1)
         C = [C, U(i,j,1) >= U_t0(i,j,1)+(b*d(i,j,1)+cijt*d(i,j,1)*m-n*T_abs*U(i,j,1))./(a+b+m*(cijt+P(i,j,t)))-X(i,j,1)]; %c15: t=1
       else       
         C = [C, U(i,j,t) >= U(i,j,t-1)+(b*d(i,j,t)+cijt*d(i,j,t)*m-n*T_abs*U(i,j,t))./(a+b+m*(cijt+P(i,j,t)))-X(i,j,t)]; %c15: t from 2
       end
     end
   end    
end

% 0424-c16-I suppose time window is M period for c16 - relax this one, hope get solution
% this also considers Xijt beyond T
for i=1:i_up %c16
    for j=1:j_up
      for t=1:t_up
        if (t==1)
          C = [C, sum(X(i,j,t:(t+M))) >= U_t0(i,j,1)+(b*d(i,j,t)+cijt*d(i,j,t)*m-n*T_abs*U(i,j,t))./(a+b+m*(cijt+P(i,j,t)))]; %c16: t=1
        else
          C = [C, sum(X(i,j,t:(t+M))) >= U(i,j,t-1)+(b*d(i,j,t)+cijt*d(i,j,t)*m-n*T_abs*U(i,j,t))./(a+b+m*(cijt+P(i,j,t)))]; %c16: t from 2
        end
      end  
    end     
end

% 0424 - c17 is in upper level, no need in lower level
% for i=1:i_up %c17
%   for j=1:j_up
%     for t=1:t_up
%       C = [C, ((b-1)*d(i,j,t)+cijt*d(i,j,t)*m)./(n*T_abs) <= U(i,j,t)];
%       C = [C, U(i,j,t) <= (b*d(i,j,t)+cijt*d(i,j,t)*m)./(n*T_abs)];
%     end
%   end
% end

% 0424-c18-here makes <= for Vi,t cuz alpha all set to 1, too relaxed, should be variaty
for i=1:i_up %c18
  for t=1:t_up
    if (t==1)
      C = [C, V(i,t) <= V_t0(i,1) + sum(alpha*X_t0(1:j_up,i,1)) - sum(X(i,1:j_up,t)) ]; %c18: t=1
    else
      C = [C, V(i,t) <= V(i,t-1) + sum(sum(alpha*X(1:j_up,i,1:(t-1)))) - sum(X(i,1:j_up,t)) ]; %c18: t from 2 (remember to use 2 sum for 2nd)
    end
  end
end

% 0508-c19-not added, cuz it looks already relaxed in math function, hidden property of our formula

for i=1:i_up %c20: Xijt
    for j=1:j_up
        for t=1:(t_up+M)
          C = [C, X(i,j,t) >= 0];
        end
    end    
end    

for i=1:i_up %c20: X_t0(i,j,1) is Xijt=0
    for j=1:j_up
        C = [C, X_t0(i,j,1) >= 0];
    end    
end    

for i=1:i_up %c20: Uijt
    for j=1:j_up
        for t=1:t_up
            C = [C, U(i,j,t) >= 0];
        end
    end    
end    

for i=1:i_up %c20: U_t0(i,j,1) is Uijt=0
    for j=1:j_up
        C = [C, U_t0(i,j,1) >= 0];
    end    
end 

for i=1:i_up %c20: Vit
    for t=1:t_up
      C = [C, V(i,t) >= 0];
    end    
end    

for i=1:i_up %c20: V_t0(i,1) is V(i,0)
  C = [C, V_t0(i,1) >= 0];
end

%1.3-solve operational level: settings: use cplex as solver

%0515-following is set optimality for MIP, but kang said make it LNContinuous-https://www.ibm.com/support/knowledgecenter/SSSA5P_12.9.0/ilog.odms.cplex.help/CPLEX/Parameters/topics/OptimalityTarget.html
%0515-'cplex.optimalitytarget'=1: Searches for a globally optimal solution to a convex model.
ops = sdpsettings ('solver','cplex','verbose',2,'cplex.optimalitytarget',1);%0508-verbose 2 shows iterations, 1 means shown normal message, 0 shows silent message
%0515-MIP limit tolerance gap to 1-https://www.ibm.com/support/knowledgecenter/SSSA5P_12.9.0/ilog.odms.cplex.help/CPLEX/Parameters/topics/EpGap.html
%ops = sdpsettings ('solver','cplex','verbose',2,'cplex.mip.tolerances.mipgap',1);%0508-verbose 2 shows iterations, 1 means shown normal message, 0 shows silent message

%0306-limit barrier iteration-https://www.ibm.com/support/knowledgecenter/SSSA5P_12.9.0/ilog.odms.cplex.help/CPLEX/Parameters/topics/BarItLim.html
%Cplex.Param.barrier.limits.iteration = 4;[not work]

%0306-limit network simplex iteration-https://www.ibm.com/support/knowledgecenter/SSSA5P_12.9.0/ilog.odms.cplex.help/CPLEX/Parameters/topics/NetItLim.html
%Cplex.Param.network.iterations = 2;[not work]

%ops.mip.limits.nodes = 4;[not work]

%0304-use the online method for out-of-memory: https://www.ilovematlab.cn/thread-266244-1-1.html
%ops.mip.strategy.file = 3; %[not work: should put in sdpsetting] this has same results 1159110 node, but reduces time for 3 min.

%0305-use online method for out-of-memory: https://www.ibm.com/developerworks/community/forums/html/topic?id=e3ea344c-7249-4edd-a47d-29e1f80fa480
%Cplex.Param.mip.limits.auxrootthreads = 2; %[not work: should put in sdpsetting] -1 or 1

%obj-operational obj z_op
z_op = sum(sum(sum(cijt.*X(1:i_up,1:j_up,1:t_up)))) + sum(sum(cijt.*X_t0(1:i_up,1:j_up,1)))...
  + sum(sum(sum(pijt.*U(1:i_up,1:j_up,1:t_up)))) + sum(sum(pijt.*U_t0(1:i_up,1:j_up,1)))...
  + sum(sum(hit.*V(1:i_up,1:t_up))) + sum(hit.*V_t0(1:i_up,1));

%1.35-express D_ij_SAV
for i=1:i_up
  for j=1:j_up
    for t=1:t_up
      D_ij_SAV(i,j,t) = (b*d(i,j,t)+cijt*d(i,j,t)*m-n*T_abs*U(i,j,t))./(a+b+m*(cijt+P(i,j,t)));
    end
  end
end

% 0604-calculate minimum k ratio covered by SAVs
for i=1:i_up
  for j=1:j_up
    for t=1:t_up
      k(i,j,t) = D_ij_SAV(i,j,t)./d(i,j,t);
    end
  end
end

% k_min = min(k(1,1,1),k(1,2,1),k(1,3,1))

% k_min = 0.99;
% % 
% if ((k_min) >= (k(1,1,1)))
%   k_min = k(1,1,1)                
% end  

% 0604-find minimum covered ratio k
% for i=1:i_up
%   for j=1:j_up
%     for t=1:t_up
%       if (k_min >= k(i,j,t))
%         k_min = k(i,j,t)                
%       end  
%     end
%   end
% end

%1.4-check if solved and output
result = optimize(C,z_op,ops); % 0508-origin
% result = solvesdp(C,z_op,ops); % 0508-this works

%result = optimize(C,z);
% if result.problem == 0 % problem =0 means solve successfully,only print Uijt here

%   fprintf('\nk value: ') %ratio k covered by SAVs
%   value(k)

  
%   fprintf('\nP_Ini_(1,1,1) = ')
%   value(P(1,1,1))

%   fprintf('\nMin cost by Lower = ')
%   value(z_op)

  all_k = [];

  for i=1:i_up
    for j=1:j_up
      for t=1:t_up
        all_k(end+1) = k(i,j,t);
      end
    end
  end
  
  min_k = min(all_k);
  
  k_avg = mean(all_k);
  
%   fprintf('\n min k: ')
%   value(min_k)
  
  max_k = max(all_k);
  
%   fprintf('\n max k: ')
%   value(max_k)
  
  % 0615-save max_k in txt

%   fid = fopen(['IniP=',num2str(2),'max_k.txt'],'a'); %save to last 'a'
%   fprintf(fid,'%.4f\t',double(max_k)); %0516-keep 4 decimals
%   fclose(fid);
  
%   fprintf('\nU(1,1,1) value: ')
%   value(U(1,1,1))
        
%   fprintf('\nX(1,1,1) value: ') 
%   value(X(1,1,1))

fprintf('\n Iteration No. = ')
value(Iter_count)



% fprintf('\n D_ij_SAV(1,1,1) value: ') %0515-as for now hide D_ij_SAV
% value(D_ij_SAV(1,1,1))
  
%   fprintf('\n k(1,1,1) ')
%   value(k(1,1,1))

%   fprintf('\n')
%     fprintf('\nV value: ')
%     value(V)
%     
%     fprintf('\nU_t0 value: ')
%     value(U_t0)
%     
%     fprintf('\nX_t0 value: ') 
%     value(X_t0)
%     
%     fprintf('\nV_t0 value: ')
%     value(V_t0)

% else
%   disp('Not solved!')
% 
% end

% 0519-save total cost z_op in txt
% fid = fopen(['Iter=',num2str(Iter_count),'_OperCost.txt'],'wt'); %creat & write z_op to txt
% fprintf(fid,'%.4f\n',double(z_op)); %0516-keep 4 decimals
% fclose(fid);

% 0605-save min ratio k covered by SAVs in txt
% fid_ratiok = fopen(['Iter=',num2str(Iter_count),'_MinRatiok.txt'],'wt'); %creat & write z_op to txt
% fprintf(fid_ratiok,'%.4f\n',double(min_k)); %0604-keep 4 decimals
% fclose(fid_ratiok);

% 0517-following is save Uijt in txt
%1.6-save Uijt as txt for upper use
% for U_t=1:t_up 
%      fid = fopen(['IniP=',num2str(IniP),'Iter=',num2str(Iter_count),'_Uijt=',num2str(U_t),'.txt'],'wt'); %creat & write Uijt to txt
%      for i=1:i_up
%          for j=1:j_up
%              if j == j_up
%                  fprintf(fid,'%.4f\n',double(U(i,j,U_t))); %0516-keep 4 decimals
%              else
%                  fprintf(fid,'%.4f\t',double(U(i,j,U_t))); %0516-keep 4 decimals
%              end     
%          end    
%      end    
%      fclose(fid);
% end

% 0607-calculate Upper Pijt upper bound
for i=1:i_up
  for j=1:j_up
    for t=1:t_up
      P(i,j,t)=(b*d(i,j,t)+cijt*d(i,j,t)*m-n*T_abs*U(i,j,t)-(a+b)*min_k*d(i,j,t))./(m*min_k*d(i,j,t))-cijt;
    end
  end
end

% fprintf('\nUpper bound P(1,1,1) = %.4f\n',double(P(1,1,1)));

%0608-----
%calculate Upper revenue
%0608-----

revenue=0;

for i=1:i_up
  for j=1:j_up
    for t=1:t_up
      if (t==1)||(t==2)||(t==9)||(t==10)||(t==17)||(t==18)||(t==31)||(t==32)||(t==39)||(t==40)
%         T_abs=2;
        revenue = revenue-(P(i,j,t).*(b.*d(i,j,t)+cijt.*d(i,j,t)*m-n*T_abs*U(i,j,t)))./(a+b+m.*(cijt+P(i,j,t))); 
      else
%         T_abs=2; %here should be 3, but temporarily set to 2, aligh with Oper level
        revenue = revenue-(P(i,j,t).*(b.*d(i,j,t)+cijt.*d(i,j,t)*m-n*T_abs*U(i,j,t)))./(a+b+m.*(cijt+P(i,j,t)));
      end
    end
  end
end

% 070320-avg. all P

 all_P = [];

  for i=1:i_up
    for j=1:j_up
      for t=1:t_up
        all_P = P(i,j,t);
      end
    end
  end

avg_P = all_P./(i_up*j_up*t_up);



% 070320-avg. all D

 all_Dijt = [];

  for i=1:i_up
    for j=1:j_up
      for t=1:t_up
        all_Dijt = D_ij_SAV(i,j,t);
      end
    end
  end

avg_Dijt = all_Dijt./(i_up*j_up*t_up);

% 070320-print results

fprintf('\n Upper bound P(1,1,1) = %.4f\n',double(P(1,1,1))); % particular P each iteration

fprintf('\n');

fprintf('\n D(1,1,1) = %.4f\n',double(D_ij_SAV(1,1,1))); % particular Dijt each iteration

fprintf('\n');

fprintf('\n Avg_P = %.4f\n',double(avg_P)); % particular avg_P each iteration

fprintf('\n');

fprintf('\n Avg_Dijt = %.4f\n',double(avg_Dijt)); % particular avg_Dijt each iteration

fprintf('\n');

fprintf('\n Max Revenue by Upper= %.4f\n',double(-revenue));

fprintf('\n');

fprintf('\n Min Cost by Lower= %.4f\n',double(z_op));

fprintf('\n');

fprintf('\n Profit = %.4f\n',double(-revenue-z_op));

fprintf('\n');

% test how many new P > old P

satis_P = 0;

for i=1:i_up
  for j=1:j_up
    for t=1:t_up
%       if (double(P_sav(i,j,t))<=double(P(i,j,t))) % 0610-check new P greater than old ones
        %0610-if new P no more than 1 cent than old, than add to satis_P terminate condition
      %if ((double(P_sav(i,j,t)))*1.0001 >= double(P(i,j,t))) 
%       if ((double(P_sav(i,j,t)))+0.01 >= double(P(i,j,t))) 
      if ((P_sav(i,j,t))+0.01 > P(i,j,t)) 
        satis_P = (satis_P)+1;
      end  
    end
  end
end

satis_P;

fprintf('\n');

% test how many new P all feasible that new_P <= new P upper bound condition


if (Iter_count==1)
  feasib_P = 0;
else
  
  feasib_P = 0;
  
  for i=1:i_up
    for j=1:j_up
      for t=1:t_up
        P_bench(i,j,t) = (b*d(i,j,t)+cijt*d(i,j,t)*m-n*T_abs*U(i,j,t)-(a+b)*min_k*d(i,j,t))./(m*min_k*d(i,j,t))-cijt;        
      end
    end
  end
  
    
  for i=1:i_up
    for j=1:j_up
      for t=1:t_up
  %       if (double(P_sav(i,j,t))<=double(P(i,j,t))) % 0610-check new P greater than old ones
          %0610-if new P no more than 1 cent than old, than add to satis_P terminate condition
        %if ((double(P_sav(i,j,t)))*1.0001 >= double(P(i,j,t))) 
        
%         if (double(P_sav(i,j,t)) <= double(P_bench(i,j,t)))
        if (P_bench(i,j,t) >= P_sav(i,j,t))
%           i,j,t
%           fprintf("\n");
          feasib_P = (feasib_P)+1;
        end  
      end
    end
  end

end

feasib_P;

fprintf('\n');

if (satis_P==i_up*j_up*t_up) && (feasib_P==i_up*j_up*t_up)
  consec_P = 1+consec_P;
else
  consec_P = 0;
end

% tik = tik-1;

end

fprintf('\n');
fprintf('\n');
fprintf('\n');


%-------case 8.2------

'-------Case 8.2-------'

IniP=10;

%initially price Pijt all the same
P=zeros(20,20,40);%create 20*20*40 matrix
for i=1:i_up 
    for j=1:j_up
        for t=1:t_up
          P(i,j,t)=IniP; %0608-initial P should be infinity, I use 2 for test
        end
    end    
end

% Pijt=2;

Iter_count = 0;

% 0517-M is c16 time window, set at 1 gives it relaxation
% 0424-M is c16 time window, set at 3 gives it relaxation
M=1; 

%0515-change all var to continours instead of MIP
X = sdpvar(i_up,j_up,t_up+M); %0424-Xi,j,t+M, cuz time window constraint c16 Xijt is beyond t
%X = intvar(i_up,j_up,t_up+M); %0424-Xi,j,t+M, cuz time window constraint c16 Xijt is beyond t
X_t0 = sdpvar(i_up,j_up,1); %0424-Xi,j,0 means moves initially for c18: inventory balance eq
%X_t0 = intvar(i_up,j_up,1); %0424-Xi,j,0 means moves initially for c18: inventory balance eq
U = sdpvar(i_up,j_up,t_up);
%U = intvar(i_up,j_up,t_up); 
U_t0 = sdpvar(i_up,j_up,1); %0424-this is for Ui,j,0 
%U_t0 = intvar(i_up,j_up,1); %0424-this is for Ui,j,0 
V = sdpvar(i_up,t_up); 
%V = intvar(i_up,t_up); 
V_t0 = sdpvar(i_up,1); %0424-this is especially for Vi,0
%V_t0 = intvar(i_up,1); %0424-this is especially for Vi,0
D_ij_SAV = sdpvar(i_up,j_up,t_up); %0425-this is to show D_ij_SAV, for futher plot, then find upper bound of Pijt
%D_ij_SAV = intvar(i_up,j_up,t_up); %0425-this is to show D_ij_SAV, for futher plot, then find upper bound of Pijt
% 0604-k ratio
k = sdpvar(i_up,j_up,t_up);

% 070220-set place for Wijt
W=zeros(20,20,40);%create 20*20*40 matrix

%1.1-parameters
% a=b=m=n=1,T=40
a=1;
b=1;
m=1;
n=1;
T=8; %070120-total 40 time periods, 5 time block, so every block is 8
T_abs=8; %070120-T_abs should be same as T
% cijt=40, pijt=10, hit=10, alphaijt=0.1 (all c,p,h,alpha are subjectively to be same)
cijt=20;
pijt=10; % 070120-increase for penalty of unmet demand
% 0425-intentionally make pijt this small to avoid all 0 in Uijt solution
hit=10;
alpha=0.5;

%%%-----0607----
%%% How to add loop?
%%%----

% 0610-here starts the loop

P_sav = zeros(20,20,40);
U_sav = zeros(20,20,40);
min_k = 0;
P_bench = zeros(20,20,40);

% tik = 2;

satis_P = 0;
feasib_P = 0;
consec_P = 0;

% while (satis_P<i_up*j_up*t_up) || (feasib_P<i_up*j_up*t_up)
% while Iter_count < 2
while (consec_P<3)

Iter_count = Iter_count+1; 


for i=1:i_up 
    for j=1:j_up
        for t=1:t_up
          P_sav(i,j,t)=P(i,j,t); %0608-initial P should be infinity, I use 2 for test
        end
    end    
end  


for i=1:i_up 
    for j=1:j_up
        for t=1:t_up
            U_sav(i,j,t)=U(i,j,t); 
        end
    end    
end


min_k_sav = min_k;


%1.2-add operational level constraints

C = [];

% c15
for i=1:i_up %c15
   for j=1:j_up
     for t=1:t_up
       if (t==1)
         C = [C, U(i,j,1) >= U_t0(i,j,1)+(b*d(i,j,1)+cijt*d(i,j,1)*m-n*T_abs*U(i,j,1))./(a+b+m*(cijt+P(i,j,t)))-X(i,j,1)]; %c15: t=1
       else       
         C = [C, U(i,j,t) >= U(i,j,t-1)+(b*d(i,j,t)+cijt*d(i,j,t)*m-n*T_abs*U(i,j,t))./(a+b+m*(cijt+P(i,j,t)))-X(i,j,t)]; %c15: t from 2
       end
     end
   end    
end

% 0424-c16-I suppose time window is M period for c16 - relax this one, hope get solution
% this also considers Xijt beyond T
for i=1:i_up %c16
    for j=1:j_up
      for t=1:t_up
        if (t==1)
          C = [C, sum(X(i,j,t:(t+M))) >= U_t0(i,j,1)+(b*d(i,j,t)+cijt*d(i,j,t)*m-n*T_abs*U(i,j,t))./(a+b+m*(cijt+P(i,j,t)))]; %c16: t=1
        else
          C = [C, sum(X(i,j,t:(t+M))) >= U(i,j,t-1)+(b*d(i,j,t)+cijt*d(i,j,t)*m-n*T_abs*U(i,j,t))./(a+b+m*(cijt+P(i,j,t)))]; %c16: t from 2
        end
      end  
    end     
end

% 0424 - c17 is in upper level, no need in lower level
% for i=1:i_up %c17
%   for j=1:j_up
%     for t=1:t_up
%       C = [C, ((b-1)*d(i,j,t)+cijt*d(i,j,t)*m)./(n*T_abs) <= U(i,j,t)];
%       C = [C, U(i,j,t) <= (b*d(i,j,t)+cijt*d(i,j,t)*m)./(n*T_abs)];
%     end
%   end
% end

% 0424-c18-here makes <= for Vi,t cuz alpha all set to 1, too relaxed, should be variaty
for i=1:i_up %c18
  for t=1:t_up
    if (t==1)
      C = [C, V(i,t) <= V_t0(i,1) + sum(alpha*X_t0(1:j_up,i,1)) - sum(X(i,1:j_up,t)) ]; %c18: t=1
    else
      C = [C, V(i,t) <= V(i,t-1) + sum(sum(alpha*X(1:j_up,i,1:(t-1)))) - sum(X(i,1:j_up,t)) ]; %c18: t from 2 (remember to use 2 sum for 2nd)
    end
  end
end

% 0508-c19-not added, cuz it looks already relaxed in math function, hidden property of our formula

for i=1:i_up %c20: Xijt
    for j=1:j_up
        for t=1:(t_up+M)
          C = [C, X(i,j,t) >= 0];
        end
    end    
end    

for i=1:i_up %c20: X_t0(i,j,1) is Xijt=0
    for j=1:j_up
        C = [C, X_t0(i,j,1) >= 0];
    end    
end    

for i=1:i_up %c20: Uijt
    for j=1:j_up
        for t=1:t_up
            C = [C, U(i,j,t) >= 0];
        end
    end    
end    

for i=1:i_up %c20: U_t0(i,j,1) is Uijt=0
    for j=1:j_up
        C = [C, U_t0(i,j,1) >= 0];
    end    
end 

for i=1:i_up %c20: Vit
    for t=1:t_up
      C = [C, V(i,t) >= 0];
    end    
end    

for i=1:i_up %c20: V_t0(i,1) is V(i,0)
  C = [C, V_t0(i,1) >= 0];
end

%1.3-solve operational level: settings: use cplex as solver

%0515-following is set optimality for MIP, but kang said make it LNContinuous-https://www.ibm.com/support/knowledgecenter/SSSA5P_12.9.0/ilog.odms.cplex.help/CPLEX/Parameters/topics/OptimalityTarget.html
%0515-'cplex.optimalitytarget'=1: Searches for a globally optimal solution to a convex model.
ops = sdpsettings ('solver','cplex','verbose',2,'cplex.optimalitytarget',1);%0508-verbose 2 shows iterations, 1 means shown normal message, 0 shows silent message
%0515-MIP limit tolerance gap to 1-https://www.ibm.com/support/knowledgecenter/SSSA5P_12.9.0/ilog.odms.cplex.help/CPLEX/Parameters/topics/EpGap.html
%ops = sdpsettings ('solver','cplex','verbose',2,'cplex.mip.tolerances.mipgap',1);%0508-verbose 2 shows iterations, 1 means shown normal message, 0 shows silent message

%0306-limit barrier iteration-https://www.ibm.com/support/knowledgecenter/SSSA5P_12.9.0/ilog.odms.cplex.help/CPLEX/Parameters/topics/BarItLim.html
%Cplex.Param.barrier.limits.iteration = 4;[not work]

%0306-limit network simplex iteration-https://www.ibm.com/support/knowledgecenter/SSSA5P_12.9.0/ilog.odms.cplex.help/CPLEX/Parameters/topics/NetItLim.html
%Cplex.Param.network.iterations = 2;[not work]

%ops.mip.limits.nodes = 4;[not work]

%0304-use the online method for out-of-memory: https://www.ilovematlab.cn/thread-266244-1-1.html
%ops.mip.strategy.file = 3; %[not work: should put in sdpsetting] this has same results 1159110 node, but reduces time for 3 min.

%0305-use online method for out-of-memory: https://www.ibm.com/developerworks/community/forums/html/topic?id=e3ea344c-7249-4edd-a47d-29e1f80fa480
%Cplex.Param.mip.limits.auxrootthreads = 2; %[not work: should put in sdpsetting] -1 or 1

%obj-operational obj z_op
z_op = sum(sum(sum(cijt.*X(1:i_up,1:j_up,1:t_up)))) + sum(sum(cijt.*X_t0(1:i_up,1:j_up,1)))...
  + sum(sum(sum(pijt.*U(1:i_up,1:j_up,1:t_up)))) + sum(sum(pijt.*U_t0(1:i_up,1:j_up,1)))...
  + sum(sum(hit.*V(1:i_up,1:t_up))) + sum(hit.*V_t0(1:i_up,1));

%1.35-express D_ij_SAV
for i=1:i_up
  for j=1:j_up
    for t=1:t_up
      D_ij_SAV(i,j,t) = (b*d(i,j,t)+cijt*d(i,j,t)*m-n*T_abs*U(i,j,t))./(a+b+m*(cijt+P(i,j,t)));
    end
  end
end

% 0604-calculate minimum k ratio covered by SAVs
for i=1:i_up
  for j=1:j_up
    for t=1:t_up
      k(i,j,t) = D_ij_SAV(i,j,t)./d(i,j,t);
    end
  end
end

% k_min = min(k(1,1,1),k(1,2,1),k(1,3,1))

% k_min = 0.99;
% % 
% if ((k_min) >= (k(1,1,1)))
%   k_min = k(1,1,1)                
% end  

% 0604-find minimum covered ratio k
% for i=1:i_up
%   for j=1:j_up
%     for t=1:t_up
%       if (k_min >= k(i,j,t))
%         k_min = k(i,j,t)                
%       end  
%     end
%   end
% end

%1.4-check if solved and output
result = optimize(C,z_op,ops); % 0508-origin
% result = solvesdp(C,z_op,ops); % 0508-this works

%result = optimize(C,z);
% if result.problem == 0 % problem =0 means solve successfully,only print Uijt here

%   fprintf('\nk value: ') %ratio k covered by SAVs
%   value(k)

  
%   fprintf('\nP_Ini_(1,1,1) = ')
%   value(P(1,1,1))

%   fprintf('\nMin cost by Lower = ')
%   value(z_op)

  all_k = [];

  for i=1:i_up
    for j=1:j_up
      for t=1:t_up
        all_k(end+1) = k(i,j,t);
      end
    end
  end
  
  min_k = min(all_k);
  
  k_avg = mean(all_k);
  
%   fprintf('\n min k: ')
%   value(min_k)
  
  max_k = max(all_k);
  
%   fprintf('\n max k: ')
%   value(max_k)
  
  % 0615-save max_k in txt

%   fid = fopen(['IniP=',num2str(2),'max_k.txt'],'a'); %save to last 'a'
%   fprintf(fid,'%.4f\t',double(max_k)); %0516-keep 4 decimals
%   fclose(fid);
  
%   fprintf('\nU(1,1,1) value: ')
%   value(U(1,1,1))
        
%   fprintf('\nX(1,1,1) value: ') 
%   value(X(1,1,1))

fprintf('\n Iteration No. = ')
value(Iter_count)



% fprintf('\n D_ij_SAV(1,1,1) value: ') %0515-as for now hide D_ij_SAV
% value(D_ij_SAV(1,1,1))
  
%   fprintf('\n k(1,1,1) ')
%   value(k(1,1,1))

%   fprintf('\n')
%     fprintf('\nV value: ')
%     value(V)
%     
%     fprintf('\nU_t0 value: ')
%     value(U_t0)
%     
%     fprintf('\nX_t0 value: ') 
%     value(X_t0)
%     
%     fprintf('\nV_t0 value: ')
%     value(V_t0)

% else
%   disp('Not solved!')
% 
% end

% 0519-save total cost z_op in txt
% fid = fopen(['Iter=',num2str(Iter_count),'_OperCost.txt'],'wt'); %creat & write z_op to txt
% fprintf(fid,'%.4f\n',double(z_op)); %0516-keep 4 decimals
% fclose(fid);

% 0605-save min ratio k covered by SAVs in txt
% fid_ratiok = fopen(['Iter=',num2str(Iter_count),'_MinRatiok.txt'],'wt'); %creat & write z_op to txt
% fprintf(fid_ratiok,'%.4f\n',double(min_k)); %0604-keep 4 decimals
% fclose(fid_ratiok);

% 0517-following is save Uijt in txt
%1.6-save Uijt as txt for upper use
% for U_t=1:t_up 
%      fid = fopen(['IniP=',num2str(IniP),'Iter=',num2str(Iter_count),'_Uijt=',num2str(U_t),'.txt'],'wt'); %creat & write Uijt to txt
%      for i=1:i_up
%          for j=1:j_up
%              if j == j_up
%                  fprintf(fid,'%.4f\n',double(U(i,j,U_t))); %0516-keep 4 decimals
%              else
%                  fprintf(fid,'%.4f\t',double(U(i,j,U_t))); %0516-keep 4 decimals
%              end     
%          end    
%      end    
%      fclose(fid);
% end

% 0607-calculate Upper Pijt upper bound
for i=1:i_up
  for j=1:j_up
    for t=1:t_up
      P(i,j,t)=(b*d(i,j,t)+cijt*d(i,j,t)*m-n*T_abs*U(i,j,t)-(a+b)*min_k*d(i,j,t))./(m*min_k*d(i,j,t))-cijt;
    end
  end
end

% fprintf('\nUpper bound P(1,1,1) = %.4f\n',double(P(1,1,1)));

%0608-----
%calculate Upper revenue
%0608-----

revenue=0;

for i=1:i_up
  for j=1:j_up
    for t=1:t_up
      if (t==1)||(t==2)||(t==9)||(t==10)||(t==17)||(t==18)||(t==31)||(t==32)||(t==39)||(t==40)
%         T_abs=2;
        revenue = revenue-(P(i,j,t).*(b.*d(i,j,t)+cijt.*d(i,j,t)*m-n*T_abs*U(i,j,t)))./(a+b+m.*(cijt+P(i,j,t))); 
      else
%         T_abs=2; %here should be 3, but temporarily set to 2, aligh with Oper level
        revenue = revenue-(P(i,j,t).*(b.*d(i,j,t)+cijt.*d(i,j,t)*m-n*T_abs*U(i,j,t)))./(a+b+m.*(cijt+P(i,j,t)));
      end
    end
  end
end

% 070320-avg. all P

 all_P = [];

  for i=1:i_up
    for j=1:j_up
      for t=1:t_up
        all_P = P(i,j,t);
      end
    end
  end

avg_P = all_P./(i_up*j_up*t_up);



% 070320-avg. all D

 all_Dijt = [];

  for i=1:i_up
    for j=1:j_up
      for t=1:t_up
        all_Dijt = D_ij_SAV(i,j,t);
      end
    end
  end

avg_Dijt = all_Dijt./(i_up*j_up*t_up);

% 070320-print results

fprintf('\n Upper bound P(1,1,1) = %.4f\n',double(P(1,1,1))); % particular P each iteration

fprintf('\n');

fprintf('\n D(1,1,1) = %.4f\n',double(D_ij_SAV(1,1,1))); % particular Dijt each iteration

fprintf('\n');

fprintf('\n Avg_P = %.4f\n',double(avg_P)); % particular avg_P each iteration

fprintf('\n');

fprintf('\n Avg_Dijt = %.4f\n',double(avg_Dijt)); % particular avg_Dijt each iteration

fprintf('\n');

fprintf('\n Max Revenue by Upper= %.4f\n',double(-revenue));

fprintf('\n');

fprintf('\n Min Cost by Lower= %.4f\n',double(z_op));

fprintf('\n');

fprintf('\n Profit = %.4f\n',double(-revenue-z_op));

fprintf('\n');

% test how many new P > old P

satis_P = 0;

for i=1:i_up
  for j=1:j_up
    for t=1:t_up
%       if (double(P_sav(i,j,t))<=double(P(i,j,t))) % 0610-check new P greater than old ones
        %0610-if new P no more than 1 cent than old, than add to satis_P terminate condition
      %if ((double(P_sav(i,j,t)))*1.0001 >= double(P(i,j,t))) 
%       if ((double(P_sav(i,j,t)))+0.01 >= double(P(i,j,t))) 
      if ((P_sav(i,j,t))+0.01 > P(i,j,t)) 
        satis_P = (satis_P)+1;
      end  
    end
  end
end

satis_P;

fprintf('\n');

% test how many new P all feasible that new_P <= new P upper bound condition


if (Iter_count==1)
  feasib_P = 0;
else
  
  feasib_P = 0;
  
  for i=1:i_up
    for j=1:j_up
      for t=1:t_up
        P_bench(i,j,t) = (b*d(i,j,t)+cijt*d(i,j,t)*m-n*T_abs*U(i,j,t)-(a+b)*min_k*d(i,j,t))./(m*min_k*d(i,j,t))-cijt;        
      end
    end
  end
  
    
  for i=1:i_up
    for j=1:j_up
      for t=1:t_up
  %       if (double(P_sav(i,j,t))<=double(P(i,j,t))) % 0610-check new P greater than old ones
          %0610-if new P no more than 1 cent than old, than add to satis_P terminate condition
        %if ((double(P_sav(i,j,t)))*1.0001 >= double(P(i,j,t))) 
        
%         if (double(P_sav(i,j,t)) <= double(P_bench(i,j,t)))
        if (P_bench(i,j,t) >= P_sav(i,j,t))
%           i,j,t
%           fprintf("\n");
          feasib_P = (feasib_P)+1;
        end  
      end
    end
  end

end

feasib_P;

fprintf('\n');

if (satis_P==i_up*j_up*t_up) && (feasib_P==i_up*j_up*t_up)
  consec_P = 1+consec_P;
else
  consec_P = 0;
end

% tik = tik-1;

end

fprintf('\n');
fprintf('\n');
fprintf('\n');


%-------case 8.3------

'-------Case 8.3-----'

IniP=100;

%initially price Pijt all the same
P=zeros(20,20,40);%create 20*20*40 matrix
for i=1:i_up 
    for j=1:j_up
        for t=1:t_up
          P(i,j,t)=IniP; %0608-initial P should be infinity, I use 2 for test
        end
    end    
end

% Pijt=2;

Iter_count = 0;

% 0517-M is c16 time window, set at 1 gives it relaxation
% 0424-M is c16 time window, set at 3 gives it relaxation
M=1; 

%0515-change all var to continours instead of MIP
X = sdpvar(i_up,j_up,t_up+M); %0424-Xi,j,t+M, cuz time window constraint c16 Xijt is beyond t
%X = intvar(i_up,j_up,t_up+M); %0424-Xi,j,t+M, cuz time window constraint c16 Xijt is beyond t
X_t0 = sdpvar(i_up,j_up,1); %0424-Xi,j,0 means moves initially for c18: inventory balance eq
%X_t0 = intvar(i_up,j_up,1); %0424-Xi,j,0 means moves initially for c18: inventory balance eq
U = sdpvar(i_up,j_up,t_up);
%U = intvar(i_up,j_up,t_up); 
U_t0 = sdpvar(i_up,j_up,1); %0424-this is for Ui,j,0 
%U_t0 = intvar(i_up,j_up,1); %0424-this is for Ui,j,0 
V = sdpvar(i_up,t_up); 
%V = intvar(i_up,t_up); 
V_t0 = sdpvar(i_up,1); %0424-this is especially for Vi,0
%V_t0 = intvar(i_up,1); %0424-this is especially for Vi,0
D_ij_SAV = sdpvar(i_up,j_up,t_up); %0425-this is to show D_ij_SAV, for futher plot, then find upper bound of Pijt
%D_ij_SAV = intvar(i_up,j_up,t_up); %0425-this is to show D_ij_SAV, for futher plot, then find upper bound of Pijt
% 0604-k ratio
k = sdpvar(i_up,j_up,t_up);

% 070220-set place for Wijt
W=zeros(20,20,40);%create 20*20*40 matrix

%1.1-parameters
% a=b=m=n=1,T=40
a=1;
b=1;
m=1;
n=1;
T=8; %070120-total 40 time periods, 5 time block, so every block is 8
T_abs=8; %070120-T_abs should be same as T
% cijt=40, pijt=10, hit=10, alphaijt=0.1 (all c,p,h,alpha are subjectively to be same)
cijt=20;
pijt=10; % 070120-increase for penalty of unmet demand
% 0425-intentionally make pijt this small to avoid all 0 in Uijt solution
hit=10;
alpha=0.5;

%%%-----0607----
%%% How to add loop?
%%%----

% 0610-here starts the loop

P_sav = zeros(20,20,40);
U_sav = zeros(20,20,40);
min_k = 0;
P_bench = zeros(20,20,40);

% tik = 2;

satis_P = 0;
feasib_P = 0;
consec_P = 0;

% while (satis_P<i_up*j_up*t_up) || (feasib_P<i_up*j_up*t_up)
% while Iter_count < 2
while (consec_P<3)

Iter_count = Iter_count+1; 


for i=1:i_up 
    for j=1:j_up
        for t=1:t_up
          P_sav(i,j,t)=P(i,j,t); %0608-initial P should be infinity, I use 2 for test
        end
    end    
end  


for i=1:i_up 
    for j=1:j_up
        for t=1:t_up
            U_sav(i,j,t)=U(i,j,t); 
        end
    end    
end


min_k_sav = min_k;


%1.2-add operational level constraints

C = [];

% c15
for i=1:i_up %c15
   for j=1:j_up
     for t=1:t_up
       if (t==1)
         C = [C, U(i,j,1) >= U_t0(i,j,1)+(b*d(i,j,1)+cijt*d(i,j,1)*m-n*T_abs*U(i,j,1))./(a+b+m*(cijt+P(i,j,t)))-X(i,j,1)]; %c15: t=1
       else       
         C = [C, U(i,j,t) >= U(i,j,t-1)+(b*d(i,j,t)+cijt*d(i,j,t)*m-n*T_abs*U(i,j,t))./(a+b+m*(cijt+P(i,j,t)))-X(i,j,t)]; %c15: t from 2
       end
     end
   end    
end

% 0424-c16-I suppose time window is M period for c16 - relax this one, hope get solution
% this also considers Xijt beyond T
for i=1:i_up %c16
    for j=1:j_up
      for t=1:t_up
        if (t==1)
          C = [C, sum(X(i,j,t:(t+M))) >= U_t0(i,j,1)+(b*d(i,j,t)+cijt*d(i,j,t)*m-n*T_abs*U(i,j,t))./(a+b+m*(cijt+P(i,j,t)))]; %c16: t=1
        else
          C = [C, sum(X(i,j,t:(t+M))) >= U(i,j,t-1)+(b*d(i,j,t)+cijt*d(i,j,t)*m-n*T_abs*U(i,j,t))./(a+b+m*(cijt+P(i,j,t)))]; %c16: t from 2
        end
      end  
    end     
end

% 0424 - c17 is in upper level, no need in lower level
% for i=1:i_up %c17
%   for j=1:j_up
%     for t=1:t_up
%       C = [C, ((b-1)*d(i,j,t)+cijt*d(i,j,t)*m)./(n*T_abs) <= U(i,j,t)];
%       C = [C, U(i,j,t) <= (b*d(i,j,t)+cijt*d(i,j,t)*m)./(n*T_abs)];
%     end
%   end
% end

% 0424-c18-here makes <= for Vi,t cuz alpha all set to 1, too relaxed, should be variaty
for i=1:i_up %c18
  for t=1:t_up
    if (t==1)
      C = [C, V(i,t) <= V_t0(i,1) + sum(alpha*X_t0(1:j_up,i,1)) - sum(X(i,1:j_up,t)) ]; %c18: t=1
    else
      C = [C, V(i,t) <= V(i,t-1) + sum(sum(alpha*X(1:j_up,i,1:(t-1)))) - sum(X(i,1:j_up,t)) ]; %c18: t from 2 (remember to use 2 sum for 2nd)
    end
  end
end

% 0508-c19-not added, cuz it looks already relaxed in math function, hidden property of our formula

for i=1:i_up %c20: Xijt
    for j=1:j_up
        for t=1:(t_up+M)
          C = [C, X(i,j,t) >= 0];
        end
    end    
end    

for i=1:i_up %c20: X_t0(i,j,1) is Xijt=0
    for j=1:j_up
        C = [C, X_t0(i,j,1) >= 0];
    end    
end    

for i=1:i_up %c20: Uijt
    for j=1:j_up
        for t=1:t_up
            C = [C, U(i,j,t) >= 0];
        end
    end    
end    

for i=1:i_up %c20: U_t0(i,j,1) is Uijt=0
    for j=1:j_up
        C = [C, U_t0(i,j,1) >= 0];
    end    
end 

for i=1:i_up %c20: Vit
    for t=1:t_up
      C = [C, V(i,t) >= 0];
    end    
end    

for i=1:i_up %c20: V_t0(i,1) is V(i,0)
  C = [C, V_t0(i,1) >= 0];
end

%1.3-solve operational level: settings: use cplex as solver

%0515-following is set optimality for MIP, but kang said make it LNContinuous-https://www.ibm.com/support/knowledgecenter/SSSA5P_12.9.0/ilog.odms.cplex.help/CPLEX/Parameters/topics/OptimalityTarget.html
%0515-'cplex.optimalitytarget'=1: Searches for a globally optimal solution to a convex model.
ops = sdpsettings ('solver','cplex','verbose',2,'cplex.optimalitytarget',1);%0508-verbose 2 shows iterations, 1 means shown normal message, 0 shows silent message
%0515-MIP limit tolerance gap to 1-https://www.ibm.com/support/knowledgecenter/SSSA5P_12.9.0/ilog.odms.cplex.help/CPLEX/Parameters/topics/EpGap.html
%ops = sdpsettings ('solver','cplex','verbose',2,'cplex.mip.tolerances.mipgap',1);%0508-verbose 2 shows iterations, 1 means shown normal message, 0 shows silent message

%0306-limit barrier iteration-https://www.ibm.com/support/knowledgecenter/SSSA5P_12.9.0/ilog.odms.cplex.help/CPLEX/Parameters/topics/BarItLim.html
%Cplex.Param.barrier.limits.iteration = 4;[not work]

%0306-limit network simplex iteration-https://www.ibm.com/support/knowledgecenter/SSSA5P_12.9.0/ilog.odms.cplex.help/CPLEX/Parameters/topics/NetItLim.html
%Cplex.Param.network.iterations = 2;[not work]

%ops.mip.limits.nodes = 4;[not work]

%0304-use the online method for out-of-memory: https://www.ilovematlab.cn/thread-266244-1-1.html
%ops.mip.strategy.file = 3; %[not work: should put in sdpsetting] this has same results 1159110 node, but reduces time for 3 min.

%0305-use online method for out-of-memory: https://www.ibm.com/developerworks/community/forums/html/topic?id=e3ea344c-7249-4edd-a47d-29e1f80fa480
%Cplex.Param.mip.limits.auxrootthreads = 2; %[not work: should put in sdpsetting] -1 or 1

%obj-operational obj z_op
z_op = sum(sum(sum(cijt.*X(1:i_up,1:j_up,1:t_up)))) + sum(sum(cijt.*X_t0(1:i_up,1:j_up,1)))...
  + sum(sum(sum(pijt.*U(1:i_up,1:j_up,1:t_up)))) + sum(sum(pijt.*U_t0(1:i_up,1:j_up,1)))...
  + sum(sum(hit.*V(1:i_up,1:t_up))) + sum(hit.*V_t0(1:i_up,1));

%1.35-express D_ij_SAV
for i=1:i_up
  for j=1:j_up
    for t=1:t_up
      D_ij_SAV(i,j,t) = (b*d(i,j,t)+cijt*d(i,j,t)*m-n*T_abs*U(i,j,t))./(a+b+m*(cijt+P(i,j,t)));
    end
  end
end

% 0604-calculate minimum k ratio covered by SAVs
for i=1:i_up
  for j=1:j_up
    for t=1:t_up
      k(i,j,t) = D_ij_SAV(i,j,t)./d(i,j,t);
    end
  end
end

% k_min = min(k(1,1,1),k(1,2,1),k(1,3,1))

% k_min = 0.99;
% % 
% if ((k_min) >= (k(1,1,1)))
%   k_min = k(1,1,1)                
% end  

% 0604-find minimum covered ratio k
% for i=1:i_up
%   for j=1:j_up
%     for t=1:t_up
%       if (k_min >= k(i,j,t))
%         k_min = k(i,j,t)                
%       end  
%     end
%   end
% end

%1.4-check if solved and output
result = optimize(C,z_op,ops); % 0508-origin
% result = solvesdp(C,z_op,ops); % 0508-this works

%result = optimize(C,z);
% if result.problem == 0 % problem =0 means solve successfully,only print Uijt here

%   fprintf('\nk value: ') %ratio k covered by SAVs
%   value(k)

  
%   fprintf('\nP_Ini_(1,1,1) = ')
%   value(P(1,1,1))

%   fprintf('\nMin cost by Lower = ')
%   value(z_op)

  all_k = [];

  for i=1:i_up
    for j=1:j_up
      for t=1:t_up
        all_k(end+1) = k(i,j,t);
      end
    end
  end
  
  min_k = min(all_k);
  
  k_avg = mean(all_k);
  
%   fprintf('\n min k: ')
%   value(min_k)
  
  max_k = max(all_k);
  
%   fprintf('\n max k: ')
%   value(max_k)
  
  % 0615-save max_k in txt

%   fid = fopen(['IniP=',num2str(2),'max_k.txt'],'a'); %save to last 'a'
%   fprintf(fid,'%.4f\t',double(max_k)); %0516-keep 4 decimals
%   fclose(fid);
  
%   fprintf('\nU(1,1,1) value: ')
%   value(U(1,1,1))
        
%   fprintf('\nX(1,1,1) value: ') 
%   value(X(1,1,1))

fprintf('\n Iteration No. = ')
value(Iter_count)



% fprintf('\n D_ij_SAV(1,1,1) value: ') %0515-as for now hide D_ij_SAV
% value(D_ij_SAV(1,1,1))
  
%   fprintf('\n k(1,1,1) ')
%   value(k(1,1,1))

%   fprintf('\n')
%     fprintf('\nV value: ')
%     value(V)
%     
%     fprintf('\nU_t0 value: ')
%     value(U_t0)
%     
%     fprintf('\nX_t0 value: ') 
%     value(X_t0)
%     
%     fprintf('\nV_t0 value: ')
%     value(V_t0)

% else
%   disp('Not solved!')
% 
% end

% 0519-save total cost z_op in txt
% fid = fopen(['Iter=',num2str(Iter_count),'_OperCost.txt'],'wt'); %creat & write z_op to txt
% fprintf(fid,'%.4f\n',double(z_op)); %0516-keep 4 decimals
% fclose(fid);

% 0605-save min ratio k covered by SAVs in txt
% fid_ratiok = fopen(['Iter=',num2str(Iter_count),'_MinRatiok.txt'],'wt'); %creat & write z_op to txt
% fprintf(fid_ratiok,'%.4f\n',double(min_k)); %0604-keep 4 decimals
% fclose(fid_ratiok);

% 0517-following is save Uijt in txt
%1.6-save Uijt as txt for upper use
% for U_t=1:t_up 
%      fid = fopen(['IniP=',num2str(IniP),'Iter=',num2str(Iter_count),'_Uijt=',num2str(U_t),'.txt'],'wt'); %creat & write Uijt to txt
%      for i=1:i_up
%          for j=1:j_up
%              if j == j_up
%                  fprintf(fid,'%.4f\n',double(U(i,j,U_t))); %0516-keep 4 decimals
%              else
%                  fprintf(fid,'%.4f\t',double(U(i,j,U_t))); %0516-keep 4 decimals
%              end     
%          end    
%      end    
%      fclose(fid);
% end

% 0607-calculate Upper Pijt upper bound
for i=1:i_up
  for j=1:j_up
    for t=1:t_up
      P(i,j,t)=(b*d(i,j,t)+cijt*d(i,j,t)*m-n*T_abs*U(i,j,t)-(a+b)*min_k*d(i,j,t))./(m*min_k*d(i,j,t))-cijt;
    end
  end
end

% fprintf('\nUpper bound P(1,1,1) = %.4f\n',double(P(1,1,1)));

%0608-----
%calculate Upper revenue
%0608-----

revenue=0;

for i=1:i_up
  for j=1:j_up
    for t=1:t_up
      if (t==1)||(t==2)||(t==9)||(t==10)||(t==17)||(t==18)||(t==31)||(t==32)||(t==39)||(t==40)
%         T_abs=2;
        revenue = revenue-(P(i,j,t).*(b.*d(i,j,t)+cijt.*d(i,j,t)*m-n*T_abs*U(i,j,t)))./(a+b+m.*(cijt+P(i,j,t))); 
      else
%         T_abs=2; %here should be 3, but temporarily set to 2, aligh with Oper level
        revenue = revenue-(P(i,j,t).*(b.*d(i,j,t)+cijt.*d(i,j,t)*m-n*T_abs*U(i,j,t)))./(a+b+m.*(cijt+P(i,j,t)));
      end
    end
  end
end

% 070320-avg. all P

 all_P = [];

  for i=1:i_up
    for j=1:j_up
      for t=1:t_up
        all_P = P(i,j,t);
      end
    end
  end

avg_P = all_P./(i_up*j_up*t_up);



% 070320-avg. all D

 all_Dijt = [];

  for i=1:i_up
    for j=1:j_up
      for t=1:t_up
        all_Dijt = D_ij_SAV(i,j,t);
      end
    end
  end

avg_Dijt = all_Dijt./(i_up*j_up*t_up);

% 070320-print results

fprintf('\n Upper bound P(1,1,1) = %.4f\n',double(P(1,1,1))); % particular P each iteration

fprintf('\n');

fprintf('\n D(1,1,1) = %.4f\n',double(D_ij_SAV(1,1,1))); % particular Dijt each iteration

fprintf('\n');

fprintf('\n Avg_P = %.4f\n',double(avg_P)); % particular avg_P each iteration

fprintf('\n');

fprintf('\n Avg_Dijt = %.4f\n',double(avg_Dijt)); % particular avg_Dijt each iteration

fprintf('\n');

fprintf('\n Max Revenue by Upper= %.4f\n',double(-revenue));

fprintf('\n');

fprintf('\n Min Cost by Lower= %.4f\n',double(z_op));

fprintf('\n');

fprintf('\n Profit = %.4f\n',double(-revenue-z_op));

fprintf('\n');

% test how many new P > old P

satis_P = 0;

for i=1:i_up
  for j=1:j_up
    for t=1:t_up
%       if (double(P_sav(i,j,t))<=double(P(i,j,t))) % 0610-check new P greater than old ones
        %0610-if new P no more than 1 cent than old, than add to satis_P terminate condition
      %if ((double(P_sav(i,j,t)))*1.0001 >= double(P(i,j,t))) 
%       if ((double(P_sav(i,j,t)))+0.01 >= double(P(i,j,t))) 
      if ((P_sav(i,j,t))+0.01 > P(i,j,t)) 
        satis_P = (satis_P)+1;
      end  
    end
  end
end

satis_P;

fprintf('\n');

% test how many new P all feasible that new_P <= new P upper bound condition


if (Iter_count==1)
  feasib_P = 0;
else
  
  feasib_P = 0;
  
  for i=1:i_up
    for j=1:j_up
      for t=1:t_up
        P_bench(i,j,t) = (b*d(i,j,t)+cijt*d(i,j,t)*m-n*T_abs*U(i,j,t)-(a+b)*min_k*d(i,j,t))./(m*min_k*d(i,j,t))-cijt;        
      end
    end
  end
  
    
  for i=1:i_up
    for j=1:j_up
      for t=1:t_up
  %       if (double(P_sav(i,j,t))<=double(P(i,j,t))) % 0610-check new P greater than old ones
          %0610-if new P no more than 1 cent than old, than add to satis_P terminate condition
        %if ((double(P_sav(i,j,t)))*1.0001 >= double(P(i,j,t))) 
        
%         if (double(P_sav(i,j,t)) <= double(P_bench(i,j,t)))
        if (P_bench(i,j,t) >= P_sav(i,j,t))
%           i,j,t
%           fprintf("\n");
          feasib_P = (feasib_P)+1;
        end  
      end
    end
  end

end

feasib_P;

fprintf('\n');

if (satis_P==i_up*j_up*t_up) && (feasib_P==i_up*j_up*t_up)
  consec_P = 1+consec_P;
else
  consec_P = 0;
end

% tik = tik-1;

end

fprintf('\n');
fprintf('\n');
fprintf('\n');
