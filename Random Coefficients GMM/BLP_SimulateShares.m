function [sim_shares] = simShareCalc(delta_guess,data,sigma_guess,v)
%Takes in the original data matrix, a guess at delta and a guess of sigma
%calculates the simulated shares using these guesses
%the v is pulled each time, but the seed is fixed so the values will not
%change

len = length(data);
len_v = length(v);

%Initialize vectors to be used below
sim_shares = nan(len,1);
num = nan(len_v,1);
denom = nan(len_v,1);
ratio=nan(len_v,1);

%Basic idea is to isolate the 5 unique cases 
for i=1:len
    if data(i,2) == 1 && data(i,3) ==1
        
        for j=1:len_v
        num(j) = exp(delta_guess(i)+data(i,6)*sigma_guess*v(j));
        denom(j) = 1 + exp(delta_guess(i)+data(i,6)*sigma_guess*v(j)) + exp(delta_guess(i+1)+data(i+1,6)*sigma_guess*v(j));
        ratio(j) = num(j)/denom(j);
        end
        
        sim_shares(i) = mean(ratio);
        
    elseif data(i,2) == 1 && data(i,3) == 2
        
        for j=1:len_v
        num(j) = exp(delta_guess(i)+data(i,6)*sigma_guess*v(j));
        denom(j) = 1 + exp(delta_guess(i-1)+data(i-1,6)*sigma_guess*v(j)) + exp(delta_guess(i)+data(i,6)*sigma_guess*v(j));
        ratio(j) = num(j)/denom(j);
        end
        
        sim_shares(i) = mean(ratio);
        
    elseif data(i,2) == 2 && data(i,3) == 1
        
        for j=1:len_v
        num(j) = exp(delta_guess(i)+data(i,6)*sigma_guess*v(j));
        denom(j) = 1 + exp(delta_guess(i)+data(i,6)*sigma_guess*v(j)) + exp(delta_guess(i+1)+data(i+1,6)*sigma_guess*v(j)) + exp(delta_guess(i+2)+data(i+2,6)*sigma_guess*v(j)); 
        ratio(j) = num(j)/denom(j);
        end
        
        sim_shares(i) = mean(ratio);        
        
    elseif data(i,2) == 2 && data(i,3) == 2
        
        for j=1:len_v
        num(j) = exp(delta_guess(i)+data(i,6)*sigma_guess*v(j));
        denom(j) = 1 + exp(delta_guess(i-1)+data(i-1,6)*sigma_guess*v(j)) + exp(delta_guess(i)+data(i,6)*sigma_guess*v(j)) + exp(delta_guess(i+1)+data(i+1,6)*sigma_guess*v(j)); 
        ratio(j) = num(j)/denom(j);
        end
        
        sim_shares(i) = mean(ratio);
        
    elseif data(i,2) == 2 && data(i,3) == 3
        
        for j=1:len_v
        num(j) = exp(delta_guess(i)+data(i,6)*sigma_guess*v(j));
        denom(j) = 1 + exp(delta_guess(i-2)+data(i-2,6)*sigma_guess*v(j)) + exp(delta_guess(i-1)+data(i-1,6)*sigma_guess*v(j)) + exp(delta_guess(i)+data(i,6)*sigma_guess*v(j)); 
        ratio(j) = num(j)/denom(j);
        end
        
        sim_shares(i) = mean(ratio);
    
    else
        
    end    
end
