function lamda=DP_HMM(x)   %x is the (1*64) weight parameter vector, this is different from the output of loglikelihood as the last three components are added
global U;
global v;
global N_mode;

if length(x)~=4*(N_mode^2+2*N_mode)+1+N_mode   
    'Warning: dimension of input of loglikelihood is wrong'
    a %used to terminate the program
end

for i=1:N_mode
    for j=1:N_mode
        w{i,j}=x(1,N_mode*4*i+4*j-4*N_mode-3:4*N_mode*i+4*j-4*N_mode);
    end
end

beta=x(4*(N_mode^2+2*N_mode)+1);

for i=1:N_mode
    w_miu(i,1:4)=x(1,4*i+4*N_mode^2-3:4*i+4*N_mode^2);
    w_sigma(i,1:4)=x(1,4*i+4*N_mode^2+4*N_mode-3:4*i+4*N_mode^2+4*N_mode);
end


for N_e=1:length(v)   % # of events

N=length(v{N_e});


%the following parameters have to be given 
%U=zeros(N,4);   %N equals the total number of time instances, the 4th column corresponds to constant 1
%w=cell(3,3);   %each cell(i,j) store the weight of the transition matrix component M(i,j)
%v=zeros(N,1);  %v is the observed velocity
%beta=1.5;
%w_miu=zeros(3,4);   %weight for miu
%w_sigma=zeros(3,4);  %weight for sigma

global_log=zeros(N_mode,N);

for l=N:-1:1
    if l>=2
        V=beta*(v{N_e}(l)-v{N_e}(l-1))+v{N_e}(l-1);
    else
        V=v{N_e}(l);   %l=1, V=v(1)
    end
for k=1:N_mode
    miu(k,1)=w_miu(k,:)*U{N_e}(l,:).';
    sigma(k,1)=exp(w_sigma(k,:)*U{N_e}(l,:).');
end
for k=1:N_mode
    local_log(k,l)=-0.5*(V-miu(k))^2/sigma(k,1)-0.5*log(sigma(k,1));
end

if l<N   %for intermediate steps, calculate transition_cost
    for k=1:N_mode    %from state k to state j 
    for j=1:N_mode
    a(j,k)=w{j,k}*U{N_e}(l,:).';
    end
%     for j=1:3
%     tran_log(j,k)=log(exp(a(j,k))./sum(exp(a(:,k))));
%     end

     tran_log(:,k)=log(softmax(a(:,k)));

    end
end

    if l==N
        global_log(:,N)=local_log(:,N);
    end
    
    if l<N
        for k=1:N_mode
            [max_value,max_index]=max(tran_log(:,k)+global_log(:,l+1));
            global_log(k,l)=local_log(k,l)+max_value;
            path_index(k,l)=max_index;   %this records that at the l-th time step, if the optimal path goes through k-th state, then at l+1 time step, the optimal path is state path_index(k,l)
        end
    end
end

%add the initial cost to the global cost
for k=1:N_mode
    global_log(k,1)=global_log(k,1)+log(x(k+4*(N_mode^2+2*N_mode)+1));
end

[global_max,initial_index]=max(global_log(:,1));  %initial_index is the first state of the optimal path
optimal_path(1,1)=initial_index;   %trace back to find the optimal path
for k=1:N-1
    optimal_path(1,k+1)=path_index(optimal_path(1,k),k);
end

lamda{N_e}=optimal_path.';

for k=1:length(v{N_e})    %prescibe the abnormal velocities
if v{N_e}(k)<0
    lamda{N_e}(k)=1;
end
if v{N_e}(k)>4
    lamda{N_e}(k)=N_mode;
end
end

clearvars -except lamda U v x w beta w_miu w_sigma N_e N_mode
end   %end for N_e

end
    

    

