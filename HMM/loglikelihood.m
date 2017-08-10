function [p,g]=loglikelihood(x)
global U;
global lamda;   %here the lamda is the output of DP, make it a global variable shared with the main program
global v;
global N_mode;   %specify the number of state of the latent variable

if length(x)~=4*(N_mode^2+2*N_mode)+1     
    'Warning: dimension of input of loglikelihood is wrong'
    a %used to terminate the program
end

grad=cell(N_mode,N_mode);   %initialize gradient

for i=1:N_mode
    for j=1:N_mode
        w{i,j}=x(1,N_mode*4*i+4*j-4*N_mode-3:4*N_mode*i+4*j-4*N_mode);
        grad{i,j}=zeros(1,4);
    end
end

beta=x(4*(N_mode^2+2*N_mode)+1);

for i=1:N_mode
    w_miu(i,1:4)=x(1,4*i+4*N_mode^2-3:4*i+4*N_mode^2);
    w_sigma(i,1:4)=x(1,4*i+4*N_mode^2+4*N_mode-3:4*i+4*N_mode^2+4*N_mode);
end

p=0;  %initialize the loglikelihood value

grad_beta=0;
grad_wmiu=zeros(N_mode,4);
grad_wsigma=zeros(N_mode,4);

for N_e=1:length(v)   % # of events

N=length(v{N_e});

%w=cell(3,3);   %each cell(i,j) store the weight of the transition matrix component M(i,j)
%U=zeros(N,4);   %N equals the total number of time instances, the 4th column corresponds to constant 1
%lamda=zeros(N,1);   %lamda is the latent variable value at each time instiances
%v=zeros(N,1);  %v is the observed velocity
%w_miu=zeros(3,4);   %weight for miu
%w_sigma=zeros(3,4);  %weight for sigma
%here above are the input variables to calculate the gradients for the
%parameters

%calculate gradient
%MM=[];

for m=2:N   %may need to change
%     if m<=1
%     p=p+log(x(lamda{N_e}(1)+61));
%     else
    k=lamda{N_e}(m-1);
    j=lamda{N_e}(m);
    for l=1:N_mode
        a(l,1)=w{l,k}*U{N_e}(m,:).';   %calculate a(jk) where j~1:3, k=landa(m-1)
    end
%     for l=1:3
%         M(l,1)=exp(a(l))./sum(exp(a));   %calculate M(l,k), where l~1:3, k=landa(m-1)
%     end
      M(:,1)=softmax(a(:,1));  %use the softmax function instead of directly calculate it to avoid numerical issue
      
    for i=1:N_mode
        if i==j
        grad{i,k}=grad{i,k}+U{N_e}(m,:)*(1-M(i,1));
        else
        grad{i,k}=grad{i,k}-U{N_e}(m,:)*M(i,1);    %grad is row vector while M and a are temporary column vector.
        end
    end
    p=p+log(M(j,1));
%     end   %end if m<=1
%MM=[MM,log(M(:,1))];

end   %finish batch update of the gradient matrix for the transition prob

%then calculate gradient for beta,w_miu,w_sigma
%first initialize the gradient

for m=1:N
    
    for k=1:N_mode
    miu(k,1)=w_miu(k,:)*U{N_e}(m,:).';
    sigma(k,1)=exp(w_sigma(k,:)*U{N_e}(m,:).');
    end
    j=lamda{N_e}(m);
    
    if m>=2
    V=(v{N_e}(m)-v{N_e}(m-1))*beta+v{N_e}(m-1);
    else
    V=v{N_e}(m);   %if m=1, use the velocity at the current step
    end

    p=p-(V-miu(j))^2/2/sigma(j)-0.5*log(sigma(j));  %update the value of the loglikelihood
    
    if m>=2
    grad_beta=grad_beta-(V-miu(j))*(v{N_e}(m)-v{N_e}(m-1))/sigma(j);
    end
    
    grad_wmiu(j,:)=grad_wmiu(j,:)+(V-miu(j))*U{N_e}(m,:)/sigma(j);
    grad_wsigma(j,:)=grad_wsigma(j,:)+((V-miu(j))^2/2/sigma(j)-1/2)*U{N_e}(m,:);
end   %finish the gradient for beta,wmiu and wsigma of the current event

clearvars -except U lamda v x N_e grad w w_miu w_sigma beta grad_beta grad_wmiu grad_wsigma p g N_mode

end   %end for N_e

%here we explicitly penalize a beta that is close to 1 to avoid overfitting
p=p+1000*log(beta-1);
grad_beta=grad_beta+1000/(beta-1);

    for i=1:N_mode
        for j=1:N_mode
            g(1,N_mode*4*i+4*j-4*N_mode-3:4*N_mode*i+4*j-4*N_mode)=grad{i,j};
        end
            g(1,4*i+4*N_mode^2-3:4*i+4*N_mode^2)=grad_wmiu(i,1:4);
            g(1,4*i+4*N_mode^2+4*N_mode-3:4*i+4*N_mode^2+4*N_mode)=grad_wsigma(i,1:4);
    end
    g(1,4*(N_mode^2+2*N_mode)+1)=grad_beta;
    
    g=g.';   %make g into a column vector
    
    %add a penalty for the weight in front of U_longitude, to reduce
    %overfitting?
%     for k=1:9
%         p=p-50*x(4*k-3)^2;
%         g(4*k-3)=g(4*k-3)-100*x(4*k-3);
%     end
    
    
    p=-p;
    g=-g;    %here we want to maximize the log probability
    
    if isnan(p)
        bbbb
        'Warning:loglikelihood is NaN, line 116'
    end
    
end   %end function

