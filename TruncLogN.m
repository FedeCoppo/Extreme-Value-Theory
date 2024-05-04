% Truncated LogNormal Density

function f=TruncLogN(x,mu,sigma,t)
f=zeros(length(x),1);
indici2=find(x>=t(1) & x<=t(2));
x2=x(indici2);
den=logncdf(t(2),mu,sigma)-logncdf(t(1),mu,sigma);
f(indici2)=lognpdf(x2,mu,sigma)./den; % truncated logn density
end
