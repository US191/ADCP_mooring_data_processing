function y = medfilt1_perso(x,w)
%MEDFILT1_PERSO 
%Gere la non presence de signal processing toolbox
%
%odd w
if(mod(w,1)==1)
    w=w+1;
end
%
n=length(x);
for i=1:n
    if(i<=w/2)
       y(i)=median(x(1:i+w/2)); 
    elseif(i>=n-w/2)
       y(i)=median(x(i-w/2:end)); 
    else
       y(i)=median(x(i-w/2:i+w/2)); 
    end
end

end

