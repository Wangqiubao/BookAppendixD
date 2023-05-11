function [f,g] = sparseGalerkin(yin,ahat,polyorder,usesine,delay,usee)
nVars =length(yin);
yout = [];
for i=1:nVars
    yout = [yout,yin(:,i)];
end

if(polyorder>=2)
    % poly order 2
    for i=1:nVars
        for j=i:nVars
            yout =  [yout,yin(:,i).*yin(:,j)];
        end
    end
end

if(polyorder>=3)
    % poly order 3
    for i=1:nVars
        for j=i:nVars
            for k=j:nVars
                yout =  [yout,yin(:,i).*yin(:,j).*yin(:,k)];
            end
        end
    end
end

if(polyorder>=4)
    % poly order 4
    for i=1:nVars
        for j=i:nVars
            for k=j:nVars
                for l=k:nVars
                    yout= [yout,yin(:,i).*yin(:,j).*yin(:,k).*yin(:,l)];
                    ind = ind+1;
                end
            end
        end
    end
end

if(polyorder>=5)
    % poly order 5
    for i=1:nVars
        for j=i:nVars
            for k=j:nVars
                for l=k:nVars
                    for m=l:nVars
                        yout = [yout,yin(:,i).*yin(:,j).*yin(:,k).*yin(:,l).*yin(:,m)];
                        ind = ind+1;
                    end
                end
            end
        end
    end
end
if(usesine)
    for k=1:10;
        yout = [yout sin(k*yin) cos(k*yin)];
    end
end
%%
syms x1tau GSN1 GSN2
if(usee)
    f = ([1,yout,exp(yin(1)),exp(yin(2))]*ahat(:,1));
    g = ([1,yout,exp(yin(1)),exp(yin(2))]*ahat(:,2));
else
    f =  ([1,yout]*ahat(:,1));
    g =  ([1,yout]*ahat(:,2));
end
f = eval(['@(x1tau,GSN1,GSN2,x1,x2)',vectorize(f)]);
g = eval(['@(x1tau,GSN1,GSN2,x1,x2)',vectorize(g)]);