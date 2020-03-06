function [xyzt,out] = triangulate360(xyt1,xyt2,t,R,v)
%TRIANGULATE360 Summary of this function goes here
%   Detailed explanation goes here

% initialize
tmin = min(xyt1(:,3));
tmax = max(xyt1(:,3));

m = cell(1,tmax);
xyzt = [];
points1 = [];
points2 = [];
alpha1t = [];
alpha2t = [];
beta2t = [];
xyz1 = cell(1,tmax);
xyz2 = cell(1,tmax);
xyz = cell(1,tmax);

for i = tmin:tmax
    
    xy1 = xyt1(xyt1(:,3) == i,:);
    xy2 = xyt2(xyt2(:,3) == i,:);
        
    if i > tmin
        warning off
    end
    
    if nargin == 5
        alpha1 = xy2alpha(xy1,v);
        alpha2 = xy2alpha(xy2,v);
    else
        alpha1 = xy2alpha(xy1);
        alpha2 = xy2alpha(xy2);
    end
            
    warning on
    
    m{i} = matchPointsTR(alpha1,alpha2,t,R);    
    
    for j = 1:length(m{i})
        
        if m{i}(j)>0
            a1 = alpha1(j,:);
            a2 = alpha2(m{i}(j),:);
            b2 = (R'*a2')';
            C = [a1' -b2'];
            
            rr = lsqnonneg(C,t);
            r1 = rr(1);
            r2 = rr(2);
            
            P1 = r1*alpha1(j,:);
            P2 = r2.*(R'*alpha2(m{i}(j),:)')' + t(:)';
            
            P = (P1+P2)/2;
            
            xyzt = vertcat(xyzt,[P i]);
            
            points1 = vertcat(points1,xy1(j,:));
            points2 = vertcat(points2,xy2(m{i}(j),:));
            alpha1t = vertcat(alpha1t,[a1 i]);
            alpha2t = vertcat(alpha2t,[a2 i]);
            beta2t = vertcat(beta2t,[b2 i]);
            
            xyz1{i} = vertcat(xyz1{i},P1);
            xyz2{i} = vertcat(xyz2{i},P2);
            xyz{i} = vertcat(xyz{i},P);
                
        end
    
    end

    w = waitbar((i-tmin)/(tmax-tmin));
    
end
close(w)

out.xyz1 = xyz1;
out.xyz2 = xyz2;
out.xyz = xyz;
out.m = m;
out.xyt1 = points1;
out.xyt2 = points2;
out.alpha1t = alpha1t;
out.alpha2t = alpha2t;
out.beta2t = beta2t;

warning on

end

