function [alpha] = xy2alpha(xy,v)
%XY2ALPHA Converts the xy coordinates in the equirectangular frame to
%   spherical coordinates alpha
%
%   xy is a list of positions (Nx2, subsequent colums are not considered)
%   v informs about size of equirectangular frames; it is either: 
%       - a movie structure 
%       - an array containing [width height]
%
%   alpha is Nx3 array
%   
%
% Raphael Sarfati, 03/2020
% Peleg Lab, University of Colorado Boulder

[theta,phi] = xy2thetaphi(xy,v);

alpha = horzcat(cos(theta).*sin(phi),sin(theta).*sin(phi),cos(phi));

end

