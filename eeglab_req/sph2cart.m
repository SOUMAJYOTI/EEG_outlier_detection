% Copyright (C) 2000  Kai Habel
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
%{
## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

## -*- texinfo -*-
## @deftypefn {Function File} {} [@var{X,Y,Z}] = sph2cart (@var{Theta},@var{Phi},@var{R})
## transforms spherical to cartesian coordinates.
## @var{X},@var{Y} and @var{Z} must be of same shape.
## @var{Theta} describes the angle relative to the x - axis.
## @var{Phi} is the angle relative to the xy - plane.
## @var{R} is the distance to the origin (0,0,0).
## @end deftypefn
## @seealso{pol2cart,cart2pol,cart2sph}

## Author:	Kai Habel <kai.habel@gmx.de>
%}

function [X, Y, Z] = sph2cart (Theta, Phi, R)

  for i=1:length(Theta)
      X(i) = R(i)*cos (Phi(i))*cos (Theta(i));
      Y(i) = R(i)*cos (Phi(i))*sin (Theta(i));
      Z(i) = R(i)*sin (Phi(i));
  end

end