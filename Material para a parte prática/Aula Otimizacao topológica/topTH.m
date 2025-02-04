%%%% A 99 LINE TOPOLOGY OPTIMIZATION CODE BY OLE SIGMUND, JANUARY 2000 %%%
%%%% CODE MODIFIED FOR INCREASED SPEED, September 2002, BY OLE SIGMUND %%%
%%%%
%%%% ALTERADO PARA PROBLEMA DE CONDUÇÃO DE CALOR
%%%%
%function top(nelx,nely,volfrac,penal,rmin);
function top;
%
% DATA
nelx    = 60;
nely    = 60;
volfrac = 0.3;
penal   = 3.0;
rmin    = 1.25;
%
% INITIALIZE
x(1:nely,1:nelx) = volfrac;
% CHANGE 3
    %for ely = 1:nely
    %    for elx = 1:nelx
    %        if sqrt((ely-nely/2.)^2+(elx-nelx/3.)^2) <  nely/3.
    %        passive(ely,elx) = 1;
    %        x(ely,elx) = 0.001;
    %        else
    %        passive(ely,elx) = 0;
    %        end
    %    end
    %end
%
loop = 0; 
change = 1.;
% START ITERATION
while change > 0.01  
  loop = loop + 1;
  xold = x;
% FE-ANALYSIS
  [U]=FE(nelx,nely,x,penal);  
% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
  [KE] = lk;
  c = 0.;
  for ely = 1:nely
    for elx = 1:nelx
    %disp([elx,ely]);
      n1 = (nely+1)*(elx-1) + ely; 
      n2 = (nely+1)* elx    + ely;
      %Ue = U([2*n1-1;2*n1; 2*n2-1;2*n2; 2*n2+1;2*n2+2; 2*n1+1;2*n1+2],1);
      Ue = U([n1;n2;n2+1;n1+1],1);
      c = c + x(ely,elx)^penal*Ue'*KE*Ue;
      dc(ely,elx) = -penal*x(ely,elx)^(penal-1)*Ue'*KE*Ue;
      % CHANGE 2
      %  dc(ely,elx) = 0.;
      %  a=0.5;
      %  for i = 1:2
      %      if (i == 2);a= 1-a;end;
      %      Ue = a*U([2*n1-1;2*n1;2*n2-1;2*n2;2*n2+1;2*n2+2;2*n1+1;2*n1+2],i);
      %      c = c + x(ely,elx)^penal*Ue'*KE*Ue;
      %      dc(ely,elx) = dc(ely,elx) - penal*x(ely,elx)^(penal-1)*Ue'*KE*Ue;
      %  end
      %
    end
  end
% FILTERING OF SENSITIVITIES
  [dc]   = check(nelx,nely,rmin,x,dc);    
%disp(dc);pause;
% DESIGN UPDATE BY THE OPTIMALITY CRITERIA METHOD
    %x(find(passive)) = 0.001; % CHANGE 3 
  [x]    = OC(nelx,nely,x,volfrac,dc); 
% PRINT RESULTS
  change = max(max(abs(x-xold)));
  disp([' It.: ' sprintf('%4i',loop) ' Obj.: ' sprintf('%10.4f',c) ...
       ' Vol.: ' sprintf('%6.3f',sum(sum(x))/(nelx*nely)) ...
        ' ch.: ' sprintf('%6.3f',change )])
% PLOT DENSITIES  
  colormap(gray); imagesc(-x); axis equal; axis off;pause(1e-6);
  % CHANGE 4 - TopPlot
  %colormap(flipud(gray));
  %v=[0 volfrac];
  %A=contourf(x,v); set(gca,'YDir','reverse');axis equal; axis off;pause(1e-6)     
end 
%%%%%%%%%% OPTIMALITY CRITERIA UPDATE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xnew]=OC(nelx,nely,x,volfrac,dc)  
l1 = 0; l2 = 100000; move = 0.2;
while (l2-l1 > 1e-4)
  lmid = 0.5*(l2+l1);
  xnew = max(0.001,max(x-move,min(1.,min(x+move,x.*sqrt(dc./lmid)))));
  if sum(sum(xnew)) - volfrac*nelx*nely > 0; 
    l1 = lmid;
  else
    l2 = lmid;
  end
%  disp(lmid);pause;
%  disp(xnew);disp('x');pause;
end
%%%%%%%%%% MESH-INDEPENDENCY FILTER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dcn]=check(nelx,nely,rmin,x,dc)
dcn=zeros(nely,nelx);
for i = 1:nelx
  for j = 1:nely
    sum=0.0; 
    for k = max(i-floor(rmin),1):min(i+floor(rmin),nelx)
      for l = max(j-floor(rmin),1):min(j+floor(rmin),nely)
        fac = rmin-sqrt((i-k)^2+(j-l)^2);
        sum = sum+max(0,fac);
        dcn(j,i) = dcn(j,i) + max(0,fac)*x(l,k)*dc(l,k);
      end
    end
    dcn(j,i) = dcn(j,i)/(x(j,i)*sum);
  end
end
%%%%%%%%%% FE-ANALYSIS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [U]=FE(nelx,nely,x,penal)
[KE] = lk; 
K = sparse((nelx+1)*(nely+1),(nelx+1)*(nely+1));
F = sparse((nely+1)*(nelx+1),1); U = zeros((nely+1)*(nelx+1),1);
% K = sparse(2*(nelx+1)*(nely+1), 2*(nelx+1)*(nely+1));
% F = sparse(2*(nely+1)*(nelx+1),1); U = zeros(2*(nely+1)*(nelx+1),1);
    % CHANGE 2
    % F = sparse(2*(nely+1)*(nelx+1),2);U = sparse(2*(nely+1)*(nelx+1),2);
    %
for elx = 1:nelx
  for ely = 1:nely
    n1 = (nely+1)*(elx-1) + ely; 
    n2 = (nely+1)* elx    + ely;
    % edof = [2*n1-1; 2*n1; 2*n2-1; 2*n2; 2*n2+1; 2*n2+2; 2*n1+1; 2*n1+2];
    edof = [n1;n2;n2+1;n1+1];    
    K(edof,edof) = K(edof,edof) + x(ely,elx)^penal*KE;
  end
end
% DEFINE LOADS AND SUPPORTS (HALF MBB-BEAM)
%F((nelx+1)*(nely+1),1) = -1;
%fixeddofs   = [1,(nelx+1)];
F(1:(nelx+1))=-1;
F((nelx+1):(nelx+1):(nelx+1)*(nely+1))=-1;
fixeddofs   = [(nelx+1)*(nely+1)];
%fixeddofs   = [1,(nelx+1)];
%fixeddofs   = [1:1:(nelx+1)];
%fixeddofs = [(nelx+1)*(nely+1)]
%pause
%
    % CHANGE 1
    %F(2*(nelx+1)*(nely+1),1) = -1;
    %fixeddofs   = [1:2*(nely+1)];
%
    % CHANGE BEFORE 2 
    %F(2*(nelx+1)*(nely+1),1) = -1.;F(2*(nelx)*(nely+1)+2,1) = 1.;
%
    % CHANGE 2
    %F(2*(nelx+1)*(nely+1),1) = -1.;
    %F(2*(nelx)*(nely+1)+2,2) = 1.;
%
%alldofs     = [1:2*(nely+1)*(nelx+1)];
alldofs     = [1:(nely+1)*(nelx+1)];
freedofs    = setdiff(alldofs,fixeddofs);
% SOLVING
U(freedofs,:) = K(freedofs,freedofs) \ F(freedofs,:);      
U(fixeddofs,:)= 0;
%%%%%%%%%% ELEMENT STIFFNESS MATRIX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [KE]=lk
E = 1.; 
nu = 0.3;
%k=[ 1/2-nu/6   1/8+nu/8 -1/4-nu/12 -1/8+3*nu/8 ... 
%   -1/4+nu/12 -1/8-nu/8  nu/6       1/8-3*nu/8];
%KE = E/(1-nu^2)*[ k(1) k(2) k(3) k(4) k(5) k(6) k(7) k(8)
%                  k(2) k(1) k(8) k(7) k(6) k(5) k(4) k(3)
%                  k(3) k(8) k(1) k(6) k(7) k(4) k(5) k(2)
%                  k(4) k(7) k(6) k(1) k(8) k(3) k(2) k(5)
%                  k(5) k(6) k(7) k(8) k(1) k(2) k(3) k(4)
%                  k(6) k(5) k(4) k(3) k(2) k(1) k(8) k(7)
%                  k(7) k(4) k(5) k(2) k(3) k(8) k(1) k(6)
%                  k(8) k(3) k(2) k(5) k(4) k(7) k(6) k(1)];
%
KE = [-0.6666666666667  0.1666666666667  0.3333333333333  0.1666666666667     
       0.1666666666667 -0.6666666666667  0.1666666666667  0.3333333333333     
       0.3333333333333  0.1666666666667 -0.6666666666667  0.1666666666667     
       0.1666666666667  0.3333333333333  0.1666666666667 -0.6666666666667];
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This Matlab code was written by Ole Sigmund, Department of Solid         %
% Mechanics, Technical University of Denmark, DK-2800 Lyngby, Denmark.     %
% Please sent your comments to the author: sigmund@fam.dtu.dk              %
%                                                                          %
% The code is intended for educational purposes and theoretical details    %
% are discussed in the paper                                               %
% "A 99 line topology optimization code written in Matlab"                 %
% by Ole Sigmund (2001), Structural and Multidisciplinary Optimization,    %
% Vol 21, pp. 120--127.                                                    %
%                                                                          %
% The code as well as a postscript version of the paper can be             %
% downloaded from the web-site: http://www.topopt.dtu.dk                   %
%                                                                          %
% Disclaimer:                                                              %
% The author reserves all rights but does not guaranty that the code is    %
% free from errors. Furthermore, he shall not be liable in any event       %
% caused by the use of the program.                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
