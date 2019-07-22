% Author :  Qida Yu                                                        
% This programe is implemented in matlab 2018a      
% Address:  Nanjing University of Aeronautics and Astronautics              


function [ R_cw,t_cw ] = OPnL_main( p1,p2,W1,W2 )

nLine = length(p1);
p1 = [p1; ones(1,nLine)];
p2 = [p2; ones(1,nLine)];

%计算平面法向量;
n_c=xnorm(cross(p1,p2));

%用来计算旋转矩阵和平移矩阵关系的系数;
Q=zeros(2*nLine,10);
N=zeros(2*nLine,3);
nx=n_c(1,:)'; ny=n_c(2,:)'; nz=n_c(3,:)';
Px=W1(1,:)'; Py=W1(2,:)'; Pz=W1(3,:)';
Q(1:2:end,:)=[Px.*nx+Py.*ny+Pz.*nz,...
          2*Py.*nz-2*Pz.*ny,...
          2*Pz.*nx-2*Px.*nz,...
          2*Px.*ny-2*Py.*nx,...
          Px.*nx-Py.*ny-Pz.*nz,...
          2.*Px.*ny+2*Py.*nx,...
          2.*Px.*nz+2*Pz.*nx,...
          Py.*ny-Px.*nx-Pz.*nz,...
          2*Py.*nz+2*Pz.*ny,...
          Pz.*nz-Px.*nx-Py.*ny];
Px=W2(1,:)'; Py=W2(2,:)'; Pz=W2(3,:)';
Q(2:2:end,:)=[Px.*nx+Py.*ny+Pz.*nz,...
          2*Py.*nz-2*Pz.*ny,...
          2*Pz.*nx-2*Px.*nz,...
          2*Px.*ny-2*Py.*nx,...
          Px.*nx-Py.*ny-Pz.*nz,...
          2.*Px.*ny+2*Py.*nx,...
          2.*Px.*nz+2*Pz.*nx,...
          Py.*ny-Px.*nx-Pz.*nz,...
          2*Py.*nz+2*Pz.*ny,...
          Pz.*nz-Px.*nx-Py.*ny];
N(1:2:end,:)=-n_c.';
N(2:2:end,:)=-n_c.';
%计算最后的矩阵;
C=((N.'*N)\N.')*Q;
E=Q-N*C;
G=E.'*E;    %the size of G is 10*10;

g12 = G(1,2); g13 = G(1,3); g14 = G(1,4); g15 = G(1,5); g16 = G(1,6); g17 = G(1,7); g18 = G(1,8); g19 = G(1,9); g110 = G(1,10);
g22 = G(2,2); g23 = G(2,3); g24 = G(2,4); g25 = G(2,5); g26 = G(2,6); g27 = G(2,7); g28 = G(2,8); g29 = G(2,9); g210 = G(2,10);
g33 = G(3,3); g34 = G(3,4); g35 = G(3,5); g36 = G(3,6); g37 = G(3,7); g38 = G(3,8); g39 = G(3,9); g310 = G(3,10);
g44 = G(4,4); g45 = G(4,5); g46 = G(4,6); g47 = G(4,7); g48 = G(4,8); g49 = G(4,9); g410 = G(4,10);
g55 = G(5,5); g56 = G(5,6); g57 = G(5,7); g58 = G(5,8); g59 = G(5,9); g510 = G(5,10);
g66 = G(6,6); g67 = G(6,7); g68 = G(6,8); g69 = G(6,9); g610 = G(6,10);
g77 = G(7,7); g78 = G(7,8); g79 = G(7,9); g710 = G(7,10);
g88 = G(8,8); g89 = G(8,9); g810 = G(8,10);
g99 = G(9,9); g910 = G(9,10);
g1010 = G(10,10);

M(1,:)=[2*g55,3*g56,3*g57,3*g25,g68,2*g58+g66,g69+g78,g28+g36,g710,g77+2*g510,g79+g610,g47+g210,...
        2*g59+2*g67,2*g26+2*g35,2*g27+2*g45,g29+g37+g46,2*g15+g22,g16+g23,g17+g24,g12];
M(2,:)=[g56,2*g58+g66,g59+g67,g26+g35,2*g88,3*g68,3*g89,3*g38,g910,g79+g610,g99+2*g810,g310+g49,2*g69+2*g78,...
        2*g28+2*g36,g29+g37+g46,2*g39+2*g48,g16+g23,2*g18+g33,g19+g34,g13];
M(3,:)=[g57,g59+g67,g77+2*g510,g27+g45,g89,g69+g78,2*g810+g99,g39+g48,2*g1010,3*g710,3*g910,3*g410,...
        2*g79+2*g610,g29+g37+g46,2*g47+2*g210,2*g49+2*g310,g17+g24,g19+g34,g44+2*g110,g14];
    
solution = OPnL_solver(M);

for i=1:size(solution,2)
    s1=solution(1,i);
    s2=solution(2,i);
    s3=solution(3,i);
    Rs=[1+s1^2-s2^2-s3^2,2*s1*s2-2*s3,2*s1*s3+2*s2;
       2*s1*s2+2*s3,1-s1^2+s2^2-s3^2,2*s2*s3-2*s1;
       2*s1*s3-2*s2,2*s2*s3+2*s1,1-s1^2-s2^2+s3^2];
   factor=1/(1+s1^2+s2^2+s3^2);
   s=[1,s1,s2,s3,s1^2,s1*s2,s1*s3,s2^2,s2*s3,s3^2].';
   R_cw(:,:,i)=Rs*factor;
   t_cw(:,i)=C*s*factor;
end

end

