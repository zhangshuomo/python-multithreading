function [tau, omega, omega_d, v, a, a_c] = rnea(leg,q,dq,ddq)

tau=zeros(3,1);

if leg<=2
    front=1;
else
    front=-1;
end
if mod(leg,2)
    left=1;
else
    left=-1;
end

% compute joints & link-CoMs transform first
% saved in joint_trans and com_trans
joint_trans=zeros(4,4,4);
com_trans=zeros(4,4,4);

x1=0.2407*front;
y1=0.051*left;
y2=0.0868*left;
z3=-0.25;
z4=-0.25;

% hip joint transform
hip_trans=[1, 0, 0, x1;...
         0,cos(q(1)),-sin(q(1)),y1;...
         0,sin(q(1)),cos(q(1)),0;...
         0,0,0,1];
joint_trans(:,:,1)=hip_trans;

% thigh joint transform
thigh_trans=[cos(q(2)) 0 sin(q(2)) 0;...
         0 1 0 y2;...
         -sin(q(2)) 0 cos(q(2)) 0;...
         0 0 0 1];
joint_trans(:,:,2)=joint_trans(:,:,1)*thigh_trans;

% calf joint transform
calf_trans=[cos(q(3)) 0 sin(q(3)) 0;...
         0 1 0 0;...
         -sin(q(3)) 0 cos(q(3)) z3;...
         0 0 0 1];
joint_trans(:,:,3)=joint_trans(:,:,2)*calf_trans;

% foot position fixed joint
foot_trans=[1 0 0 0;...
            0 1 0 0;...
            0 0 1 z4;...
            0 0 0 1];
joint_trans(:,:,4)=joint_trans(:,:,3)*foot_trans; 

% setting com_trans as the transform from 
% corresponding joint origin to CoM, first
com_trans(:,:,1)=[1 0 0 -0.022191*front;...
                  0 1 0 0.015144*left;...
                  0 0 1 -1.5e-05;...
                  0 0 0 1];
com_trans(:,:,2)=[1 0 0 -0.005607;...
                  0 1 0 -0.003877*left;...
                  0 0 1 -0.048199;...
                  0 0 0 1];
com_trans(:,:,3)=[1 0 0 0.002781;...
                  0 1 0 6.3e-05;...
                  0 0 1 -0.142518;...
                  0 0 0 1];
com_trans(:,:,4)=[1 0 0 0;...
                  0 1 0 0;...
                  0 0 1 0;...
                  0 0 0 1];
for link=1:4
    com_trans(:,:,link)=joint_trans(:,:,link)*com_trans(:,:,link);
end

omega=zeros(3,4);
omega_d=zeros(3,4);
v=zeros(3,4);
a=zeros(3,4);
a_c=zeros(3,4);
% forward iteration
for link=1:4
    if link==1
        axis=1;
        omega(:,link)=dq(link)*joint_trans(1:3,axis,link);
        omega_d(:,link)=ddq(link)*joint_trans(1:3,axis,link);
        p_c=com_trans(1:3,4,link)-joint_trans(1:3,4,link);
        v(:,link)=zeros(3,1);
%         a(:,link)=[0,0,9.81];
        a_c(:,link)=a(:,link)+cross(omega_d(:,link),p_c)+cross(omega(:,link),cross(omega(:,link),p_c));
    else 
      if link==2||link==3
        axis=2;
        omega(:,link)=omega(:,link-1)+dq(link)*joint_trans(1:3,axis,link);
        omega_d(:,link)=omega_d(:,link-1)+ddq(link)*joint_trans(1:3,axis,link)+dq(link)*cross(omega(:,link-1),joint_trans(1:3,axis,link));
        p=joint_trans(1:3,4,link)-joint_trans(1:3,4,link-1);
        p_c=com_trans(1:3,4,link)-joint_trans(1:3,4,link);
        v(:,link)=v(:,link-1)+cross(omega(:,link-1),p);
        a(:,link)=a(:,link-1)+cross(omega_d(:,link-1),p)+cross(omega(:,link-1),cross(omega(:,link-1),p));
        a_c(:,link)=a(:,link)+cross(omega_d(:,link),p_c)+cross(omega(:,link),cross(omega(:,link),p_c));
      else % link==4
        omega(:,link)=omega(:,link-1);
        omega_d(:,link)=omega_d(:,link-1);
        p=joint_trans(1:3,4,link)-joint_trans(1:3,4,link-1);
        p_c=com_trans(1:3,4,link)-joint_trans(1:3,4,link);
        v(:,link)=v(:,link-1)+cross(omega(:,link-1),p);
        a(:,link)=a(:,link-1)+cross(omega_d(:,link-1),p)+cross(omega(:,link-1),cross(omega(:,link-1),p));
        a_c(:,link)=a(:,link)+cross(omega_d(:,link),p_c)+cross(omega(:,link),cross(omega(:,link),p_c));
      end
    end
end

mass=[1.993, 0.639, 0.207, 0.06];
inertials=zeros(3,3,4);
inertials(:,:,1)=[0.002903894 -0.000071850*left*front -0.000001262*front;
                  -0.000071850*left*front 0.004907517 -0.00000175*left;
                  -0.000001262*front -0.00000175*left 0.005586944];
inertials(:,:,2)=[0.005666803 0.000003597*left 0.000491446;
                  0.000003597*left 0.005847229 0.000010086*left;
                  0.000491446 0.000010086*left 0.000369811];
inertials(:,:,3)=[0.006341369 -0.000000003 -0.000087951;
                  -0.000000003 0.006355157 -0.000001336;
                  -0.000087951 -0.000001336 0.000039188];
inertials(:,:,4)=[1.6854e-05 0.0 0.0;
                  0.0 1.6854e-05 0.0;
                  0.0 0.0 1.6854e-05];
forces=zeros(3,5);
moment=zeros(3,5);

% backward iteration
for link=4:-1:1
    forces(:,link)=forces(:,link+1)+mass(link)*a_c(:,link)-mass(link)*[0;0;-9.81];
    if link==4
        p_f=zeros(3,1);
    else
        p_f=joint_trans(1:3,4,link+1)-com_trans(1:3,4,link);
    end
    p_b=joint_trans(1:3,4,link)-com_trans(1:3,4,link);
    inertial_w=com_trans(1:3,1:3,link)*inertials(:,:,link)*com_trans(1:3,1:3,link)';
    moment(:,link)=moment(:,link+1)-cross(p_b,forces(:,link))+cross(p_f,forces(:,link+1))+...
        inertial_w*omega_d(:,link)+cross(omega(:,link),inertial_w*omega(:,link));
end
for joint=1:3
    if joint==1
        axis=1;
    else
        axis=2;
    end
    tau(joint)=moment(:,joint)'*joint_trans(1:3,axis,joint);
end
end

