clc
clear all
close all



b = 8;
h = 6;
n_b = 4;
n_h = 3;
G = 200;
theta = 10;

x_coord = repmat(linspace(0,b,n_b+1),1,n_h+1);
y_coord = repelem(linspace(0,h,n_h+1), n_b+1);
n_info = [x_coord.' , y_coord.']

e_matrix=[0 0 0 0];  %initialize e_info, will delete this row at the end
for i=1:n_h  %Loop through cells vertically starting on top
  for j=1:n_b  %Loop through horizotally starting on left each time
    top_lef=(i-1)*(n_b+1)+j;
    top_rig=top_lef+1;
    bot_lef=top_lef+(n_b+1);
    bot_rig=top_rig+(n_b+1);
    e_matrix=[e_matrix;bot_lef bot_rig top_rig top_lef];
  end
end
e_matrix(1,:)=[];
e_info = fliplr(e_matrix)


all_node = [1:(n_b+1)*(n_h+1)]'
all_matrix = flipud(reshape(all_node,n_b+1,n_h+1).')
all_matrix(1,:) = 0;
all_matrix(end,:) = 0;
all_matrix(:,1) = 0;
all_matrix(:,end) = 0
id_f = sort(nonzeros(all_matrix))
id_s = setdiff(all_node, id_f)

num_node = (n_b+1)*(n_h+1);
K = zeros(num_node, num_node);
num_e=size(e_info,1);
for e=1:num_e
  i=e_info(e,1);
  j=e_info(e,2);
  k=e_info(e,3);
  l=e_info(e,4);
  Lx=n_info(j,1)-n_info(i,1);  %Horizontal element size, is constant for this problem but this would work for any
  Ly=n_info(k,2)-n_info(j,2);   %Vertical element size
  r=Ly/Lx;
  k11=(1+r^2)/(3*r);
  k12=1/(6*r)-r/3;
  k13=-(1+r^2)/(6*r);
  k14=(-2+r^2)/(6*r);
  ke=(1/G)*[k11 k12 k13 k14;    
        k12 k11 k14 k13;
        k13 k14 k11 k12;
        k14 k13 k12 k11];

  K(i,i)=K(i,i)+ke(1,1);
  K(i,j)=K(i,j)+ke(1,2);
  K(i,k)=K(i,k)+ke(1,3);
  K(i,l)=K(i,l)+ke(1,4);

  K(j,i)=K(j,i)+ke(2,1);
  K(j,j)=K(j,j)+ke(2,2);
  K(j,k)=K(j,k)+ke(2,3);
  K(j,l)=K(j,l)+ke(2,4);

  K(k,i)=K(k,i)+ke(3,1);
  K(k,j)=K(k,j)+ke(3,2);
  K(k,k)=K(k,k)+ke(3,3);
  K(k,l)=K(k,l)+ke(3,4);

  K(l,i)=K(l,i)+ke(4,1);
  K(l,j)=K(l,j)+ke(4,2);
  K(l,k)=K(l,k)+ke(4,3);
  K(l,l)=K(l,l)+ke(4,4);
end
K

Pfef = zeros(num_node, 1);
num_elm = size(e_info,1);
for e = 1:num_elm
    i_node = e_info(e,1);
    j_node = e_info(e,2);
    k_node = e_info(e,3);    
    l_node = e_info(e,4);
    Lx=n_info(j_node,1)-n_info(i_node,1);
    Ly=n_info(k_node,2)-n_info(j_node,2);  
    Pfefe = Lx*Ly*2 *theta*[1/4;1/4;1/4;1/4];
    Pfef(i_node) = Pfef(i_node)+Pfefe(1);
    Pfef(j_node) = Pfef(j_node)+Pfefe(2);
    Pfef(l_node) = Pfef(l_node)+Pfefe(3);
    Pfef(k_node) = Pfef(k_node)+Pfefe(4);
end
Pfef

P = zeros(num_node,1);
Pf= P(id_f);
Pfeff = Pfef(id_f);
Kff=K(id_f,id_f);
df = inv(Kff) * (Pfeff);
Pfeff
Kff
df
d = zeros(num_node,1);
d(id_f) = df;

phi = d;
area_elm = (b/n_b)*(h/n_h)
vol = sum(area_elm.*phi);
J = 2*vol/(G*theta)













