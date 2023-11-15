clc
close all
clear all

nfig=0;

a=1;   % side in mm



unit_nodes_x=[0,a,a,0];
unit_nodes_y=[0,0,a,a];
unit_nodes_z=[a,a,a,a];


cell_nodes=[
    0,0,0;
    a,0,0;
    a,a,0;
    0,a,0;
    a/2,a/2,0;
    0,a/2,a/2;
    a/2,0,a/2;
    a/2,a,a/2;
    a,a/2,a/2;
    0,0,a;
    a,0,a;
    a,a,a;
    0,a,a;
    a/2,a/2,a;
];

beam_cell=[
    1,5;
    1,6;
    1,7;
    2,5;
    2,7;
    2,9;
    3,5;
    3,8;
    3,9;
    4,5;
    4,6;
    4,8;
    5,6;
    5,7;
    5,8;
    5,9;
    6,7;
    6,8;
    7,9;
    8,9;
    6,14;
    7,14;
    8,14;
    9,14;
    10,14;
    10,6;
    10,7;
    11,14;
    11,7;
    11,9;
    12,14;
    12,8;
    12,9;
    13,14;
    13,6;
    13,8;
    

    

%     3,4
    ];








size=30;

box_x=[0,5];
%3*sqrt(5)];

box_y=[0,3*sqrt(5)];

box_z=[0,3*sqrt(5)];



%omega=atand(1/2);
omega=0;	% angle de rotation du crystal autour de l'axe X (à 0 pour l'instant)

M=[1,0,0;
    0,cosd(omega),-sind(omega);
    0,sind(omega),cosd(omega);
];


flag=0;
num_node=0;

% start_i=1;
% start_j=5;
% start_k=1;
% 
% in_x=9;
% in_y=10;
% in_z=7;

start_i=1;
start_j=1;
start_k=5;

in_x=5;
in_y=10;
in_z=10;


for i=1:in_x
    for j=1:in_y
        for k=1:in_z
            
            if flag==0
                nodes=cell_nodes;
                
                nodes(:,1)=nodes(:,1)+(i-start_i)*a;
                nodes(:,2)=nodes(:,2)+(j-start_j)*a;
                nodes(:,3)=nodes(:,3)+(k-start_k)*a;
                
                beams=beam_cell;
                flag=1; 
                
            else
 
            
            
            
            new_nodes=cell_nodes;
            new_nodes(:,1)=new_nodes(:,1)+(i-start_i)*a;
            new_nodes(:,2)=new_nodes(:,2)+(j-start_j)*a;
            new_nodes(:,3)=new_nodes(:,3)+(k-start_k)*a;
            
            

            
            
            
            new_beams_old_id=beam_cell;
            new_beams_new_id=beam_cell+length(nodes);
            
            [Lia, Locb] = ismember(nodes,new_nodes, 'rows');
            [Lic, Loc_new_in_old] = ismember(new_nodes,nodes, 'rows');
            
            Loc_new_in_old(Loc_new_in_old==0)=[];
            for z=1:length(Loc_new_in_old)
                [row,col]=find(new_beams_old_id(:,1)==Loc_new_in_old(z));
                
%                 for zz=1:length(row)
%                     new_beams_new_id(row(zz),col(zz))=Loc_new_in_old(z);
%                 end
                   
                
                
                
            end
%             return
%             new_beams_new_id=
            
%             Locb(Locb==0)=[];            
%             new_nodes(Locb(Locb~=0),:)=[];
%             new_nodes(Locb,:)=0;
            
            nodes=[nodes;new_nodes];
            beams=[beams;new_beams_new_id];
            end
        end
    end
end

   
nodes=M*nodes';
nodes=nodes';

%% Cut
epsilon=0.01;
for i=1:length(nodes)
    if ((nodes(i,1)>=(box_x(1))) && (nodes(i,1)<=(box_x(2))) && (nodes(i,2)>=(box_y(1))) && (nodes(i,2)<=(box_y(2))) && (nodes(i,3)>=(box_z(1))) && (nodes(i,3)<=(box_z(2))))

%         disp('in')
        flag_in(i,1)=1;
    else
%         disp('out')
        flag_in(i,1)=0;
    end
%     if (i==4)
%         return
%     end
    
end


nfig=nfig+1;
fig=figure(nfig);
hold on
% scatter3(nodes(:,1),nodes(:,2),nodes(:,3),80,'b','filled');


id=1:length(nodes);
id(flag_in==0)=[];

    
[La, Loc]=ismember(beams(:,1),id);
[Lb, Loc]=ismember(beams(:,2),id);
L=La+Lb;

beams(L~=2,:)=[];
   

id_nodes=transpose(1:length(nodes));
k=0;
for i=1:length(id_nodes)
    if flag_in(i)==1
        k=k+1;
        new_id_nodes(i,1)=k;
    else
        new_id_nodes(i,1)=-1;
    end
    
    
end
nodes(flag_in==0,:)=[];
id_nodes(flag_in==0,:)=[];
% 
% [Lic, Loc_new_in_old] = ismember(new_nodes,nodes, 'rows');

% return

for i=1:length(beams)
    beams(i,1)=new_id_nodes(beams(i,1));
    beams(i,2)=new_id_nodes(beams(i,2));
end
% 
% new_id_nodes(flag_in==0,:)=[];
% [La, Loc]=ismember(new_id_nodes,beams(:,1));
% [Lb, Loc]=ismember(new_id_nodes,beams(:,2));
% 
% L=La+Lb;
% 
% beams(L~=2,:)=[];


% [Lb, Loc]=ismember(beams(:,2),id);
% 
% for i=1:length(new_id_nodes)
%     
%     [La, Loc]=ismember(beams(:,1),id);
%     
% end

%% Plot






for i=1:length(beams)
    node_1=beams(i,1);
    node_2=beams(i,2);
    
%     p=
%     id_node_1=find(id_nodes==beams(i,1));
%     id_node_2=find(id_nodes==beams(i,2));
    
    if node_2==107
        disp('here')
    end
    
    node_pos=[nodes(node_1,1),nodes(node_1,2),nodes(node_1,3);
        nodes(node_2,1),nodes(node_2,2),nodes(node_2,3)
        ];
        
    
%     plot3([cell(node_1,1),cell(node_2,1)],[cell(node_1,2),cell(node_2,2)],[cell(node_1,3),cell(node_2,3)],'k','LineWidth',2);
    s=plot3(node_pos(:,1),node_pos(:,2),node_pos(:,3),'LineWidth',4,'Color',[0,0,0,0.5]);


end


lnodes=length(nodes);
lbeams=length(beams);
lbeams0=lbeams;
epsilon=0.05;
for i=1:length(nodes)
    if ((nodes(i,3)<=epsilon) && (nodes(i,3)>=-epsilon))
        lnodes=lnodes+1;
        nodes(lnodes,:)=nodes(i,:);
        nodes(lnodes,3)=-0.2;
       
        lbeams=lbeams+1;
        beams(lbeams,:)=[lnodes,i];
    end
         
end


scatter3(nodes(:,1),nodes(:,2),nodes(:,3),80,'r','filled');

% scatter3(nodes(107,1),nodes(107,2),nodes(107,3),80,'g','filled');
% scatter3(nodes(108,1),nodes(108,2),nodes(108,3),80,'m','filled');

% scatter3(nodes(flag_in==1,1),nodes(flag_in==1,2),nodes(flag_in==1,3),80,'r','filled');


plot3([0,0],[0,0],[box_z(1),box_z(2)],'-r','LineWidth',4);
% scatter3(new_nodes(:,1),new_nodes(:,2),new_nodes(:,3),80,'r','filled');

pbaspect([1 1 1])
xlim([-1 8]);
ylim([-1 8]);
zlim([-1 8]);
xlabel('x');
ylabel('y');
zlabel('z');

view([90,0]);

r_node=0.05;
r_beam=0.05;

r_beams=r_beam*ones(length(beams),1);
r_beams(lbeams0+1:end)=0.01;
beams=beams;
create_file(nodes,beams,r_node,r_beams);



%%
% [pos_and_angle] = beam_pos_and_angle(nodes,beams);

