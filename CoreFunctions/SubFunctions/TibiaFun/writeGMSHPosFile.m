function Mat3D_final = writeGMSHPosFile(Bone,Implt,CortSize,IntfcSize, maxSize,Nxp,Oxp)
%UNTITLED Write a file that will contain a discrete field that control the
% element size of the mesh of the implant
%   Goal is to refine mesh on cortical surface and aroun the implant
Pts_bone = Bone.Points;
Pts_implt = Implt.Points;


% Remove points too close to intfc to get a uniform filed at the interface
Pts_bone_on_intfc = abs((Pts_bone*Nxp-Oxp*Nxp))<4;
Pts_bone = Pts_bone(~Pts_bone_on_intfc,:);

% Keep only implnat points on the interface with the bone
PtsOnIntfc = (Pts_implt*Nxp-Oxp*Nxp)<0.1;
Pts_implt = Pts_implt(PtsOnIntfc,:);

% Rmove a layer of points under the plate
PtsOnIntfc = (Pts_implt*Nxp-Oxp*Nxp)>-0.2 | (Pts_implt*Nxp-Oxp*Nxp)<-4.2;
Pts_implt = Pts_implt(PtsOnIntfc,:);


X_implt=Pts_implt(:,1); Y_implt=Pts_implt(:,2); Z_implt=Pts_implt(:,3); 
X_bone=Pts_bone(:,1); Y_bone=Pts_bone(:,2); Z_bone=Pts_bone(:,3);

Xmin=min([min(X_implt) min(X_bone)]);
Ymin=min([min(Y_implt) min(Y_bone)]);
Zmin=min([min(Z_implt) min(Z_bone)]);

Xmax=max([max(X_implt) max(X_bone)]);
Ymax=max([max(Y_implt) max(Y_bone)]);
Zmax=max([max(Z_implt) max(Z_bone)]);

X_t_implt = X_implt-Xmin;
Y_t_implt = Y_implt-Ymin;
Z_t_implt = Z_implt-Zmin;

X_t_bone = X_bone-Xmin;
Y_t_bone = Y_bone-Ymin;
Z_t_bone = Z_bone-Zmin;

% Initial discrete field
Mat3D=ones(ceil(Xmax-Xmin)+20,ceil(Ymax-Ymin)+20,ceil(Zmax-Zmin)+20);
Mat3D = (1./maxSize)*Mat3D;

Mat3D_implt = Mat3D;
invIntfcSize = 4/IntfcSize;
for l=1:length(Pts_implt)
    Mat3D_implt(round(X_t_implt(l))+10,round(Y_t_implt(l))+10,round(Z_t_implt(l))+10)=invIntfcSize;
end

Mat3D_bone = Mat3D;
invCortSize = 4/CortSize;
for l=1:length(Pts_bone)
    Mat3D_bone(round(X_t_bone(l))+10,round(Y_t_bone(l))+10,round(Z_t_bone(l))+10)=...
        invCortSize/max(1,(4.6-4.0*Z_t_bone(l)/ceil(Zmax-Zmin)));
end

Mat3D_blurred_implt = imgaussfilt3(Mat3D_implt,5);
Mat3D_blurred_bone = imgaussfilt3(Mat3D_bone,5.5);

Mat3D_blurred_max = max(Mat3D_blurred_implt,Mat3D_blurred_bone);

Mat3D_final=Mat3D_blurred_max.^-1;

%% Ecriture d'un fichier .pos (post processing field de gmsh) pour utiliser
% en tant que background mesh, pour controler la taille des elements en
% periprothetique

fID2=fopen('outGMSHField.pos','w');
fprintf(fID2,'View "background mesh" {\r\n');

%format du fichier .pos : SH pour scalar hexahedron
formatSpec = 'SH(%4.2f,%4.2f,%4.2f,%4.2f,%4.2f,%4.2f,%4.2f,%4.2f,%4.2f,%4.2f,%4.2f,%4.2f,%4.2f,%4.2f,%4.2f,%4.2f,%4.2f,%4.2f,%4.2f,%4.2f,%4.2f,%4.2f,%4.2f,%4.2f){%4.2f,%4.2f,%4.2f,%4.2f,%4.2f,%4.2f,%4.2f,%4.2f};\r\n';
tic
'Writting discrete field to control element size'

Lx=[0,1,1,0,0,1,1,0]*2;
Ly=[0,0,1,1,0,0,1,1]*2;
Lz=[0,0,0,0,1,1,1,1]*2;
A=zeros(1,24); B=zeros(1,8);
for x=1:2:ceil(Xmax-Xmin)+18
    x/(ceil(range(X_t_implt))+18)*100
    for y=1:2:ceil(Ymax-Ymin)+18
        for z=1:2:ceil(Zmax-Zmin)+18
            for i=1:8
                A(3*i-2:3*i)=[x+Lx(i)-10+Xmin,...
                    y+Ly(i)-10+Ymin,...
                    z+Lz(i)-10+Zmin];
                B(i)=Mat3D_final(x+Lx(i),y+Ly(i),z+Lz(i));
            end
            AB=[A,B];
            fprintf(fID2,formatSpec,AB);
        end
    end
end
toc
fprintf(fID2,'};')
fclose('all')


end

