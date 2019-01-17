function [ Prosthesis, StemCenter, Thickness, size ] = SelectImplantSize(CWD, PC_ML_Width, PC_AP_Width, Prosth_Type, LongStem, LegSideName , reduceSize )
%Select prosthesis based on antero-Post & Medio-Lat widths
%



switch Prosth_Type
    case 1
        Type = {'S4','S5','S6'};
        Widths = [66.25 47; 75 47;75 61];
        Thickness = 4;
    case 2
        Type = {'S4','S5','S6','S7'};
        Widths = [67 42.8; 67 46.4;67 46.8;76 50.5;76 55.5];
        Thickness = 2.5;
        
    case 3
        Type = {'SizeC','SizeD','SizeE','SizeF','SizeG'};
        Dims = [44.9	39.5	63.8;
            47.2	41.8	67;
            50.2	44.6	71;
            53.3	47.4	75.1;
            56.5	50.2	79];
        Widths = [Dims(:,3), 0.8*Dims(:,2) + 0.2*Dims(:,1)] ;
        Thickness = 3.8;
        
        
        
end
% Choose inferior size if it's too large MedioLaterally
k = dsearchn(Widths,[PC_ML_Width PC_AP_Width]);
if Widths(k,1)-PC_ML_Width<0 && k>1
    k=k-1;
end

if nargin == 7
    warning(strcat('Implant size have been reduced: old size :',char(Type(k))));
    k=k-1;
    k=max(1,k);
    warning(strcat('Implant size have been reduced: new size :',char(Type(k))));
end

size = char(Type(k));

if Prosth_Type ==3
    size = strcat(size,'_',LegSideName);
end

if LongStem
    size = strcat(size,'_WLS');
end

% ImplantFile = py.txt2mtlb.find_ProsthFile(CWD,size,Prosth_Type);
% ImplantFile = py.txt2mtlb.find_ProsthFile(CWD,size,Prosth_Type);
% ImplantFile = py.txt2mtlb.find_ProsthFile(pyargs('directory',CWD,'Pname',size,'Ptype',Prosth_Type));
% ImplantFile = char(ImplantFile);
% XYZ = py.txt2mtlb.read_nodesGMSH(ImplantFile);
% Pts_prosthesis_init = [cell2mat(cell(XYZ{'X'}))' cell2mat(cell(XYZ{'Y'}))' cell2mat(cell(XYZ{'Z'}))'];

ProsthFileName = py.txt2mtlb.find_ProsthFile(pyargs('directory',CWD,'Pname',size,'Ptype',Prosth_Type));
ProsthFileName = char(ProsthFileName);
% Pts_prosthesis_init = [cell2mat(cell(XYZELMTS{'X'}))' cell2mat(cell(XYZELMTS{'Y'}))' cell2mat(cell(XYZELMTS{'Z'}))'];
% Elmts2D = double([cell2mat(cell(XYZELMTS{'N1'}))' cell2mat(cell(XYZELMTS{'N2'}))' cell2mat(cell(XYZELMTS{'N3'}))']);
% Elmts2D = fixNormals(Pts_prosthesis_init, Elmts2D );
% Prosthesis = triangulation(Elmts2D, Pts_prosthesis_init);

[Prosthesis] = ReadMesh(ProsthFileName);



% Thickness of the slice in distal points for stem tip center calculation
ht = 5;

switch Prosth_Type
    case 1
        [~,Izmax] = max(Prosthesis.Points(:,3));
        Zmax = Prosthesis.Points(Izmax,3);
        StemCenter = mean(Prosthesis.Points(Prosthesis.Points(:,3)>Zmax-ht,:))+...
            [0 0 ht/2];
    case 2
        [~,Izmin] = min(Prosthesis.Points(:,3));
        Zmin = Prosthesis.Points(Izmin,3);
        StemCenter = mean(Prosthesis.Points(Prosthesis.Points(:,3)<Zmin+ht,:))+...
            [0 0 -ht/2];
        
    case 3
        [~,Izmin] = min(Prosthesis.Points(:,3));
        Zmin = Prosthesis.Points(Izmin,3);
        StemCenter = mean(Prosthesis.Points(Prosthesis.Points(:,3)<Zmin+ht,:))+...
            [0 0 -ht/2];
        Centroid = mean(Prosthesis.Points);
        Centroid(3) = 0;
        Centroid(2) = 0.8*range(Prosthesis.Points(:,2));
        Prosthesis = triangulation(Prosthesis.ConnectivityList,bsxfun(@minus,Prosthesis.Points, Centroid));
end

end

