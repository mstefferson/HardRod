ParamSt = 'Params4.txt';
ParamFile    =  importdata(ParamSt);

trial = ParamFile(1);
Nx   = ParamFile(2);
Ny   = ParamFile(3);
Nm   = ParamFile(4);
trec = ParamFile(9);
Nrec = ParamFile(19);

ConcSt = sprintf('Conc%d.txt',trial);
PoSt   = sprintf('PO%d.txt',trial);
NoSt   = sprintf('NO%d.txt',trial);
DistSt = sprintf('Dist%d.txt',trial);
AmpSt  = sprintf('Amp%d.txt',trial);

ConcFile    =  importdata(ConcSt);
PoRecFile      =  importdata(PoSt);
NoRecFile      =  importdata(NoSt);
DistFile       =  importdata(DistSt);
FTampFile      =  importdata(AmpSt);

ConcRec   = zeros(Nx,Ny,Nrec);
NoRec     = zeros(Nx,Ny,Nrec);
PoRec     = zeros(Nx,Ny,Nrec);
DistRec   = zeros(Nm,Nrec);
RhoFTRec1 = zeros(1,Nrec);
RhoFTRec2 = zeros(1,Nrec);
RhoFTRec3 = zeros(1,Nrec);
RhoFTRec4 = zeros(1,Nrec);


for t = 1:Nrec
   RhoFTRec1(t) = FTampFile(t,1) + sqrt(-1) * FTampFile(t,2);
   RhoFTRec2(t) = FTampFile(t,3) + sqrt(-1) * FTampFile(t,4);
   RhoFTRec3(t) = FTampFile(t,5) + sqrt(-1) * FTampFile(t,6);
   RhoFTRec4(t) = FTampFile(t,7) + sqrt(-1) * FTampFile(t,8);
   DistRec(:,t) = DistFile( (t-1) * Nm + 1 : (t-1) * Nm + Nm );
   
    for i = 1:Nx 
        for j = 1:Ny
            ConcRec(i,j,t) = ConcFile(i + (t-1) * Nx,j);
            PoRec(i,j,t) = PoRecFile(i + (t-1) * Nx,j);
            NoRec(i,j,t) = NoRecFile(i + (t-1) * Nx,j);
            
                       
        end
    end
end






        


