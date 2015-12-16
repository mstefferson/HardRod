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

RhoFTRec  = zeros(8,Nrec);


for t = 1:Nrec
   RhoFTRec(1,t) = FTampFile(t,1) + sqrt(-1) * FTampFile(t,2);
   RhoFTRec(2,t) = FTampFile(t,3) + sqrt(-1) * FTampFile(t,4);
   RhoFTRec(3,t) = FTampFile(t,5) + sqrt(-1) * FTampFile(t,6);
   RhoFTRec(4,t) = FTampFile(t,7) + sqrt(-1) * FTampFile(t,8);
   RhoFTRec(5,t) = FTampFile(t,9) + sqrt(-1) * FTampFile(t,10);
   RhoFTRec(6,t) = FTampFile(t,11) + sqrt(-1) * FTampFile(t,12);
   RhoFTRec(7,t) = FTampFile(t,13) + sqrt(-1) * FTampFile(t,14);
   RhoFTRec(8,t) = FTampFile(t,14) + sqrt(-1) * FTampFile(t,16);

   DistRec(:,t) = DistFile( (t-1) * Nm + 1 : (t-1) * Nm + Nm );
   
    for i = 1:Nx 
        for j = 1:Ny
            ConcRec(i,j,t) = ConcFile(i + (t-1) * Nx,j);
            PoRec(i,j,t) = PoRecFile(i + (t-1) * Nx,j);
            NoRec(i,j,t) = NoRecFile(i + (t-1) * Nx,j);
            
                       
        end
    end
end






        


