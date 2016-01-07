if( param_name.compare( "seed" ) == 0 ) {
  params.seed =  atoi( param_value.c_str() );
  std::cout << " seed = " << param_value << std::endl;
}

if( param_name.compare( "trial" ) == 0 ) {
  params.trial =  atoi( param_value.c_str() );
  std::cout << " trial = " << param_value << std::endl;
}

if( param_name.compare( "Nx" ) == 0 ) {
  params.Nx =  atoi( param_value.c_str() );
  std::cout << " Nx = " << param_value << std::endl;
}

if( param_name.compare( "Ny" ) == 0 ) {
  params.Ny =  atoi( param_value.c_str() );
  std::cout << " Ny =  " << param_value << std::endl;
}

if( param_name.compare( "Nm" ) == 0 ) {
  params.Nm =  atoi( param_value.c_str() );
  std::cout << " Nm = " << param_value << std::endl;
}

if( param_name.compare( "StepFlag" ) == 0 ) {
  params.StepFlag =  atoi( param_value.c_str() );
  std::cout << " StepFlag = " << param_value << std::endl;
}

if( param_name.compare( "IsoDiffFlag" ) == 0 ) {
  params.IsoDiffFlag =  atoi( param_value.c_str() );
  std::cout << " IsoDiffFlag = " << param_value << std::endl;
}

if( param_name.compare( "IcFlag" ) == 0 ) {
  params.IcFlag =  atoi( param_value.c_str() );
  std::cout << " IcFlag = " << param_value << std::endl;
}

if( param_name.compare( "RandPerbAmpFlag" ) == 0 ) {
  params.RandPerbAmpFlag =  atoi( param_value.c_str() );
  std::cout << " RandPerbAmpFlag = " << param_value << std::endl;
}

if( param_name.compare( "wrtRhoFlag" ) == 0 ) {
  params.wrtRhoFlag =  atoi( param_value.c_str() );
  std::cout << " wrtRhoFlag = " << param_value << std::endl;
}

if( param_name.compare( "wrtOpFlag" ) == 0 ) {
  params.wrtOpFlag =  atoi( param_value.c_str() );
  std::cout << " wrtOpFlag = " << param_value << std::endl;
}

if( param_name.compare( "pXmodes" ) == 0 ) {
  params.pXmodes =  atoi( param_value.c_str() );
  std::cout << " pXmodes = " << param_value << std::endl;
}

if( param_name.compare( "pYmodes" ) == 0 ) {
  params.pYmodes =  atoi( param_value.c_str() );
  std::cout << " pYmodes = " << param_value << std::endl;
}

if( param_name.compare( "pMmodes" ) == 0 ) {
  params.pMmodes =  atoi( param_value.c_str() );
  std::cout << " pMmodes = " << param_value << std::endl;
}

if( param_name.compare( "dt" ) == 0 ) {
  params.dt =  atof( param_value.c_str() );
  std::cout << " dt = " << param_value << std::endl;
}

if( param_name.compare( "trec" ) == 0 ) {
  params.trec =  atof( param_value.c_str() );
  std::cout << " trec = " << param_value << std::endl;
}

if( param_name.compare( "tend" ) == 0 ) {
  params.tend =  atof( param_value.c_str() );
  std::cout << " tend = " << param_value << std::endl;
}

if( param_name.compare( "bc" ) == 0 ) {
  params.bc =  atof( param_value.c_str() );
  std::cout << " bc = " << param_value << std::endl;
}

if( param_name.compare( "vD" ) == 0 ) {
  params.vD =  atof( param_value.c_str() );
  std::cout << " vD = " << param_value << std::endl;
}

if( param_name.compare( "Lx" ) == 0 ) {
  params.Lx =  atof( param_value.c_str() );
  std::cout << " Lx = " << param_value << std::endl;
}

if( param_name.compare( "Ly" ) == 0 ) {
  params.Ly =  atof( param_value.c_str() );
  std::cout << " Ly = " << param_value << std::endl;
}

if( param_name.compare( "Lrod" ) == 0 ) {
  params.Lrod =  atof( param_value.c_str() );
  std::cout << " Lrod = " << param_value << std::endl;
}

if( param_name.compare( "Dpar" ) == 0 ) {
  params.Dpar =  atof( param_value.c_str() );
  std::cout << " Dpar = " << param_value << std::endl;
}

if( param_name.compare( "Dperp" ) == 0 ) {
  params.Dperp =  atof( param_value.c_str() );
  std::cout << " Dperp = " << param_value << std::endl;
}

if( param_name.compare( "Dr" ) == 0 ) {
  params.Dr =  atof( param_value.c_str() );
  std::cout << " Dr = " << param_value << std::endl;
}

if( param_name.compare( "PertrbAmp" ) == 0 ) {
  params.PertrbAmp =  atof( param_value.c_str() );
  std::cout << " PertrbAmp = " << param_value << std::endl;
}


