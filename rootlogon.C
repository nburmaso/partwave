void rootlogon(){
  gSystem->SetIncludePath("-I. -I$ALICE_ROOT/include -I$ALICE_PHYSICS/include -I$GEANT3DIR/TGeant3/include -I$ROOTSYS/include/Math");
  gSystem->Load("libMathMore");
  
}
