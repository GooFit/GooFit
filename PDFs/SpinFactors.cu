/*
04/05/2016 Christoph Hasse
The Spinfactors are an adaptation from the MINT implementation, by Jonas Rademacker.

DISCLAIMER:
This code is not sufficently tested yet and still under heavy development!
*/
#include "SpinFactors.hh"
#include "SpinHelper.hh"

#define ZEMACH 1

EXEC_TARGET fptype DtoV1V2_V1toP1P2_V2toP3P4_S (fptype* Vecs, unsigned int* indices) {
  unsigned int p1          = indices[2];
  unsigned int p2          = indices[3];
  unsigned int p3          = indices[4];
  unsigned int p4          = indices[5];
  gpuLVec P1(Vecs[0 + 4*p1], Vecs[1 + 4*p1], Vecs[2 + 4*p1], Vecs[3 + 4*p1]);
  gpuLVec P2(Vecs[0 + 4*p2], Vecs[1 + 4*p2], Vecs[2 + 4*p2], Vecs[3 + 4*p2]);
  gpuLVec P3(Vecs[0 + 4*p3], Vecs[1 + 4*p3], Vecs[2 + 4*p3], Vecs[3 + 4*p3]);
  gpuLVec P4(Vecs[0 + 4*p4], Vecs[1 + 4*p4], Vecs[2 + 4*p4], Vecs[3 + 4*p4]);

  // printf("vec%i %.5g, %.5g, %.5g, %.5g\n",0, P1.getX(), P1.getY(), P1.getZ(),P1.getE());
  // printf("vec%i %.5g, %.5g, %.5g, %.5g\n",1, P2.getX(), P2.getY(), P2.getZ(),P2.getE());
  // printf("vec%i %.5g, %.5g, %.5g, %.5g\n",2, P3.getX(), P3.getY(), P3.getZ(),P3.getE());
  // printf("vec%i %.5g, %.5g, %.5g, %.5g\n",3, P4.getX(), P4.getY(), P4.getZ(),P4.getE());

  gpuLVec pV1 = P1 + P2;
  gpuLVec qV1 = P1 - P2;
  gpuLVec pV2 = P3 + P4;
  gpuLVec qV2 = P3 - P4;
  
  fptype MV1 = SQRT(pV1.Dot(pV1));
  fptype MV2 = SQRT(pV2.Dot(pV2));

  fptype returnVal = (qV1.Dot(qV2) 
                   - qV1.Dot(pV1) * pV1.Dot(qV2) / (MV1*MV1)
                   - qV1.Dot(pV2) * pV2.Dot(qV2) / (MV2*MV2)
                   + qV1.Dot(pV1) * pV1.Dot(pV2) * pV2.Dot(qV2) 
                   / (MV1*MV1 * MV2*MV2));
  // printf("s1 %.5g; %i,%i,%i,%i\n",returnVal, indices[2], indices[3], indices[4], indices[5]);
  return returnVal;
}

EXEC_TARGET fptype DtoV1V2_V1toP1P2_V2toP3P4_P (fptype* Vecs, unsigned int* indices) {
  unsigned int p1          = indices[2];
  unsigned int p2          = indices[3];
  unsigned int p3          = indices[4];
  unsigned int p4          = indices[5];
  gpuLVec P1(Vecs[0 + 4*p1], Vecs[1 + 4*p1], Vecs[2 + 4*p1], Vecs[3 + 4*p1]);
  gpuLVec P2(Vecs[0 + 4*p2], Vecs[1 + 4*p2], Vecs[2 + 4*p2], Vecs[3 + 4*p2]);
  gpuLVec P3(Vecs[0 + 4*p3], Vecs[1 + 4*p3], Vecs[2 + 4*p3], Vecs[3 + 4*p3]);
  gpuLVec P4(Vecs[0 + 4*p4], Vecs[1 + 4*p4], Vecs[2 + 4*p4], Vecs[3 + 4*p4]);

  gpuLVec pV1 = P1 + P2;
  gpuLVec qV1 = P1 - P2;
  gpuLVec pV2 = P3 + P4;
  gpuLVec qV2 = P3 - P4;
  
  gpuLVec pD = pV1 + pV2;
  gpuLVec qD = pV1 - pV2;

  if (ZEMACH) {
    fptype MV1 = SQRT(pV1.Dot(pV1));
    fptype MV2 = SQRT(pV2.Dot(pV2));
    
    ZTspin1 LD(qD,pD,pD.M());
    ZTspin1 LV1(qV1,pV1,MV1);
    ZTspin1 LV2(qV2,pV2,MV2);
    
    return LeviCivita(pD,LD,LV1).Dot(LV2);
  }  

  return LeviCevita(pD, qD, qV1, qV2);
}

EXEC_TARGET fptype DtoV1V2_V1toP1P2_V2toP3P4_D (fptype* Vecs, unsigned int* indices) {
  unsigned int p1          = indices[2];
  unsigned int p2          = indices[3];
  unsigned int p3          = indices[4];
  unsigned int p4          = indices[5];
  gpuLVec P1(Vecs[0 + 4*p1], Vecs[1 + 4*p1], Vecs[2 + 4*p1], Vecs[3 + 4*p1]);
  gpuLVec P2(Vecs[0 + 4*p2], Vecs[1 + 4*p2], Vecs[2 + 4*p2], Vecs[3 + 4*p2]);
  gpuLVec P3(Vecs[0 + 4*p3], Vecs[1 + 4*p3], Vecs[2 + 4*p3], Vecs[3 + 4*p3]);
  gpuLVec P4(Vecs[0 + 4*p4], Vecs[1 + 4*p4], Vecs[2 + 4*p4], Vecs[3 + 4*p4]);

  gpuLVec pV1 = P1 + P2;
  gpuLVec qV1 = P1 - P2;
  gpuLVec pV2 = P3 + P4;
  gpuLVec qV2 = P3 - P4;


  if (ZEMACH) {
    gpuLVec pD = pV1 + pV2;
    gpuLVec qD = pV1 - pV2;
    double mD = pD.M();
    
    ZTspin1 tV1(qV1, pV1, pV1.M());
    ZTspin1 tV2(qV2, pV2, pV2.M());
    ZTspin2 tD(qD, pD, mD);
    
    return tV1.Contract(tD.Contract(tV2));
  } 
  
  // printf("%f, %f, %f, %f\n",P1.getX(), P1.getY(), P1.getZ(), P1.getE() );
  // printf("%f, %f, %f, %f\n",P2.getX(), P2.getY(), P2.getZ(), P2.getE() );
  // printf("%f, %f, %f, %f\n",P3.getX(), P3.getY(), P3.getZ(), P3.getE() );
  // printf("%f, %f, %f, %f\n",P4.getX(), P4.getY(), P4.getZ(), P4.getE() );

  fptype MV1 = SQRT(pV1.Dot(pV1));
  fptype MV2 = SQRT(pV2.Dot(pV2));
  fptype returnVal = (  qV1.Dot(pV2) - qV1.Dot(pV1) * pV1.Dot(pV2)/(MV1*MV1)
                     )*( 
                     qV2.Dot(pV1) - qV2.Dot(pV2) * pV2.Dot(pV1)/(MV2*MV2)
                     );
  return returnVal;
}

EXEC_TARGET fptype DtoV1P1_V1toV2P2_V2toP3P4 (fptype* Vecs, unsigned int* indices) {
  unsigned int p1          = indices[2];
  unsigned int p2          = indices[3];
  unsigned int p3          = indices[4];
  unsigned int p4          = indices[5];
  gpuLVec P1(Vecs[0 + 4*p1], Vecs[1 + 4*p1], Vecs[2 + 4*p1], Vecs[3 + 4*p1]);
  gpuLVec P2(Vecs[0 + 4*p2], Vecs[1 + 4*p2], Vecs[2 + 4*p2], Vecs[3 + 4*p2]);
  gpuLVec P3(Vecs[0 + 4*p3], Vecs[1 + 4*p3], Vecs[2 + 4*p3], Vecs[3 + 4*p3]);
  gpuLVec P4(Vecs[0 + 4*p4], Vecs[1 + 4*p4], Vecs[2 + 4*p4], Vecs[3 + 4*p4]);

  gpuLVec pV1 = P2 + P3 + P4;
  gpuLVec pV2 = P3 + P4;
  gpuLVec qV1 = (P3 + P4) - P2;
  gpuLVec qV2 = P3 - P4;

  if (ZEMACH) {
    double MV1 = pV1.M();
    double MV2 = pV2.M();
    
    gpuLVec pD = pV1 + P1;
    gpuLVec qD = pV1 - P1;
    
    ZTspin1 LD(qD,pD,pD.M());
    ZTspin1 LV1(qV1,pV1,MV1);
    ZTspin1 LV2(qV2,pV2,MV2);
    SpinSumV PV1(pV1,MV1);
    
    gpuLVec tmp = PV1.Dot(LD);

    return LeviCivita(pV1,LV1,LV2).Dot(tmp);
  }   
  
  // printf("%f, %f, %f, %f\n",P1.getX(), P1.getY(), P1.getZ(), P1.getE() );
  // printf("%f, %f, %f, %f\n",P2.getX(), P2.getY(), P2.getZ(), P2.getE() );
  // printf("%f, %f, %f, %f\n",P3.getX(), P3.getY(), P3.getZ(), P3.getE() );
  // printf("%f, %f, %f, %f\n",P4.getX(), P4.getY(), P4.getZ(), P4.getE() );

  fptype returnVal = LeviCevita(pV1, qV1, P1, qV2);
  return returnVal;
}

EXEC_TARGET fptype DtoVS_VtoP1P2_StoP3P4 (fptype* Vecs, unsigned int* indices) {
  unsigned int p1          = indices[2];
  unsigned int p2          = indices[3];
  unsigned int p3          = indices[4];
  unsigned int p4          = indices[5];
  gpuLVec P1(Vecs[0 + 4*p1], Vecs[1 + 4*p1], Vecs[2 + 4*p1], Vecs[3 + 4*p1]);
  gpuLVec P2(Vecs[0 + 4*p2], Vecs[1 + 4*p2], Vecs[2 + 4*p2], Vecs[3 + 4*p2]);
  gpuLVec P3(Vecs[0 + 4*p3], Vecs[1 + 4*p3], Vecs[2 + 4*p3], Vecs[3 + 4*p3]);
  gpuLVec P4(Vecs[0 + 4*p4], Vecs[1 + 4*p4], Vecs[2 + 4*p4], Vecs[3 + 4*p4]);

  gpuLVec pS =  P3 + P4;
  gpuLVec pV =  P1 + P2;
  gpuLVec qV =  P1 - P2;
  

  if(ZEMACH){
      gpuLVec pD= pV + pS;
      gpuLVec qD= pV - pS;
      ZTspin1 LD(qD,pD,pD.M());
      ZTspin1 LV(qV,pV,pV.M());
      return (LD.Dot(LV));
  }

  // printf("%f, %f, %f, %f\n",P1.getX(), P1.getY(), P1.getZ(), P1.getE() );
  // printf("%f, %f, %f, %f\n",P2.getX(), P2.getY(), P2.getZ(), P2.getE() );
  // printf("%f, %f, %f, %f\n",P3.getX(), P3.getY(), P3.getZ(), P3.getE() );
  // printf("%f, %f, %f, %f\n",P4.getX(), P4.getY(), P4.getZ(), P4.getE() );

  fptype MV = SQRT(pV.Dot(pV));

  fptype returnVal = (pS.Dot(qV) - pS.Dot(pV) * pV.Dot(qV) / (MV*MV));
  return returnVal;
}

EXEC_TARGET fptype DtoAP1_AtoSP2_StoP3P4 (fptype* Vecs, unsigned int* indices) {
  unsigned int p1          = indices[2];
  unsigned int p2          = indices[3];
  unsigned int p3          = indices[4];
  unsigned int p4          = indices[5];
  gpuLVec P1(Vecs[0 + 4*p1], Vecs[1 + 4*p1], Vecs[2 + 4*p1], Vecs[3 + 4*p1]);
  gpuLVec P2(Vecs[0 + 4*p2], Vecs[1 + 4*p2], Vecs[2 + 4*p2], Vecs[3 + 4*p2]);
  gpuLVec P3(Vecs[0 + 4*p3], Vecs[1 + 4*p3], Vecs[2 + 4*p3], Vecs[3 + 4*p3]);
  gpuLVec P4(Vecs[0 + 4*p4], Vecs[1 + 4*p4], Vecs[2 + 4*p4], Vecs[3 + 4*p4]);

  gpuLVec pS =  P3 + P4;
  gpuLVec pA =  P2 + pS;
  gpuLVec qA =  pS - P2;
  gpuLVec pD =  pA + P1;
  gpuLVec qD =  pA - P1;

  // printf("%f, %f, %f, %f\n",P1.GetX(), P1.GetY(), P1.GetZ(), P1.GetE() );
  // printf("%f, %f, %f, %f\n",P2.GetX(), P2.GetY(), P2.GetZ(), P2.GetE() );
  // printf("%f, %f, %f, %f\n",P3.GetX(), P3.GetY(), P3.GetZ(), P3.GetE() );
  // printf("%f, %f, %f, %f\n",P4.GetX(), P4.GetY(), P4.GetZ(), P4.GetE() );

  if(ZEMACH){
    ZTspin1 LD(qD,pD,pD.M());
    ZTspin1 LA(qA,pA,pA.M());
    return (LD.Dot(LA));
  } 

  fptype MA = SQRT(pA.Dot(pA));
  fptype returnVal = (P1.Dot(qA) - P1.Dot(pA) * pA.Dot(qA) / (MA*MA));
  return returnVal;
}

EXEC_TARGET fptype DtoAP1_AtoVP2_VtoP3P4 (fptype* Vecs, unsigned int* indices) {
  unsigned int p1          = indices[2];
  unsigned int p2          = indices[3];
  unsigned int p3          = indices[4];
  unsigned int p4          = indices[5];
  gpuLVec P1(Vecs[0 + 4*p1], Vecs[1 + 4*p1], Vecs[2 + 4*p1], Vecs[3 + 4*p1]);
  gpuLVec P2(Vecs[0 + 4*p2], Vecs[1 + 4*p2], Vecs[2 + 4*p2], Vecs[3 + 4*p2]);
  gpuLVec P3(Vecs[0 + 4*p3], Vecs[1 + 4*p3], Vecs[2 + 4*p3], Vecs[3 + 4*p3]);
  gpuLVec P4(Vecs[0 + 4*p4], Vecs[1 + 4*p4], Vecs[2 + 4*p4], Vecs[3 + 4*p4]);

  gpuLVec pV = P3 + P4;
  gpuLVec qV = P3 - P4;
  gpuLVec pA = P2 + pV;
  gpuLVec p0 = P1;  
  gpuLVec pD = P1 + pA;
  gpuLVec qD = pA - P1;
  
  if(ZEMACH){  
      ZTspin1 LB(qD,pD,pD.M());
      ZTspin1 LV(qV,pV,pV.M());
      SpinSumV PA(pA,pA.M());  
      
      gpuLVec tmp= PA.Dot(LV);
      return (LB.Dot(tmp));
  }

  fptype MA = SQRT(pA.Dot(pA));
  fptype MV = SQRT(pV.Dot(pV));
  fptype returnVal =  P1.Dot(qV)
      -   p0.Dot(pA) * pA.Dot(qV) / (MA*MA)
      -   p0.Dot(pV) * pV.Dot(qV) / (MV*MV)
      +   p0.Dot(pA) * pA.Dot(pV) * pV.Dot(qV) / (MA*MA * MV*MV);
  // printf("spin %.7g\n",returnVal );
  return returnVal;
}

EXEC_TARGET fptype DtoAP1_AtoVP2Dwave_VtoP3P4 (fptype* Vecs, unsigned int* indices) {
  unsigned int p1          = indices[2];
  unsigned int p2          = indices[3];
  unsigned int p3          = indices[4];
  unsigned int p4          = indices[5];
  gpuLVec P1(Vecs[0 + 4*p1], Vecs[1 + 4*p1], Vecs[2 + 4*p1], Vecs[3 + 4*p1]);
  gpuLVec P2(Vecs[0 + 4*p2], Vecs[1 + 4*p2], Vecs[2 + 4*p2], Vecs[3 + 4*p2]);
  gpuLVec P3(Vecs[0 + 4*p3], Vecs[1 + 4*p3], Vecs[2 + 4*p3], Vecs[3 + 4*p3]);
  gpuLVec P4(Vecs[0 + 4*p4], Vecs[1 + 4*p4], Vecs[2 + 4*p4], Vecs[3 + 4*p4]);

  gpuLVec pV = P3 + P4;
  gpuLVec qV = P3 - P4;
  gpuLVec pA = P2 + pV;
  gpuLVec qA = P2 - pV;
  gpuLVec pD = P1 + pA;
  gpuLVec qD = pA - P1;

  if(ZEMACH){  
    ZTspin1 LD(qD,pD,pD.M());
    ZTspin1 LV(qV,pV,pV.M());
    ZTspin2 LA(qA,pA,pA.M()); 
    gpuLVec tmp= LA.Contract(LV);

    return (LD.Dot(tmp));
   }  
  
  fptype MA = SQRT(pA.Dot(pA));
  fptype MV = SQRT(pV.Dot(pV));

  fptype returnVal = P1.Dot((qA - qA.Dot(pA)*pA * (1./(MA*MA)))) * (qV - qV.Dot(pV)*pV * (1./(MV*MV))).Dot(pA);

  // printf("spin %.7g\n",returnVal );
  return returnVal;
}



MEM_DEVICE spin_function_ptr ptr_to_DtoV1V2_V1toP1P2_V2toP3P4_S = DtoV1V2_V1toP1P2_V2toP3P4_S;
MEM_DEVICE spin_function_ptr ptr_to_DtoV1V2_V1toP1P2_V2toP3P4_P = DtoV1V2_V1toP1P2_V2toP3P4_P;
MEM_DEVICE spin_function_ptr ptr_to_DtoV1V2_V1toP1P2_V2toP3P4_D = DtoV1V2_V1toP1P2_V2toP3P4_D;
MEM_DEVICE spin_function_ptr ptr_to_DtoVS_VtoP1P2_StoP3P4       = DtoVS_VtoP1P2_StoP3P4;
MEM_DEVICE spin_function_ptr ptr_to_DtoV1P1_V1toV2P2_V2toP3P4   = DtoV1P1_V1toV2P2_V2toP3P4;
MEM_DEVICE spin_function_ptr ptr_to_DtoAP1_AtoSP2_StoP3P4       = DtoAP1_AtoSP2_StoP3P4;
MEM_DEVICE spin_function_ptr ptr_to_DtoAP1_AtoVP2_VtoP3P4       = DtoAP1_AtoVP2_VtoP3P4;
MEM_DEVICE spin_function_ptr ptr_to_DtoAP1_AtoVP2Dwave_VtoP3P4  = DtoAP1_AtoVP2Dwave_VtoP3P4;




SpinFactor::SpinFactor (std::string name, SF_4Body SF, unsigned int P0, unsigned int P1, unsigned int P2, unsigned int P3)
 : GooPdf(0,name), _SF(SF), _P0(P0), _P1(P1), _P2(P2), _P3(P3){
  vector<unsigned int> pindices; 
  pindices.push_back(0); //dummy for index to constants.
  pindices.push_back(P0);
  pindices.push_back(P1);
  pindices.push_back(P2);
  pindices.push_back(P3);
  switch(SF){
    case SF_4Body::DtoV1V2_V1toP1P2_V2toP3P4_S:
      GET_FUNCTION_ADDR(ptr_to_DtoV1V2_V1toP1P2_V2toP3P4_S);
      break;
    case SF_4Body::DtoV1V2_V1toP1P2_V2toP3P4_P:
      GET_FUNCTION_ADDR(ptr_to_DtoV1V2_V1toP1P2_V2toP3P4_P);
      break;
    case SF_4Body::DtoV1V2_V1toP1P2_V2toP3P4_D:
      GET_FUNCTION_ADDR(ptr_to_DtoV1V2_V1toP1P2_V2toP3P4_D);
      break;
    case SF_4Body::DtoAP1_AtoVP2_VtoP3P4:
      GET_FUNCTION_ADDR(ptr_to_DtoAP1_AtoVP2_VtoP3P4);
      break;
    case SF_4Body::DtoAP1_AtoVP2Dwave_VtoP3P4:
      GET_FUNCTION_ADDR(ptr_to_DtoAP1_AtoVP2Dwave_VtoP3P4);
      break;
    case SF_4Body::DtoVS_VtoP1P2_StoP3P4:
      GET_FUNCTION_ADDR(ptr_to_DtoVS_VtoP1P2_StoP3P4);
      break;
    case SF_4Body::DtoV1P1_V1toV2P2_V2toP3P4:
      GET_FUNCTION_ADDR(ptr_to_DtoV1P1_V1toV2P2_V2toP3P4);
      break;    
    case SF_4Body::DtoAP1_AtoSP2_StoP3P4:
      GET_FUNCTION_ADDR(ptr_to_DtoAP1_AtoSP2_StoP3P4);
      break;
        
    default:
      std::cout << "No Spinfunction implemented for that kind." << std::endl;
      exit(0);
      break;
  }
  
  initialise(pindices);
}

