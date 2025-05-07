//============================================================================
 // PlayterSim.cs
 // Simulation of a Playter Doll.
 //
 // This is where students write their code.
 //============================================================================
 using System;
 
 public partial class PlayterSim : Simulator
 {
     // // Parameters
     // double mA;     // dimensionless arm mass
     // double rho;    // dimless radius if gyration for moment of inertia, I_Gx
     // double gammaY; // ratio I_Gy/I_Gx
     // double gammaZ; // ratio I_Gz/I_Gx
     // double h;      // dimless vertical (b_y) distance of shoulder from body CG
     // double L;      // dimless distance of arm mass from shoulder
     // double k;      // dimless torsional stiffness of shoulder spring
     // double c;      // dimless torsional damping coeff in shoulder damper
     // double phi;    // angle of arm swing plane relative to vertical
     // double cosPhi;
     // double sinPhi;
 
     // // generalized speeds
     // double omegaX;
     // double omegaY;
     // double omegaZ;
     // double omegaFL;  // time derivative of thetaL
     // double omegaFR;  // time derivative of thetaR
     // double vx;
     // double vy;
     // double vz;
 
     // // generalized coordinates
     // double q0;      // quaternion coords
     // double q1;
     // double q2;
     // double q3;
     // double thetaL;  // left arm angle
     // double thetaR;  // right arm angle
     // double xG;      // coordinates of body's center of mass
     // double yG;
     // double zG;
 
     // Some extra stuff a student might want (feel free to define more.)
     LinSysEq sys;    // linear algebraic equation solver
     double[,] Amat;  // some arrays that might be handy 
     double[] Bmat;
 
 
 
     //------------------------------------------------------------------------
     // StudentInit: Student might want to allocate arrays and the like 
     //              before simulation begins
     //------------------------------------------------------------------------
     private void StudentInit()
     {
           sys   = new LinSysEq(8);
    Amat  = new double[8,8];
    Bmat  = new double[8];

    // --- Initialize ALL the Playter parameters! ---
    mA       = 0.107;             // arm mass ratio
    rho      = 2.21;              // radius of gyration
    gammaY   = 0.091;             // IGy/IGx
    gammaZ   = 1.05;              // IGz/IGx
    h        = 1.56;              // nondimensional shoulder height
    L        = 1.65;              // nondimensional arm length
    k        = 6;              // torsional spring stiffness (test)
    c        =  .25;              // torsional damping  (test)
    phi      = 0;         // swing‐plane angle
    sinPhi   = Math.Sin(phi);
    cosPhi   = Math.Cos(phi);
 
     }
 
     //------------------------------------------------------------------------
     // RHSFuncPlayter:  Evaluates the right sides of the differential
     //                   equations for the Playter Doll
     //------------------------------------------------------------------------
     private void RHSFuncPlayter(double[] xx, double t, double[] ff)
     {
         omegaX  = xx[0];
         omegaY  = xx[1];
         omegaZ  = xx[2];
         omegaFL = xx[3];
         omegaFR = xx[4];
         vx      = xx[5];
         vy      = xx[6];
         vz      = xx[7];
 
         q0      = xx[8];
         q1      = xx[9];
         q2      = xx[10];
         q3      = xx[11];
         thetaL  = xx[12];
         thetaR  = xx[13];
         xG      = xx[14];
         yG      = xx[15];
         zG      = xx[16];
 
double rho2 = rho*rho; 
 

                  
//ff[0] = -AngVelCrossAngMo.x/rho2;
//ff[1] = -AngVelCrossAngMo.y/(rho2*gammaY);
//ff[2] = -AngVelCrossAngMo.z/(rho2*gammaZ);

         
       

// 1) Unit basis vectors
Vex bX = new Vex(q0*q0 + q1*q1 - q2*q2 - q3*q3,
                 2*(q1*q2 + q0*q3),
                 2*(q1*q3 - q0*q2));

Vex bY = new Vex(2*(q1*q2 - q0*q3),
                 q0*q0 - q1*q1 + q2*q2 - q3*q3,
                 2*(q2*q3 + q0*q1));

Vex bZ = new Vex(2*(q1*q3 + q0*q2),
                 2*(q2*q3 - q0*q1),
                 q0*q0 - q1*q1 - q2*q2 + q3*q3);

Vex vG_vx = (q0*q0 + q1*q1 - q2*q2 - q3*q3)*bX
            +(-2*q0*q3   + 2*q1*q2)*bY
            +(2*q0*q2   + 2*q1*q3)*bZ;

Vex vG_vy = (2*q0*q3   + 2*q1*q2)*bX
            +(q0*q0 - q1*q1 + q2*q2 - q3*q3)*bY
            +(-2*q0*q1   + 2*q2*q3)*bZ;

Vex vG_vz = (-2*q0*q2   + 2*q1*q3)*bX
           +(2*q0*q1   + 2*q2*q3)*bY
           +(q0*q0 - q1*q1 - q2*q2 + q3*q3)*bZ;

Vex vG_B = (vG_vx * vx + vG_vy * vy + vG_vz * vz);



Vex rFL_SL = L * (Math.Cos(thetaL) * bX +
                  Math.Sin(thetaL) * cosPhi * bY +
                  Math.Sin(thetaL) * sinPhi * bZ);

Vex rFR_SR = -L * (Math.Cos(thetaR) * bX +
                   Math.Sin(thetaR) * cosPhi * bY +
                   Math.Sin(thetaR) * sinPhi * bZ);


Vex rSL_G   = bX + (h * bY);
Vex rSR_G   = (h * bY) -  bX ;

// 4) Full arm mass locations (B‐frame)
Vex rFL_G = rSL_G + rFL_SL;
Vex rFR_G = rSR_G + rFR_SR;

// 5) Body and hinge angular speeds
// Omega N>B
Vex omegaB = (omegaX*bX + omegaY*bY + omegaZ*bZ);
// Omega B>F (sZ)
Vex sZ = cosPhi * bZ - sinPhi * bY;
Vex sZn = -1 * sZ;
//Vex sZ = new Vex(0, -sinPhi , cosPhi);

// Omega B>F_L
Vex omegaBFL = omegaFL*cosPhi*bZ;
Vex omegaBFR = omegaFR*cosPhi*bZ;

// Omega B>F_R
Vex omegaFL_N = omegaX*bX + omegaY*bY + omegaZ*bZ + omegaFL*bZ;
Vex omegaFR_N = omegaX*bX + omegaY*bY + omegaZ*bZ + omegaFR*bZ;


//test
// Body's angular momentum
Vex H_body = new Vex(rho2 * omegaX, rho2 * gammaY * omegaY, rho2 * gammaZ * omegaZ);

// Arms' contributions
Vex H_L = mA * Vex.Cross(rFL_G, Vex.Cross(omegaB, rFL_G)) +
          mA * L * L * omegaFL * sZ;

Vex H_R = mA * Vex.Cross(rFR_G, Vex.Cross(omegaB, rFR_G)) +
          mA * L * L * omegaFR * sZ;

// Total angular momentum
Vex H_total = H_body + H_L + H_R;
Vex AngVelCrossAngMo = Vex.Cross(omegaB, H_total);


// 6) Partial velocities (Eq. 27)
//  Left arm:
Vex vFL_wx = Vex.Cross(bX, rFL_G);
Vex vFL_wy = Vex.Cross(bY, rFL_G);
Vex vFL_wz = Vex.Cross(bZ, rFL_G);
Vex vFL_wL = Vex.Cross(bZ, rFL_SL);
//Vex vFL_wL = Vex.Cross(sZn, rFL_SL);

//  Right arm:
Vex vFR_wx = Vex.Cross(bX, rFR_G);
Vex vFR_wy = Vex.Cross(bY, rFR_G);
Vex vFR_wz = Vex.Cross(bZ, rFR_G);
Vex vFR_wR = Vex.Cross(bZ, rFR_SR);
//Vex vFR_wR = Vex.Cross(sZn, rFR_SR);

// 7) Left‐arm transport/Coriolis terms
Vex term1   = Vex.Cross(omegaB, Vex.Cross(omegaB,   rSL_G));
Vex term2   = Vex.Cross(Vex.Cross(omegaB, omegaBFL), rFL_SL);
Vex term3   = Vex.Cross(omegaFL_N, Vex.Cross(omegaFL_N, rFL_SL));

double term1s  = Vex.Dot(term1, vFL_wL);
double term2s  = Vex.Dot(term2, vFL_wL);
double term3s  = Vex.Dot(term3, vFL_wL);


Vex transportSum   = term1 + term2 + term3;
double Q_L = Vex.Dot(transportSum, vFL_wL);

// 8) Right‐arm transport/Coriolis terms
Vex term1_R = Vex.Cross(omegaB, Vex.Cross(omegaB,   rSR_G));
Vex term2_R = Vex.Cross(Vex.Cross(omegaB, omegaBFR), rFR_SR);
Vex term3_R = Vex.Cross(omegaFR_N, Vex.Cross(omegaFR_N, rFR_SR));

double term1Rs  = Vex.Dot(term1_R, vFR_wR);
double term2Rs  = Vex.Dot(term2_R, vFR_wR);
double term3Rs  = Vex.Dot(term3_R, vFR_wR);

Vex transportSum_R = term1_R + term2_R + term3_R;
double Q_R = Vex.Dot(transportSum_R, vFR_wR);




//torque
double TtildeL = -mA * L * L * (k * thetaL + c * omegaFL);
double TtildeR = -mA * L * L * (k * thetaR + c * omegaFR);





// 2) Fill A (inertia/mass matrix)

// -- Body‐spin inertia (EQUATION 1)

Vex omegaCrossH = Vex.Cross( new Vex(omegaX,omegaY,omegaZ),
                              new Vex(rho2*omegaX,
                                      rho2*gammaY*omegaY,
                                      rho2*gammaZ*omegaZ) );



double mTotal = 1.0 + 2.0*mA;


// -----------------------------
// Rigid Body Inertia Terms
// -----------------------------

double A00 = rho2 +   mA * Vex.Dot(vFL_wx, vFL_wx) + mA * Vex.Dot(vFR_wx, vFR_wx) ;
double A01 =          mA * Vex.Dot(vFL_wx, vFL_wy) + mA * Vex.Dot(vFR_wx, vFR_wy);
double A02 =          mA * Vex.Dot(vFL_wx, vFL_wz) + mA * Vex.Dot(vFR_wx, vFR_wz);



double A03 =          mA * Vex.Dot(vFL_wx, vFL_wL); 
double A04 =          mA * Vex.Dot(vFR_wx, vFR_wR); 

//tests




double T1 = mA * Vex.Dot(vFL_wx, vFR_wR); 
double T2 = mA * Vex.Dot(vFR_wx, vFL_wL); 
double T3 = vFL_wL.z;
double T4 = vFR_wR.x;
double T5 = vFR_wR.y;
double T6 = vFR_wR.z;

double A05 = mA * (vFL_wx.x + vFR_wx.x);
double A06 = mA * (vFL_wx.y + vFR_wx.y);
double A07 = mA * (vFL_wx.z + vFR_wx.z);


double A10 = A01;
double A11 = rho2 * gammaY + mA * Vex.Dot(vFL_wy , vFL_wy) + mA * Vex.Dot(vFR_wy , vFR_wy) ;
double A12 = mA * Vex.Dot(vFL_wy, vFL_wz) + mA * Vex.Dot(vFR_wy, vFR_wz);       

double A13 =  mA * Vex.Dot(vFL_wy, vFL_wL);     
double A14 =  mA * Vex.Dot(vFR_wy, vFR_wR);      

double A15 = mA * (vFL_wy.x + vFR_wy.x);
double A16 = mA * (vFL_wy.y + vFR_wy.y);
double A17 = mA * (vFL_wy.z + vFR_wy.z);


double A20 = A02;
double A21 = A12;
double A22 = rho2 * gammaZ + mA * Vex.Dot(vFL_wz , vFL_wz) + mA * Vex.Dot(vFR_wz , vFR_wz);       

double A23 =  mA * Vex.Dot(vFL_wz, vFL_wL);     
double A24 =  mA * Vex.Dot(vFR_wz, vFR_wR);      

double A25 = mA * (vFL_wz.x + vFR_wz.x);
double A26 = mA * (vFL_wz.y + vFR_wz.y);
double A27 = mA * (vFL_wz.z + vFR_wz.z);


double A30 = A03;
double A31 = A13;
double A32 = A23;      

double A33 =  mA * Vex.Dot(vFL_wL, vFL_wL);     
double A34 =  0;      

double A35 = mA * vFL_wL.x;
double A36 = mA * vFL_wL.y;
double A37 = mA * vFL_wL.z;



double A40 = A04;
double A41 = A14;
double A42 = A24;      

double A43 = A34;     
double A44 =  0;      

double A45 = mA * vFR_wR.x;
double A46 = mA * vFR_wR.y;
double A47 = mA * vFR_wR.z;



double A50 = A05;
double A51 = A15;
double A52 = A25;      

double A53 = A35;     
double A54 = A45;      

double A55 = mTotal;
double A56 = 0;
double A57 = 0;



double A60 = A06;
double A61 = A16;
double A62 = A26;      

double A63 = A36;     
double A64 = A46;      

double A65 = 0;
double A66 = mTotal;
double A67 = 0;



double A70 = A07;
double A71 = A17;
double A72 = A27;      

double A73 = A37;     
double A74 = A47;      

double A75 = 0;
double A76 = 0;
double A77 = mTotal;

sys.SetA(0, 0, A00);      
sys.SetA(0, 1, A01); 
sys.SetA(0, 2, A02); 
sys.SetA(0, 3, A03); 
sys.SetA(0, 4, A04); 
sys.SetA(0, 5, A05); 
sys.SetA(0, 6, A06); 
sys.SetA(0, 7, A07);  

sys.SetA(1, 0, A10); 
sys.SetA(1, 1, A11); 
sys.SetA(1, 2, A12); 
sys.SetA(1, 3, A13); 
sys.SetA(1, 4, A14); 
sys.SetA(1, 5, A15); 
sys.SetA(1, 6, A16); 
sys.SetA(1, 7, A17); 

sys.SetA(2, 0, A20); 
sys.SetA(2, 1, A21); 
sys.SetA(2, 2, A22); 
sys.SetA(2, 3, A23); 
sys.SetA(2, 4, A24); 
sys.SetA(2, 5, A25); 
sys.SetA(2, 6, A26); 
sys.SetA(2, 7, A27); 

sys.SetA(3, 0, A30); 
sys.SetA(3, 1, A31); 
sys.SetA(3, 2, A32); 
sys.SetA(3, 3, A33); 
sys.SetA(3, 4, A34); 
sys.SetA(3, 5, A35); 
sys.SetA(3, 6, A36); 
sys.SetA(3, 7, A37); 

sys.SetA(4, 0, A40); 
sys.SetA(4, 1, A41); 
sys.SetA(4, 2, A42); 
sys.SetA(4, 3, A43); 
sys.SetA(4, 4, A44); 
sys.SetA(4, 5, A45); 
sys.SetA(4, 6, A46); 
sys.SetA(4, 7, A47);

sys.SetA(5, 0, A50); 
sys.SetA(5, 1, A51); 
sys.SetA(5, 2, A52); 
sys.SetA(5, 3, A53); 
sys.SetA(5, 4, A54); 
sys.SetA(5, 5, A55); 
sys.SetA(5, 6, A56); 
sys.SetA(5, 7, A57); 

sys.SetA(6, 0, A60); 
sys.SetA(6, 1, A61); 
sys.SetA(6, 2, A62); 
sys.SetA(6, 3, A63); 
sys.SetA(6, 4, A64); 
sys.SetA(6, 5, A65); 
sys.SetA(6, 6, A66); 
sys.SetA(6, 7, A67); 

sys.SetA(7, 0, A70); 
sys.SetA(7, 1, A71); 
sys.SetA(7, 2, A72); 
sys.SetA(7, 3, A73); 
sys.SetA(7, 4, A74); 
sys.SetA(7, 5, A75); 
sys.SetA(7, 6, A76); 
sys.SetA(7, 7, A77); 


Double B0= -omegaCrossH.x ;
Double B1= -omegaCrossH.y ;
Double B2= -omegaCrossH.z ;
Double B3= -Q_L + TtildeL ;
Double B4= -Q_R + TtildeR ;
Double B5= vx ;
Double B6= vy ;
Double B7= vz ;

sys.SetB(0, B0);
sys.SetB(1, B1);  
sys.SetB(2, B2);  
sys.SetB(3, B3);  
sys.SetB(4, B4);  
sys.SetB(5, B5);  
sys.SetB(6, B6);  
sys.SetB(7, B7);        




    
  
 



// 3) Fill B (generalized forces)

// Left‐arm P’s:
double P_Ly = Vex.Dot( Vex.Cross(bY,   rFL_G), sZ );    // row ω̇y
double P_Lz = Vex.Dot( Vex.Cross(bZ,   rFL_G), sZ );    // row ω̇z
double P_LL = Vex.Dot( Vex.Cross(sZ,   rFL_SL),sZ );    // row ω̇FL (will be zero)
// Right‐arm P’s:
double P_Ry = Vex.Dot( Vex.Cross(bY,   rFR_G), sZ );    // row ω̇y
double P_Rz = Vex.Dot( Vex.Cross(bZ,   rFR_G), sZ );    // row ω̇z
double P_RR = Vex.Dot( Vex.Cross(sZ,   rFR_SR),sZ );    // row ω̇FR (zero again)








// 4) Solve the 8×8 system
sys.SolveGauss();

ff[0]= sys.Sol(0);
ff[1]= sys.Sol(1);
ff[2]= sys.Sol(2);
ff[3]= sys.Sol(3);
ff[4]= sys.Sol(4);
ff[5]= sys.Sol(5);
ff[6]= sys.Sol(6);
ff[7]= sys.Sol(7);


       

  

ff[8] = .5*(-q1*omegaX - q2*omegaY - q3*omegaZ);
ff[9] = .5*(q0*omegaX - q3*omegaY + q2*omegaZ);
ff[10] = .5*(q3*omegaX + q0*omegaY - q1*omegaZ);
ff[11] = .5*(-q2*omegaX + q1*omegaY + q0*omegaZ);
ff[12] = omegaFL;   // θ̇L = ωFL
ff[13] = omegaFR;   // θ̇R = ωFR
ff[14] = 0;   // ẋG = vx
ff[15] = 0;   // ẏG = vy
ff[16] = 0;   // żG = vz


SetDebugVal(0, term1s); 
SetDebugVal(1, term2s);  
SetDebugVal(2, term3s);   
SetDebugVal(3, term1Rs); 
SetDebugVal(4, term2Rs);  
SetDebugVal(5, term3Rs);   
SetDebugVal(6, transportSum.x);
SetDebugVal(7, transportSum.y);
SetDebugVal(8, transportSum.z);
SetDebugVal(9, transportSum_R.x);
SetDebugVal(10, transportSum_R.y);
SetDebugVal(11, transportSum_R.z);
   
//SetDebugVal(6, thetaR);  
//SetDebugVal(7, 0);

       
        //SetDebugVal(3, omegaFL);      // Expected: oscillates around 0, damps over time
       // SetDebugVal(4, TtildeL);      // Should be negative when thetaL is positive (restoring)
        //SetDebugVal(5, Q_L);          // Transport/Coriolis term, usually smaller than TtildeL
       // SetDebugVal(6, TtildeL + Q_L);// Total net driving torque on the arm
       // SetDebugVal(7, Vex.Dot(vFL_wL, sZ));  // Should be near 1.0 (it's the partial angular velocity projection)
       // SetDebugVal(8, Vex.Dot(Vex.Cross(sZ, rFL_SL), sZ));  // Should also be near 1.0 (for PLL)


        //SetDebugVal(2,  ff[2]);  // ff[2] = ω̇z
      //  SetDebugVal(3,  ff[3]);  // ff[3] = ω̇FL
       // SetDebugVal(4,  ff[4]);  // ff[4] = ω̇FR
      //  SetDebugVal(5,  ff[5]);  // ff[0] = ω̇x
       // SetDebugVal(6,  ff[6]);  // ff[1] = ω̇y
       // SetDebugVal(7,  ff[7]);  // ff[2] = ω̇z
        //SetDebugVal(8,  ff[8]);  // ff[3] = ω̇FL
   //     SetDebugVal(9,  omegaX);  // ff[4] = ω̇FR
     //   SetDebugVal(10,  omegaY);  // ff[0] = ω̇x
       // SetDebugVal(11,  omegaZ);  // ff[1] = ω̇y
        SetDebugVal(12,  H_L.x);  // ff[2] = ω̇z
        SetDebugVal(13,  H_R.x);  // ff[3] = ω̇FL
        SetDebugVal(14,  H_body.x);  // ff[4] = ω̇FR
        SetDebugVal(15,  ff[15]);  // ff[0] = ω̇x
        SetDebugVal(16,  ff[16]);  // ff[1] = ω̇y
        







   
 
         // COMMENT THESE OUT OR REMOVE WHEN READY
         //ff[0] = ff[1] = ff[2] = 0.0;   // derivs of body angular velocities set to zero
         //ff[3] = ff[4] = 0.0;           // derivs of arm angular velocities
         //ff[5] = ff[6] = ff[7] = 0.0;   // derivs of cener of mass velocities
         //ff[8] = ff[9] = ff[10] = ff[11] = 0.0;  // derivs of quaternion coords
         //ff[12] = ff[13] = 0.0;         // derivs of arm angles
         //ff[14] = ff[15] = ff[16] = 0.0;  // derivs of CG coordinates
     }
 
 } // end class