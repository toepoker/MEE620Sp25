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
 
         Vex AngVel = new Vex(omegaX, omegaY,omegaZ);
         Vex AngMo = new Vex(rho2*omegaX, rho2*gammaY*omegaY, rho2*gammaZ*omegaZ);
         Vex AngVelCrossAngMo = Vex.Cross(AngVel, AngMo);
         
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

Vex vG_vx = new Vex(
            (q0*q0 + q1*q1 - q2*q2 - q3*q3)*bX
            +(-2*q0*q3   + 2*q1*q2)*bY
            +(2*q0*q2   + 2*q1*q3)*bZ
        );
Vex vG_vy = new Vex(
            (2*q0*q3   + 2*q1*q2)*bX
            +(q0*q0 - q1*q1 + q2*q2 - q3*q3)*bY
            +(-2*q0*q1   + 2*q2*q3)*bZ
        );
Vex vG_vz = new Vex(
           (-2*q0*q2   + 2*q1*q3)*bX
           +(2*q0*q1   + 2*q2*q3)*bY
           +(q0*q0 - q1*q1 - q2*q2 + q3*q3)*bZ
        );

Vex vG_B = vG_vx * vx + vG_vy * vy + vG_vz * vz;

// 2) Hinge‐axis unit in B‐frame
Vex sZ = new Vex(0, -sinPhi, cosPhi )  ;

Vex rFL_SL = L * (Math.Cos(thetaL) * bX +
                  Math.Sin(thetaL) * cosPhi * bY +
                  Math.Sin(thetaL) * sinPhi * bZ);

Vex rFR_SR = -L * (Math.Cos(thetaR) * bX +
                   Math.Sin(thetaR) * cosPhi * bY +
                   Math.Sin(thetaR) * sinPhi * bZ);


Vex rSL_G   = bX + (h * bY);
Vex rSR_G   = (h * bY) -  bX ;

// 3) Shoulder & arm position vectors (B‐frame)
//Vex rSL_G   = new Vex( 1.0, h, 0.0 );
//Vex rSR_G   = new Vex(-1.0, h, 0.0 );
//Vex rFL_SL  = L * new Vex(Math.Cos(thetaL),
//                         Math.Sin(thetaL)*cosPhi,
 //                        Math.Sin(thetaL)*sinPhi);
//Vex rFR_SR  = -L * new Vex(Math.Cos(thetaR),
  //                        Math.Sin(thetaR)*cosPhi,
    //                      Math.Sin(thetaR)*sinPhi);

// 4) Full arm mass locations (B‐frame)
Vex rFL_G = rSL_G + rFL_SL;
Vex rFR_G = rSR_G + rFR_SR;

// 5) Body and hinge angular speeds
Vex omegaB      = new Vex(omegaX, omegaY, omegaZ);
//Vex omegaFL_vec = sZ * omegaFL;
//Vex omegaFR_vec = sZ * omegaFR;

//test
Vex omegaFL_vec = omegaB + sZ * omegaFL;
Vex omegaFR_vec = omegaB + sZ * omegaFR;




// 6) Partial velocities (Eq. 27)
//  Left arm:
Vex vFL_wx = Vex.Cross(bX, rFL_G);
Vex vFL_wy = Vex.Cross(bY, rFL_G);
Vex vFL_wz = Vex.Cross(bZ, rFL_G);
Vex vFL_wL = Vex.Cross(sZ, rFL_SL);
//  Right arm:
Vex vFR_wx = Vex.Cross(bX, rFR_G);
Vex vFR_wy = Vex.Cross(bY, rFR_G);
Vex vFR_wz = Vex.Cross(bZ, rFR_G);
Vex vFR_wR = Vex.Cross(sZ, rFR_SR);

// 7) Left‐arm transport/Coriolis terms
Vex term1   = Vex.Cross(omegaB, Vex.Cross(omegaB,   rSL_G));
Vex term2   = Vex.Cross(Vex.Cross(omegaB, omegaFL_vec), rFL_SL);
Vex term3   = Vex.Cross(omegaFL_vec, Vex.Cross(omegaFL_vec, rFL_SL));
Vex transportSum   = term1 + term2 + term3;
double Q_L = Vex.Dot(transportSum, vFL_wL);

// 8) Right‐arm transport/Coriolis terms
Vex term1_R = Vex.Cross(omegaB, Vex.Cross(omegaB,   rSR_G));
Vex term2_R = Vex.Cross(Vex.Cross(omegaB, omegaFR_vec), rFR_SR);
Vex term3_R = Vex.Cross(omegaFR_vec, Vex.Cross(omegaFR_vec, rFR_SR));
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

sys.SetA(0, 0,   (rho2) )           ; // I_Gx * ω̇x
sys.SetA(1, 1,   (rho2 * gammaY) )   ;// I_Gy * ω̇y 
sys.SetA(2, 2,   (rho2 * gammaZ) );// I_Gz * ω̇z

// -- Left hinge inertia (row 3) --
double mArm = mA;
sys.SetA(3, 0,  mArm * Vex.Dot(vFL_wx, vFL_wL));   // coupling ωx → ω̇FL
sys.SetA(3, 1,  mArm * Vex.Dot(vFL_wy, vFL_wL));   // coupling ωy → ω̇FL
sys.SetA(3, 2,  mArm * Vex.Dot(vFL_wz, vFL_wL));   // coupling ωz → ω̇FL
sys.SetA(3, 3,  mArm * Vex.Dot(vFL_wL, vFL_wL));   // hinge inertia term
sys.SetA(3, 4,  0.0);
sys.SetA(3, 5,  mArm * Vex.Dot(vG_vx,  vFL_wL));   // coupling vx → ω̇FL
sys.SetA(3, 6,  mArm * Vex.Dot(vG_vy,  vFL_wL));   // coupling vy → ω̇FL
sys.SetA(3, 7,  mArm * Vex.Dot(vG_vz,  vFL_wL));   // coupling vz → ω̇FL

// row 4: ω̇FR equation
sys.SetA(4, 0,  mArm * Vex.Dot(vFR_wx, vFR_wR));
sys.SetA(4, 1,  mArm * Vex.Dot(vFR_wy, vFR_wR));
sys.SetA(4, 2,  mArm * Vex.Dot(vFR_wz, vFR_wR));
sys.SetA(4, 3,  0.0);
sys.SetA(4, 4,  mArm * Vex.Dot(vFR_wR, vFR_wR));
sys.SetA(4, 5,  mArm * Vex.Dot(vG_vx,  vFR_wR));
sys.SetA(4, 6,  mArm * Vex.Dot(vG_vy,  vFR_wR));
sys.SetA(4, 7,  mArm * Vex.Dot(vG_vz,  vFR_wR));

// -- CG (rows 5–7) --
double mTotal = 1.0 + 2.0*mA;
sys.SetA(5, 5,    mTotal     ); // m_total * v̇x
sys.SetA(6, 6,    mTotal    ); // m_total * v̇y
sys.SetA(7, 7,    mTotal    ); // m_total * v̇z

// 3) Fill B (generalized forces)

// Left‐arm P’s:
double P_Ly = Vex.Dot( Vex.Cross(bY,   rFL_G), sZ );    // row ω̇y
double P_Lz = Vex.Dot( Vex.Cross(bZ,   rFL_G), sZ );    // row ω̇z
double P_LL = Vex.Dot( Vex.Cross(sZ,   rFL_SL),sZ );    // row ω̇FL (will be zero)
// Right‐arm P’s:
double P_Ry = Vex.Dot( Vex.Cross(bY,   rFR_G), sZ );    // row ω̇y
double P_Rz = Vex.Dot( Vex.Cross(bZ,   rFR_G), sZ );    // row ω̇z
double P_RR = Vex.Dot( Vex.Cross(sZ,   rFR_SR),sZ );    // row ω̇FR (zero again)




sys.SetB(0, -omegaCrossH.x );  // Q_x = −(ω×H)_x
sys.SetB(1, -omegaCrossH.y );  // Q_y = −(ω×H)_y
sys.SetB(2, -omegaCrossH.z );  // Q_z = −(ω×H)_z
sys.SetB(3,  -Q_L + TtildeL );        // Q_L = spring/damper + transport projection
sys.SetB(4, -Q_R + TtildeR );        // Q_R
// -- No external CG force (rows 5–7) --
sys.SetB(5, 0.0);  // Q_vx
sys.SetB(6, 0.0);  // Q_vy
sys.SetB(7, 0.0);  // Q_vz









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


       

      //  ff[3] = vFL_wL.x
       // ff[4] = vFL_wR.x
     //   ff[5] = vG_vx.x
       // ff[6] = vG_vy.y
        //ff[7] = vG_vz.z

        ff[8] = .5*(-q1*omegaX - q2*omegaY - q3*omegaZ);
        ff[9] = .5*(q0*omegaX - q3*omegaY + q2*omegaZ);
        ff[10] = .5*(q3*omegaX + q0*omegaY - q1*omegaZ);
        ff[11] = .5*(-q2*omegaX + q1*omegaY + q0*omegaZ);
        ff[12] = omegaFL;   // θ̇L = ωFL
        ff[13] = omegaFR;   // θ̇R = ωFR
        ff[14] = 0;   // ẋG = vx
        ff[15] = 0;   // ẏG = vy
        ff[16] = 0;   // żG = vz

        SetDebugVal(0,  ff[0]);  // ff[0] = ω̇x
        SetDebugVal(1,  ff[1]);  // ff[1] = ω̇y
        SetDebugVal(2, thetaL);       // Expected: starts at 0, swings up, dampens back down
        SetDebugVal(3, omegaFL);      // Expected: oscillates around 0, damps over time
        SetDebugVal(4, TtildeL);      // Should be negative when thetaL is positive (restoring)
        SetDebugVal(5, Q_L);          // Transport/Coriolis term, usually smaller than TtildeL
        SetDebugVal(6, TtildeL + Q_L);// Total net driving torque on the arm
        SetDebugVal(7, Vex.Dot(vFL_wL, sZ));  // Should be near 1.0 (it's the partial angular velocity projection)
        SetDebugVal(8, Vex.Dot(Vex.Cross(sZ, rFL_SL), sZ));  // Should also be near 1.0 (for PLL)


        //SetDebugVal(2,  ff[2]);  // ff[2] = ω̇z
      //  SetDebugVal(3,  ff[3]);  // ff[3] = ω̇FL
       // SetDebugVal(4,  ff[4]);  // ff[4] = ω̇FR
      //  SetDebugVal(5,  ff[5]);  // ff[0] = ω̇x
       // SetDebugVal(6,  ff[6]);  // ff[1] = ω̇y
       // SetDebugVal(7,  ff[7]);  // ff[2] = ω̇z
        //SetDebugVal(8,  ff[8]);  // ff[3] = ω̇FL
        SetDebugVal(9,  ff[9]);  // ff[4] = ω̇FR
        SetDebugVal(10,  ff[10]);  // ff[0] = ω̇x
        SetDebugVal(11,  ff[11]);  // ff[1] = ω̇y
        SetDebugVal(12,  ff[12]);  // ff[2] = ω̇z
        SetDebugVal(13,  ff[13]);  // ff[3] = ω̇FL
        SetDebugVal(14,  ff[14]);  // ff[4] = ω̇FR
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