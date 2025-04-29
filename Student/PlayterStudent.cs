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
         sys = new LinSysEq(8);
         Amat = new double[8,8];  // an 8 by 8 array of doubles
         Bmat = new double[8];    // an 8 by 1 array of doubles

           // Set phi and precompute trigonometric values once:
         phi = Math.PI / 4;   // example: ϕ = 45 degrees, adjust as necessary
        sinPhi = Math.Sin(phi);
        cosPhi = Math.Cos(phi);
 
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
         
         ff[0] = -AngVelCrossAngMo.x/rho2;
         ff[1] = -AngVelCrossAngMo.y/(rho2*gammaY);
         ff[2] = -AngVelCrossAngMo.z/(rho2*gammaZ);

         
        Vex vG_vx = new Vex(
            q0*q0 + q1*q1 - q2*q2 - q3*q3,
           -2*q0*q3   + 2*q1*q2,
            2*q0*q2   + 2*q1*q3
        );
        Vex vG_vy = new Vex(
            2*q0*q3   + 2*q1*q2,
            q0*q0 - q1*q1 + q2*q2 - q3*q3,
           -2*q0*q1   + 2*q2*q3
        );
        Vex vG_vz = new Vex(
           -2*q0*q2   + 2*q1*q3,
            2*q0*q1   + 2*q2*q3,
            q0*q0 - q1*q1 - q2*q2 + q3*q3
        );
        Vex vG_B = vx * vG_vx + vy * vG_vy + vz * vG_vz;

        Vex bX = new Vex(1,0,0);
        Vex bY = new Vex(0,1,0);
        Vex bZ = new Vex(0,0,1);
        // hinge axis in the shoulder frame
        Vex sZ = new Vex(0.0, -sinPhi, cosPhi);

        // 2) Position vectors (you may already have these)
        Vex rSL_G  = new Vex( 1.0, h, 0.0 );                           // left shoulder → CG
        Vex rFL_SL = L * new Vex(Math.Cos(thetaL), Math.Sin(thetaL)*cosPhi, Math.Sin(thetaL)*sinPhi);
        Vex rFL_G  = rSL_G + rFL_SL;                                   // CG → left arm mass

        Vex rSR_G  = new Vex(-1.0, h, 0.0 );                           // right shoulder → CG
        Vex rFR_SR = -L * new Vex(Math.Cos(thetaR), Math.Sin(thetaR)*cosPhi, Math.Sin(thetaR)*sinPhi);
        Vex rFR_G  = rSR_G + rFR_SR;                                   // CG → right arm mass

        // 3) Left‐arm partial velocities (Eq.27)
        //    ωx‐partial
        Vex vFL_wx = Vex.Cross(bX, rFL_G);
        //    ωy‐partial
        Vex vFL_wy = Vex.Cross(bY, rFL_G);
        //    ωz‐partial
        Vex vFL_wz = Vex.Cross(bZ, rFL_G);
        //    ωL‐partial (hinge speed)
        Vex vFL_wL = Vex.Cross(sZ, rFL_SL);

        // 4) Right‐arm partial velocities
        Vex vFR_wx = Vex.Cross(bX, rFR_G);
        Vex vFR_wy = Vex.Cross(bY, rFR_G);
        Vex vFR_wz = Vex.Cross(bZ, rFR_G);
        Vex vFR_wR = Vex.Cross(sZ, rFR_SR);
        // 1) First‐term: translational partials
        Vex aFL = vx * vG_vx
                + vy * vG_vy
                + vz * vG_vz;

        // 2) Rotational partials from body‐spin
        aFL += omegaX * vFL_wx
            + omegaY * vFL_wy
             + omegaZ * vFL_wz;

        // 3) Rotational partial from hinge speed
        aFL += omegaFL * vFL_wL;

        // 4) Transport/coriolis terms:
        //    NωB × (NωB × rSL/G)
        Vex omegaB   = new Vex(omegaX, omegaY, omegaZ);
        Vex term4a   = Vex.Cross(omegaB, Vex.Cross(omegaB, rSL_G));

        //    (NωB × BωFL) × rFL/SL
        Vex omegaB_BwFL = Vex.Cross(omegaB, sZ * omegaFL);
        Vex term4b      = Vex.Cross(omegaB_BwFL, rFL_SL);

        //    NωFL × (NωFL × rFL/SL)
        Vex omegaFL_vec = sZ * omegaFL; 
        Vex term4c      = Vex.Cross(omegaFL_vec, Vex.Cross(omegaFL_vec, rFL_SL));

        // sum transport terms
        aFL += term4a + term4b + term4c;

        // 5) Dump into debug slots so you can watch the vector
        SetDebugVal(5, aFL.x);
        SetDebugVal(6, aFL.y);
        SetDebugVal(7, aFL.z);

// 6) If you want to actually set ff for testing, you can do:
// ff[?] = aFL.x;  // but typically ff[3]–ff[7] reserved for speeds & CG
 
         SetDebugVal(0,omegaX); // use these for debugging,  displays on screen.
         SetDebugVal(1,omegaY);
         SetDebugVal(2,omegaZ); 
         SetDebugVal(3, ff[3]);
         SetDebugVal(4, ff[4]);
 
         ff[8] = .5*(-q1*omegaX - q2*omegaY - q3*omegaZ);
         ff[9] = .5*(q0*omegaX - q3*omegaY + q2*omegaZ);
         ff[10] = .5*(q3*omegaX + q0*omegaY - q1*omegaZ);
         ff[11] = .5*(-q2*omegaX + q1*omegaY + q0*omegaZ);

       

        ff[5] = vG_B.x;   // ẋG = vx
        ff[6] = vG_B.y;   // ẏG = vy
        ff[7] = vG_B.z;   // żG = vz


        ff[12] = omegaFL;   // θ̇L = ωFL
        ff[13] = omegaFR;   // θ̇R = ωFR

        ff[14] = vx;   // ẋG = vx
        ff[15] = vy;   // ẏG = vy
        ff[16] = vz;   // żG = vz




     //SetDebugVal(5,  vG_vx.x);   // TestVal_5
    // SetDebugVal(6,  vG_vx.y);   // TestVal_6
    // SetDebugVal(7,  vG_vx.z);   // TestVal_7
    SetDebugVal(8,  vG_vy.x);   // TestVal_8
    SetDebugVal(9,  vG_vy.y);   // TestVal_9
    SetDebugVal(10, vG_vy.z);   // TestVal_10
    SetDebugVal(11, vG_vz.x);   // TestVal_11
    SetDebugVal(12, vG_vz.y);   // TestVal_12
    SetDebugVal(13, vG_vz.z);   // TestVal_13
    SetDebugVal(12, vG_B.x);
    SetDebugVal(13, vG_B.y);
    SetDebugVal(14, vG_B.z);
    SetDebugVal(15, ff[15]);
    SetDebugVal(16, ff[16]);

 
         // COMMENT THESE OUT OR REMOVE WHEN READY
         //ff[0] = ff[1] = ff[2] = 0.0;   // derivs of body angular velocities set to zero
         //ff[3] = ff[4] = 0.0;           // derivs of arm angular velocities
         //ff[5] = ff[6] = ff[7] = 0.0;   // derivs of cener of mass velocities
         //ff[8] = ff[9] = ff[10] = ff[11] = 0.0;  // derivs of quaternion coords
         //ff[12] = ff[13] = 0.0;         // derivs of arm angles
         //ff[14] = ff[15] = ff[16] = 0.0;  // derivs of CG coordinates
     }
 
 } // end class