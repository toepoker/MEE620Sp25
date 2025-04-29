//============================================================================
// PlayterStudent.cs – Playter doll full Kane matrix (debug PVs added)
//-----------------------------------------------------------------------------
using System;

public partial class PlayterSim : Simulator
{
    /* one-time allocations */
    LinSysEq  sys8;
    double[,] A;
    double[]  B;

    private void StudentInit()
    {
        sys8 = new LinSysEq(8);
        A    = new double[8, 8];
        B    = new double[8];
    }

    /* ------------------------------------------------------------------ */
    /* RHS   ẋ = f(x , t)                                                */
    /* ------------------------------------------------------------------ */
    private void RHSFuncPlayter(double[] x, double t, double[] f)
    {
        /* ---------------- unpack state ---------------- */
        omegaX = x[0]; omegaY = x[1]; omegaZ = x[2];
        omegaFL = x[3]; omegaFR = x[4];
        vx = x[5]; vy = x[6]; vz = x[7];
        q0 = x[8]; q1 = x[9]; q2 = x[10]; q3 = x[11];
        thetaL = x[12]; thetaR = x[13];

        double rho2 = rho * rho;
        double Iarm = mA * L * L;                     // arm inertia (dimensionless)

        /* ---------------- geometry -------------------- */
        Vex rSL  = new Vex( 1.0, h, 0.0);
        Vex rSR  = new Vex(-1.0, h, 0.0);
        Vex rFLS = new Vex( L * Math.Cos(thetaL),
                            L * Math.Sin(thetaL) * cosPhi,
                            L * Math.Sin(thetaL) * sinPhi);
        Vex rFRS = new Vex(-L * Math.Cos(thetaR),
                           -L * Math.Sin(thetaR) * cosPhi,
                           -L * Math.Sin(thetaR) * sinPhi);
        Vex rFL  = rSL + rFLS;
        Vex rFR  = rSR + rFRS;

        /* ------------- body + arm inertia terms ------------- */
        double rsqL = Vex.Dot(rFL, rFL);
        double ILxx = mA * (rsqL - rFL.x * rFL.x);
        double ILyy = mA * (rsqL - rFL.y * rFL.y);
        double ILzz = mA * (rsqL - rFL.z * rFL.z);
        double ILxy = -mA * rFL.x * rFL.y;

        double rsqR = Vex.Dot(rFR, rFR);
        double IRxx = mA * (rsqR - rFR.x * rFR.x);
        double IRyy = mA * (rsqR - rFR.y * rFR.y);
        double IRzz = mA * (rsqR - rFR.z * rFR.z);
        double IRxy = -mA * rFR.x * rFR.y;

        double Ix  = rho2         + ILxx + IRxx;
        double Iy  = rho2 * gammaY + ILyy + IRyy;
        double Iz  = rho2 * gammaZ + ILzz + IRzz;
        double Ixy = ILxy + IRxy;

        /* ------------- gyroscopic coupling columns ------------- */
        Vex Sz    = new Vex(0.0, -sinPhi,  cosPhi);
        Vex coupL = mA * Vex.Cross(rFL, Vex.Cross(Sz, rFLS));
        Vex coupR = mA * Vex.Cross(rFR, Vex.Cross(Sz, rFRS));

        /* ------------- shoulder torques ------------------------ */
        double TL = -Iarm * (k * thetaL + c * omegaFL);
        double TR = -Iarm * (k * thetaR + c * omegaFR);

        /* ------------- clear A, B ------------------------------ */
        Array.Clear(A, 0, A.Length);
        Array.Clear(B, 0, B.Length);

        /* rotational rows (0-2) */
        A[0,0] = Ix;   A[0,1] = -Ixy; A[0,3] = coupL.x; A[0,4] = coupR.x;
        A[1,0] = -Ixy; A[1,1] = Iy;   A[1,3] = coupL.y; A[1,4] = coupR.y;
        A[2,2] = Iz;                       A[2,3] = coupL.z; A[2,4] = coupR.z;

        Vex  w   = new Vex(omegaX, omegaY, omegaZ);
        Vex  Hb  = new Vex(rho2 * omegaX, rho2 * gammaY * omegaY,
                           rho2 * gammaZ * omegaZ);
        Vex  wXH = Vex.Cross(w, Hb);
        B[0] = -wXH.x;
        B[1] = -wXH.y;
        B[2] = -wXH.z;
        B[1] += -(TL + TR) * sinPhi;   // reaction couples
        B[2] +=  (TL - TR) * cosPhi;

        /* hinge rows (3-4) */
        A[3,3] = Iarm;  B[3] = TL;
        A[4,4] = Iarm;  B[4] = TR;

        /* translational rows (5-7) – analytic COM acceleration */
        Vex wFL = omegaFL * Sz;
        Vex wFR = omegaFR * Sz;
        Vex aGconst = (-mA / (1 + 2 * mA)) * (
            Vex.Cross(w,   Vex.Cross(w,   rFL))
          + Vex.Cross(wFL, Vex.Cross(wFL, rFLS))
          + Vex.Cross(w,   Vex.Cross(w,   rFR))
          + Vex.Cross(wFR, Vex.Cross(wFR, rFRS)) );

        /* partial-velocity columns for ṽ_G */
        Vex vpx = new Vex(q0*q0 + q1*q1 - q2*q2 - q3*q3,
                          -2*q0*q3 + 2*q1*q2,
                           2*q0*q2 + 2*q1*q3);

        Vex vpy = new Vex( 2*q0*q3 + 2*q1*q2,
                           q0*q0 - q1*q1 + q2*q2 - q3*q3,
                          -2*q0*q1 + 2*q2*q3);

        Vex vpz = new Vex(-2*q0*q2 + 2*q1*q3,
                           2*q0*q1 + 2*q2*q3,
                           q0*q0 - q1*q1 - q2*q2 + q3*q3);

        A[5,5] = vpx.x; A[5,6] = vpy.x; A[5,7] = vpz.x;
        A[6,5] = vpx.y; A[6,6] = vpy.y; A[6,7] = vpz.y;
        A[7,5] = vpx.z; A[7,6] = vpy.z; A[7,7] = vpz.z;

        B[5] = aGconst.x;  B[6] = aGconst.y;  B[7] = aGconst.z;

        /* ------------ solve the 8×8 system -------------------- */
        for (int i = 0; i < 8; ++i)
        {
            for (int j = 0; j < 8; ++j) sys8.SetA(i, j, A[i, j]);
            sys8.SetB(i, B[i]);
        }
        sys8.SolveGauss();
        for (int i = 0; i < 8; ++i) f[i] = sys8.Sol(i);

        /* quaternion kinematics */
        f[8]  = 0.5 * (-q1*omegaX - q2*omegaY - q3*omegaZ);
        f[9]  = 0.5 * ( q0*omegaX - q3*omegaY + q2*omegaZ);
        f[10] = 0.5 * ( q3*omegaX + q0*omegaY - q1*omegaZ);
        f[11] = 0.5 * (-q2*omegaX + q1*omegaY + q0*omegaZ);

        /* θ̇ and COM velocity */
        f[12] = omegaFL;
        f[13] = omegaFR;
        f[14] = vx;
        f[15] = vy;
        f[16] = vz;

        /* ---------------- debug values ------------------------ */
        SetDebugVal(0,  omegaX);                               // TestVal_0
        SetDebugVal(1,  omegaY);                               // TestVal_1
        SetDebugVal(2,  omegaZ);                               // TestVal_2
        SetDebugVal(3,  Math.Sqrt(q0*q0+q1*q1+q2*q2+q3*q3)-1); // TestVal_3

        // partial velocities  (compare to professor’s screenshot)
        SetDebugVal(6,  vpx.x);  // TestVal_6
        SetDebugVal(7,  vpx.y);  // TestVal_7
        SetDebugVal(8,  vpx.z);  // TestVal_8

        SetDebugVal(9,  vpy.x);  // TestVal_9
        SetDebugVal(10, vpy.y);  // TestVal_10
        SetDebugVal(11, vpy.z);  // TestVal_11

        SetDebugVal(12, vpz.x);  // TestVal_12
        SetDebugVal(13, vpz.y);  // TestVal_13
        SetDebugVal(14, vpz.z);  // TestVal_14

        // COM-acceleration diagnostics (optional)
        SetDebugVal(15, f[5]);   // TestVal_15  aG.x
        SetDebugVal(16, f[6]);   // TestVal_16  aG.y
        SetDebugVal(17, f[7]);   // TestVal_17  aG.z
    }
}
