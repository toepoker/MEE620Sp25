//============================================================================
// PlayterStudent.cs – Playter doll full Kane matrix (no tuples)          
//-----------------------------------------------------------------------------
// 2025‑04‑26                                                             
// * Rotational rows include arm inertia + gyroscopic coupling columns.   
// * Translational rows use analytic aG_const (no COM drift).             
//============================================================================
using System;

public partial class PlayterSim : Simulator
{
    /* one‑time allocations */
    LinSysEq sys8;
    double[,] A;
    double[]  B;

    private void StudentInit()
    {
        sys8 = new LinSysEq(8);
        A = new double[8,8];
        B = new double[8];
    }

    /* ------------------------------------------------------------------ */
    /* RHS  ˣ̇ = f(x,t)                                                     */
    /* ------------------------------------------------------------------ */
    private void RHSFuncPlayter(double[] x, double t, double[] f)
    {
        /* unpack state */
        omegaX = x[0]; omegaY = x[1]; omegaZ = x[2];
        omegaFL= x[3]; omegaFR= x[4];
        vx     = x[5]; vy     = x[6]; vz     = x[7];
        q0=x[8]; q1=x[9]; q2=x[10]; q3=x[11];
        thetaL=x[12]; thetaR=x[13];

        double rho2 = rho*rho;
        double Iarm = mA*L*L;

        /* ---------------- Geometry ------------------------------------ */
        Vex rSL = new Vex( 1.0, h, 0.0);
        Vex rSR = new Vex(-1.0, h, 0.0);
        Vex rFLS= new Vex( L*Math.Cos(thetaL),  L*Math.Sin(thetaL)*cosPhi,  L*Math.Sin(thetaL)*sinPhi);
        Vex rFRS= new Vex(-L*Math.Cos(thetaR), -L*Math.Sin(thetaR)*cosPhi, -L*Math.Sin(thetaR)*sinPhi);
        Vex rFL = rSL + rFLS;
        Vex rFR = rSR + rFRS;

        /* ---------------- Arm inertia components (no tuples) ---------- */
        double rsqL = Vex.Dot(rFL, rFL);
        double ILxx = mA*(rsqL - rFL.x*rFL.x);
        double ILyy = mA*(rsqL - rFL.y*rFL.y);
        double ILzz = mA*(rsqL - rFL.z*rFL.z);
        double ILxy = -mA*rFL.x*rFL.y;

        double rsqR = Vex.Dot(rFR, rFR);
        double IRxx = mA*(rsqR - rFR.x*rFR.x);
        double IRyy = mA*(rsqR - rFR.y*rFR.y);
        double IRzz = mA*(rsqR - rFR.z*rFR.z);
        double IRxy = -mA*rFR.x*rFR.y;

        double Ix  = rho2        + ILxx + IRxx;
        double Iy  = rho2*gammaY + ILyy + IRyy;
        double Iz  = rho2*gammaZ + ILzz + IRzz;
        double Ixy = ILxy + IRxy;

        /* Gyroscopic coupling columns */
        Vex Sz = new Vex(0.0, -sinPhi, cosPhi);
        Vex coupL = mA * Vex.Cross(rFL, Vex.Cross(Sz, rFLS));
        Vex coupR = mA * Vex.Cross(rFR, Vex.Cross(Sz, rFRS));

        /* ---------------- Arm torques ---------------- */
  
        double Iarm = mA * L * L;                 // nondimensional Iarm
        double TL   = -Iarm * (k*thetaL + c*omegaFL);
        double TR   = -Iarm * (k*thetaR + c*omegaFR);


        /* ---------------- Clear A, B ----------------- */
        Array.Clear(A,0,A.Length); Array.Clear(B,0,B.Length);

        /* -------- Rotational rows 0‑2 --------------- */
        // Row ω̇x
        A[0,0]=Ix;     A[0,1]=-Ixy;  A[0,3]=coupL.x; A[0,4]=coupR.x;
        // Row ω̇y
        A[1,0]=-Ixy;   A[1,1]=Iy;    A[1,3]=coupL.y; A[1,4]=coupR.y;
        // Row ω̇z
        A[2,2]=Iz;                      A[2,3]=coupL.z; A[2,4]=coupR.z;

        Vex w = new Vex(omegaX,omegaY,omegaZ);
        Vex Hbody = new Vex(rho2*omegaX, rho2*gammaY*omegaY, rho2*gammaZ*omegaZ);
        Vex wXH = Vex.Cross(w,Hbody);
        B[0] = -wXH.x;
        B[1] = -wXH.y;
        B[2] = -wXH.z;

        /* Reaction couples from shoulder torques */
        B[1] += -(TL+TR)*sinPhi;
        B[2] +=  (TL-TR)*cosPhi;

        /* -------- Hinge rows 3-4 --------------- */
        A[3,3]=Iarm;  B[3]=TL;
        A[4,4]=Iarm;  B[4]=TR;

        /* -------- Translational rows 5‑7 (analytic aG) -------- */
        // constant part of aG
        Vex wFL = omegaFL*Sz;
        Vex wFR = omegaFR*Sz;
        Vex aGconst = (-mA/(1+2*mA))*(
            Vex.Cross(w, Vex.Cross(w, rFL)) + Vex.Cross(wFL, Vex.Cross(wFL,rFLS)) +
            Vex.Cross(w, Vex.Cross(w, rFR)) + Vex.Cross(wFR, Vex.Cross(wFR,rFRS)) );

        // partial‑velocity columns (Eq. 18 from notes)
        Vex vpx = new Vex(q0*q0+q1*q1-q2*q2-q3*q3, -2*q0*q3+2*q1*q2,  2*q0*q2+2*q1*q3);
        Vex vpy = new Vex( 2*q0*q3+2*q1*q2,        q0*q0-q1*q1+q2*q2-q3*q3, -2*q0*q1+2*q2*q3);
        Vex vpz = new Vex(-2*q0*q2+2*q1*q3,        2*q0*q1+2*q2*q3,         q0*q0-q1*q1-q2*q2+q3*q3);

        A[5,5]=vpx.x; A[5,6]=vpy.x; A[5,7]=vpz.x;
        A[6,5]=vpx.y; A[6,6]=vpy.y; A[6,7]=vpz.y;
        A[7,5]=vpx.z; A[7,6]=vpy.z; A[7,7]=vpz.z;

        B[5]=aGconst.x; B[6]=aGconst.y; B[7]=aGconst.z;

        /* -------- Load into solver & solve -------- */
        for(int i=0;i<8;i++)
        {
            for(int j=0;j<8;j++) sys8.SetA(i,j,A[i,j]);
            sys8.SetB(i,B[i]);
        }
        sys8.SolveGauss();
        for(int i=0;i<8;i++) f[i]=sys8.Sol(i);

        /* -------- Quaternion kin -------- */
        f[8]=0.5*(-q1*omegaX - q2*omegaY - q3*omegaZ);
        f[9]=0.5*( q0*omegaX - q3*omegaY + q2*omegaZ);
        f[10]=0.5*(q3*omegaX + q0*omegaY - q1*omegaZ);
        f[11]=0.5*(-q2*omegaX + q1*omegaY + q0*omegaZ);

        /* -------- θ̇, ẋ ------------------ */
        f[12]=omegaFL;
        f[13]=omegaFR;
        f[14]=vx; f[15]=vy; f[16]=vz;

        /* Debug values -------------------- */
        SetDebugVal(0,omegaX);
        SetDebugVal(1,omegaY);
        SetDebugVal(2,omegaZ);
        SetDebugVal(3,Math.Sqrt(q0*q0+q1*q1+q2*q2+q3*q3)-1.0);
        SetDebugVal(4,f[5]); SetDebugVal(5,f[6]); SetDebugVal(6,f[7]);
        SetDebugVal(7,Math.Sqrt(vx*vx+vy*vy+vz*vz));
        SetDebugVal(8, omegaFL);   // TestVal_9 in UI (indices start at 0)
        SetDebugVal(9, omegaFR);   // TestVal_10

    }
}
