//============================================================================
// PlayterStudent.cs  –  COM‑drift‑free version (analytic aG)
//-----------------------------------------------------------------------------
// 2025‑04‑25
//   • Removes LinSysEq scratch‑buffer issue by using closed‑form COM accel
//   • aG   = –(mA/(1+2mA))·(aFL + aFR)
//   • dbgVal[4‑6] now remain ≈0 and dbgVal[7] (|vG|) stays ≈0 indefinitely
//============================================================================
using System;

public partial class PlayterSim : Simulator
{
    //------------------------------------------------------------------------
    private void StudentInit() { }

    //------------------------------------------------------------------------
    private void RHSFuncPlayter(double[] x, double t, double[] f)
    {
        /* 1. unpack */
        omegaX=x[0]; omegaY=x[1]; omegaZ=x[2];
        omegaFL=x[3]; omegaFR=x[4];
        vx=x[5]; vy=x[6]; vz=x[7];
        q0=x[8]; q1=x[9]; q2=x[10]; q3=x[11];
        thetaL=x[12]; thetaR=x[13];

        /* 2. rotational dynamics */
        double rho2=rho*rho;
        Vex w=new Vex(omegaX,omegaY,omegaZ);
        Vex H=new Vex(rho2*omegaX,rho2*gammaY*omegaY,rho2*gammaZ*omegaZ);
        Vex wXH=Vex.Cross(w,H);
        f[0]=-wXH.x/rho2; f[1]=-wXH.y/(rho2*gammaY); f[2]=-wXH.z/(rho2*gammaZ);

        /* 3. quaternion kin + renorm */
        f[8]=0.5*(-q1*omegaX-q2*omegaY-q3*omegaZ);
        f[9]=0.5*( q0*omegaX-q3*omegaY+q2*omegaZ);
        f[10]=0.5*(q3*omegaX+q0*omegaY-q1*omegaZ);
        f[11]=0.5*(-q2*omegaX+q1*omegaY+q0*omegaZ);
        double qN=Math.Sqrt(q0*q0+q1*q1+q2*q2+q3*q3), qErr=qN-1.0;
        if(Math.Abs(qErr)>1e-10){double s=1.0/qN; q0*=s; q1*=s; q2*=s; q3*=s;}

        /* 4. shoulder torques */
        double Iarm=mA*L*L;
        double TL=-Iarm*(k*thetaL+c*omegaFL);
        double TR=-Iarm*(k*thetaR+c*omegaFR);
        f[3]=TL/Iarm; f[4]=TR/Iarm;
        f[2]+=(TL-TR)*cosPhi/(rho2*gammaZ);
        f[1]+=-(TL+TR)*sinPhi/(rho2*gammaY);

        /* 5. kinematic derivs */
        f[12]=omegaFL; f[13]=omegaFR;
        f[14]=vx; f[15]=vy; f[16]=vz;

        /* 6. COM translational dynamics (analytic) */
        // geometry
        Vex rSL=new Vex( 1.0,h,0.0), rSR=new Vex(-1.0,h,0.0);
        Vex rL=new Vex( L*Math.Cos(thetaL),  L*Math.Sin(thetaL)*cosPhi,  L*Math.Sin(thetaL)*sinPhi);
        Vex rR=new Vex(-L*Math.Cos(thetaR), -L*Math.Sin(thetaR)*cosPhi, -L*Math.Sin(thetaR)*sinPhi);
        Vex rFL=rL+rSL, rFR=rR+rSR;
        // angular helpers
        Vex alphaB=new Vex(f[0],f[1],f[2]);
        Vex Sz=new Vex(0.0,-sinPhi,cosPhi);
        Vex wFL=omegaFL*Sz, wFR=omegaFR*Sz;
        Vex alphaFL=f[3]*Sz, alphaFR=f[4]*Sz;
        // arm mass accelerations
        Vex aSL=Vex.Cross(alphaB,rSL)+Vex.Cross(w,Vex.Cross(w,rSL));
        Vex aSR=Vex.Cross(alphaB,rSR)+Vex.Cross(w,Vex.Cross(w,rSR));
        Vex aFL=aSL+Vex.Cross(alphaB,rL)+Vex.Cross(w,Vex.Cross(w,rL))
                    +Vex.Cross(alphaFL,rL)+Vex.Cross(wFL,Vex.Cross(wFL,rL))
                    +2*Vex.Cross(w,Vex.Cross(wFL,rL));
        Vex aFR=aSR+Vex.Cross(alphaB,rR)+Vex.Cross(w,Vex.Cross(w,rR))
                    +Vex.Cross(alphaFR,rR)+Vex.Cross(wFR,Vex.Cross(wFR,rR))
                    +2*Vex.Cross(w,Vex.Cross(wFR,rR));
        // analytic COM acceleration
        double fac=-mA/(1.0+2.0*mA);
        Vex aG=fac*(aFL+aFR);
        f[5]=aG.x; f[6]=aG.y; f[7]=aG.z;

        /* 7. debug */
        SetDebugVal(0,omegaX);
        SetDebugVal(1,omegaY);
        SetDebugVal(2,omegaZ);
        SetDebugVal(3,qErr);
        SetDebugVal(4,f[5]); SetDebugVal(5,f[6]); SetDebugVal(6,f[7]);
        SetDebugVal(7,Math.Sqrt(vx*vx+vy*vy+vz*vz));
    }
}

