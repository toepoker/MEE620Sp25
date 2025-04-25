//============================================================================
// PlayterStudent.cs  –  compile‑clean version (no ZeroA call)
//-----------------------------------------------------------------------------
// • Replaces sys3.ZeroA() with a fresh LinSysEq(3) each frame.
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
        Vex w = new Vex(omegaX,omegaY,omegaZ);
        Vex H = new Vex(rho2*omegaX,rho2*gammaY*omegaY,rho2*gammaZ*omegaZ);
        Vex wXH=Vex.Cross(w,H);
        f[0]=-wXH.x/rho2; f[1]=-wXH.y/(rho2*gammaY); f[2]=-wXH.z/(rho2*gammaZ);

        /* 3. quaternion kin + renorm */
        f[8]=0.5*(-q1*omegaX-q2*omegaY-q3*omegaZ);
        f[9]=0.5*( q0*omegaX-q3*omegaY+q2*omegaZ);
        f[10]=0.5*(q3*omegaX+q0*omegaY-q1*omegaZ);
        f[11]=0.5*(-q2*omegaX+q1*omegaY+q0*omegaZ);
        double qN=Math.Sqrt(q0*q0+q1*q1+q2*q2+q3*q3), qErr=qN-1.0;
        if(Math.Abs(qErr)>1e-10){double inv=1.0/qN; q0*=inv; q1*=inv; q2*=inv; q3*=inv;}

        /* 4. shoulders */
        double Iarm=mA*L*L;
        double TL=-Iarm*(k*thetaL+c*omegaFL);
        double TR=-Iarm*(k*thetaR+c*omegaFR);
        f[3]=TL/Iarm; f[4]=TR/Iarm;
        f[2]+=(TL-TR)*cosPhi/(rho2*gammaZ);
        f[1]+=-(TL+TR)*sinPhi/(rho2*gammaY);

        /* 5. kinematic derivs */
        f[12]=omegaFL; f[13]=omegaFR;
        f[14]=vx; f[15]=vy; f[16]=vz;

        /* 6. COM translational balance */
        Vex rSL=new Vex( 1.0,h,0.0), rSR=new Vex(-1.0,h,0.0);
        Vex rL=new Vex( L*Math.Cos(thetaL), L*Math.Sin(thetaL)*cosPhi, L*Math.Sin(thetaL)*sinPhi);
        Vex rR=new Vex(-L*Math.Cos(thetaR),-L*Math.Sin(thetaR)*cosPhi,-L*Math.Sin(thetaR)*sinPhi);
        Vex rFL=rL+rSL, rFR=rR+rSR;
        Vex alphaB=new Vex(f[0],f[1],f[2]);
        Vex Sz=new Vex(0.0,-sinPhi,cosPhi);
        Vex wFL=omegaFL*Sz, wFR=omegaFR*Sz;
        Vex aSL=Vex.Cross(alphaB,rSL)+Vex.Cross(w,Vex.Cross(w,rSL));
        Vex aSR=Vex.Cross(alphaB,rSR)+Vex.Cross(w,Vex.Cross(w,rSR));
        Vex aFL=aSL+Vex.Cross(alphaB,rL)+Vex.Cross(w,Vex.Cross(w,rL))+Vex.Cross(f[3]*Sz,rL)+Vex.Cross(wFL,Vex.Cross(wFL,rL))+2*Vex.Cross(w,Vex.Cross(wFL,rL));
        Vex aFR=aSR+Vex.Cross(alphaB,rR)+Vex.Cross(w,Vex.Cross(w,rR))+Vex.Cross(f[4]*Sz,rR)+Vex.Cross(wFR,Vex.Cross(wFR,rR))+2*Vex.Cross(w,Vex.Cross(wFR,rR));
        double Mtot=1+2*mA; Vex rhs=-mA*(aFL+aFR);
        LinSysEq A=new LinSysEq(3);
        A.SetA(0,0,Mtot); A.SetA(1,1,Mtot); A.SetA(2,2,Mtot);
        A.SetB(0,rhs.x); A.SetB(1,rhs.y); A.SetB(2,rhs.z);
        A.SolveGauss();
        f[5]=A.Sol(0); f[6]=A.Sol(1); f[7]=A.Sol(2);

        /* 7. debug vals */
        SetDebugVal(0,omegaX); SetDebugVal(1,omegaY); SetDebugVal(2,omegaZ);
        SetDebugVal(3,qErr);
        SetDebugVal(4,f[5]); SetDebugVal(5,f[6]); SetDebugVal(6,f[7]);
        SetDebugVal(7,Math.Sqrt(vx*vx+vy*vy+vz*vz));
    }
}

