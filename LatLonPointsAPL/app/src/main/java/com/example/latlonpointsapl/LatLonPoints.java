package com.example.latlonpointsapl;

public class LatLonPoints extends myPoint2D {
    final static double a=6378137;
    final static double rf=298.257222101;
    final static double m0=0.9999;
    final static double s2r=Math.PI/648000;
    final static double n=0.5/(rf-0.5);
    final static double n15=1.5*n;
    final static double anh=0.5*a/(1+n);
    final static double nsq=n*n;
    final static double e2n=2*Math.sqrt(n)/(1+n);
    final static double ra=2*anh*m0*(1+nsq/4+nsq*nsq/64);
    final static int jt=5;
    final static int jt2=2*jt;

    // 平面直角座標の座標系原点の緯度を度単位で、経度を分単位で格納
    final static double phi0[]={0,33,33,36,33,36,36,36,36,36,40,44,44,44,26,26,26,26,20,26};
    final static double lmbd0={0,7770,7860,7930,8010,8060,8160,8230,8310,8390,8450,8415,8535,8655,8520,7
            650,7440,7860,8160,9240};

    private double ep=1.0;
    private double e[];
    private double alp[];

    public LatLonPoints(double lat, double lon) {
        super(lat, lon);
        e= new double[2(jt+1)];
        alp= new double[6];
        for(int k=1; k<=jt; k++) {
            e[k]=n15/k-n;
            ep*=e[k];
            e[k+jt]=n15/(k+jt)-n;
        }

        // 展開パラメータの事前入力
        alp[1]=(1/2+(-2/3+(5/16+(41/180-127/288*n)*n)*n)*n)*n;
        alp[2]=(13/48+(-3/5+(557/1440+281/630*n)*n)*n)*nsq;
        alp[3]=(61/240+(-103/140+15061/26880*n)*n)*n*nsq;
        alp[4]=(49561/161280-179/168*n)*nsq*nsq;
        alp[5]=34729/80640*n*nsq*nsq;
    }

    public myPoint2D changeXY() {
        s =[0.0];
        t =[];
    }
    // 該当緯度の 2 倍角の入力により赤道からの子午線弧長を求める関数
    private double Merid(double phi2) {
                dc=2.0*Math.cos(phi2) ; s[1]=Math.sin(phi2)
        for(i=1; i<=jt2; i++) { s[i+1]=dc*s[i]-s[i-1] ; t[i]=(1.0/i-4.0*i)*s[i] }
        sum=0.0 ; c1=ep ; j=jt
        while(j) {
            c2=phi2 ; c3=2.0 ; l=j ; m=0
            while(l) { c2+=(c3/=e[l--])*t[++m]+(c3*=e[2*j-l])*t[++m] }
            sum+=c1*c1*c2 ; c1/=e[j--]
        }
        return anh*(sum+phi2)
    }
// 与件入力
        num=eval(prompt("座標系番号を入力してください。"))
        phi=eval(prompt("緯度を ddmmss.ssss 形式で入力してください。"))
        lmbd=eval(prompt("経度を dddmmss.ssss 形式で入力してください。"))
        phideg=Math.floor(phi/10000) ; phimin=Math.floor((phi-phideg*10000)/100)
        phirad=(phideg*3600+phimin*60+phi-phideg*10000-phimin*100)*s2r
        lmbddeg=Math.floor(lmbd/10000) ; lmbdmin=Math.floor((lmbd-lmbddeg*10000)/100)
        lmbdsec=lmbddeg*3600+lmbdmin*60+lmbd-lmbddeg*10000-lmbdmin*100
// 実際の計算実行部分
        sphi=Math.sin(phirad) ; nphi=(1-n)/(1+n)*Math.tan(phirad)
        dlmbd=(lmbdsec-lmbd0[num]*60)*s2r
        sdlmbd=Math.sin(dlmbd) ; cdlmbd=Math.cos(dlmbd)
        tchi=sinh(arctanh(sphi)-e2n*arctanh(e2n*sphi)) ; cchi=Math.sqrt(1+tchi*tchi)
        xi=xip=Math.atan2(tchi, cdlmbd) ; eta=etap=arctanh(sdlmbd/cchi) ; sgm=1 ; tau=0
        for(j=alp.length; --j; ) {
            alsin=alp[j]*Math.sin(2*j*xip) ; alcos=alp[j]*Math.cos(2*j*xip)
            xi+=alsin*cosh(2*j*etap) ; eta+=alcos*sinh(2*j*etap)
            sgm+=2*j*alcos*cosh(2*j*etap) ; tau+=2*j*alsin*sinh(2*j*etap)
        }
        x=ra*xi-m0*Merid(2*phi0[num]*3600*s2r) ; y=ra*eta
        gmm=Math.atan2(tau*cchi*cdlmbd+sgm*tchi*sdlmbd, sgm*cchi*cdlmbd-tau*tchi*sdlmbd)
        m=ra/a*Math.sqrt((sgm*sgm+tau*tau)/(tchi*tchi+cdlmbd*cdlmbd)*(1+nphi*nphi))
// ラジアン → 度分秒変換

        sgn=(gmm<0)
        gdo=Math.floor(gmm/s2r/3600)+sgn
        gfun=Math.floor((gmm/s2r-gdo*3600)/60)+sgn
        gbyou=gmm/s2r-gdo*3600-gfun*60
// 結果表示
        document.write("<h2>座標系番号： " + num + " 緯度： " + phi + " 経度： " + lmbd +
                "<br/><br/>")
        document.write("Ｘ＝" + x + "，Ｙ＝" + y + "<br/>")
        document.write("γ＝" + (sgn?"－":"＋") + Math.abs(gdo) + "°" + Math.abs(gfun) + "′
                " + Math.abs(gbyou) + "″，m＝" + m + "<br/></h2>")
// --></script>
    }

    public double sinh(x) { return 0.5*(Math.exp(x)-Math.exp(-x)); }
    public double cosh(x) { return 0.5*(Math.exp(x)+Math.exp(-x)); }
    public double arctanh(x) { return 0.5*Math.log((1+x)/(1-x)); }

}
