// This file is derived from TMinuit.h in ROOT version 5.30/06.
// It is distributed under the GNU LGPL. The original copyright notice
// reads:

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/
// ---------------------------------- minuit.h

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TMinuit                                                              //
//                                                                      //
// The MINUIT minimisation package (base class)                         //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#ifndef TMINUIT_HH
#define TMINUIT_HH

#include <string>
using namespace std;

class TMinuit {
public:
    static TMinuit*  gMinuit;
    enum {kMAXWARN=100};

    int        fNpfix;            //Number of fixed parameters
    int        fEmpty;            //Initialization flag (1 = Minuit initialized)
    int        fMaxpar;           //Maximum number of parameters
    int        fMaxint;           //Maximum number of internal parameters
    int        fNpar;             //Number of free parameters (total number of pars = fNpar + fNfix)
    int        fMaxext;           //Maximum number of external parameters
    int        fMaxIterations;    //Maximum number of iterations
    int        fMaxpar5;          // fMaxpar*(fMaxpar+1)/2
    int        fMaxcpt;
    int        fMaxpar2;          // fMaxpar*fMaxpar
    int        fMaxpar1;          // fMaxpar*(fMaxpar+1)

    double     fAmin;             //Minimum value found for FCN
    double     fUp;               //FCN+-UP defines errors (for chisquare fits UP=1)
    double     fEDM;              //Estimated vertical distance to the minimum
    double     fFval3;            //
    double     fEpsi;             //
    double     fApsi;             //
    double     fDcovar;           //Relative change in covariance matrix
    double     fEpsmac;           //machine precision for floating points:
    double     fEpsma2;           //sqrt(fEpsmac)
    double     fVlimlo;           //
    double     fVlimhi;           //
    double     fUndefi;           //Undefined number = -54321
    double     fBigedm;           //Big EDM = 123456
    double     fUpdflt;           //
    double     fXmidcr;           //
    double     fYmidcr;           //
    double     fXdircr;           //
    double     fYdircr;           //

    double*     fU;               //[fMaxpar2] External (visible to user in FCN) value of parameters
    double*     fAlim;            //[fMaxpar2] Lower limits for parameters. If zero no limits
    double*     fBlim;            //[fMaxpar2] Upper limits for parameters
    double*     fErp;             //[fMaxpar] Positive Minos errors if calculated
    double*     fErn;             //[fMaxpar] Negative Minos errors if calculated
    double*     fWerr;            //[fMaxpar] External parameters error (standard deviation, defined by UP)
    double*     fGlobcc;          //[fMaxpar] Global Correlation Coefficients
    double*     fX;               //[fMaxpar] Internal parameters values
    double*     fXt;              //[fMaxpar] Internal parameters values X saved as Xt
    double*     fDirin;           //[fMaxpar] (Internal) step sizes for current step
    double*     fXs;              //[fMaxpar] Internal parameters values saved for fixed params
    double*     fXts;             //[fMaxpar] Internal parameters values X saved as Xt for fixed params
    double*     fDirins;          //[fMaxpar] (Internal) step sizes for current step for fixed params
    double*     fGrd;             //[fMaxpar] First derivatives
    double*     fG2;              //[fMaxpar]
    double*     fGstep;           //[fMaxpar] Step sizes
    double*     fGin;             //[fMaxpar2]
    double*     fDgrd;            //[fMaxpar] Uncertainties
    double*     fGrds;            //[fMaxpar]
    double*     fG2s;             //[fMaxpar]
    double*     fGsteps;          //[fMaxpar]
    double*     fVhmat;           //[fMaxpar5] (Internal) error matrix stored as Half MATrix, since it is symmetric
    double*     fVthmat;          //[fMaxpar5] VHMAT is sometimes saved in VTHMAT, especially in MNMNOT
    double*     fP;               //[fMaxpar1]
    double*     fPstar;           //[fMaxpar2]
    double*     fPstst;           //[fMaxpar]
    double*     fPbar;            //[fMaxpar]
    double*     fPrho;            //[fMaxpar] Minimum point of parabola
    double*     fWord7;           //[fMaxpar]
    double*     fXpt;             //[fMaxcpt] X array of points for contours
    double*     fYpt;             //[fMaxcpt] Y array of points for contours

    double*     fCONTgcc;         //[fMaxpar] array used in mncont
    double*     fCONTw;           //[fMaxpar] array used in mncont
    double*     fFIXPyy;          //[fMaxpar] array used in mnfixp
    double*     fGRADgf;          //[fMaxpar] array used in mngrad
    double*     fHESSyy;          //[fMaxpar] array used in mnhess
    double*     fIMPRdsav;        //[fMaxpar] array used in mnimpr
    double*     fIMPRy;           //[fMaxpar] array used in mnimpr
    double*     fMATUvline;       //[fMaxpar] array used in mnmatu
    double*     fMIGRflnu;        //[fMaxpar] array used in mnmigr
    double*     fMIGRstep;        //[fMaxpar] array used in mnmigr
    double*     fMIGRgs;          //[fMaxpar] array used in mnmigr
    double*     fMIGRvg;          //[fMaxpar] array used in mnmigr
    double*     fMIGRxxs;         //[fMaxpar] array used in mnmigr
    double*     fMNOTxdev;        //[fMaxpar] array used in mnmnot
    double*     fMNOTw;           //[fMaxpar] array used in mnmnot
    double*     fMNOTgcc;         //[fMaxpar] array used in mnmnot
    double*     fPSDFs;           //[fMaxpar] array used in mnpsdf
    double*     fSEEKxmid;        //[fMaxpar] array used in mnseek
    double*     fSEEKxbest;       //[fMaxpar] array used in mnseek
    double*     fSIMPy;           //[fMaxpar] array used in mnsimp
    double*     fVERTq;           //[fMaxpar] array used in mnvert
    double*     fVERTs;           //[fMaxpar] array used in mnvert
    double*     fVERTpp;          //[fMaxpar] array used in mnvert
    double*     fCOMDplist;       //[fMaxpar] array used in mncomd
    double*     fPARSplist;       //[fMaxpar] array used in mnpars

    int*        fNvarl;           //[fMaxpar2] parameters flag (-1=undefined, 0=constant..)
    int*        fNiofex;          //[fMaxpar2] Internal parameters number, or zero if not currently variable
    int*        fNexofi;          //[fMaxpar] External parameters number for currently variable parameters
    int*        fIpfix;           //[fMaxpar] List of fixed parameters
    int        fNu;               //
    int        fIsysrd;           //standardInput unit
    int        fIsyswr;           //standard output unit
    int        fIsyssa;           //
    int        fNpagwd;           //Page width
    int        fNpagln;           //Number of lines per page
    int        fNewpag;           //
    int        fIstkrd[10];       //
    int        fNstkrd;           //
    int        fIstkwr[10];       //
    int        fNstkwr;           //
    int        fISW[7];           //Array of switches
    int        fIdbg[11];         //Array of internal debug switches
    int        fNblock;           //Number of Minuit data blocks
    int        fIcomnd;           //Number of commands
    int        fNfcn;             //Number of calls to FCN
    int        fNfcnmx;           //Maximum number of calls to FCN
    int        fNfcnlc;           //
    int        fNfcnfr;           //
    int        fItaur;            //
    int        fIstrat;           //
    int        fNwrmes[2];        //
    int        fNfcwar[20];       //
    int        fIcirc[2];         //
    int        fStatus;           //Status flag for the last called Minuit function
    int        fKe1cr;            //
    int        fKe2cr;            //
    bool       fLwarn;            //true if warning messges are to be put out (default=true)
    bool       fLrepor;           //true if exceptional conditions are put out (default=false)
    bool       fLimset;           //true if a parameter is up against limits (for MINOS)
    bool       fLnolim;           //true if there are no limits on any parameters (not yet used)
    bool       fLnewmn;           //true if the previous process has unexpectedly improved FCN
    bool       fLphead;           //true if a heading should be put out for the next parameter definition
    //bool       fGraphicsMode;     //true if graphics mode on (default)
    char*         fChpt;            //!Character to be plotted at the X,Y contour positions
    string*      fCpnam;           //[fMaxpar2] Array of parameters names
    string      fCfrom;            //
    string      fCstatu;           //
    string      fCtitl;            //
    string      fCword;            //
    string      fCundef;           //
    string      fCvrsn;            //
    string      fCovmes[4];        //
    string      fOrigin[kMAXWARN]; //
    string      fWarmes[kMAXWARN]; //
    //TObject      *fObjectFit;       //Pointer to object being fitted
    //TObject      *fPlot;            //Pointer to TGraph object created by mncont
    //TMethodCall  *fMethodCall;      //Pointer to MethodCall in case of interpreted function
    void (*fFCN)(int& npar, double* gin, double& f, double* u, int flag);         //!

    // methods performed on TMinuit class
public:

    TMinuit(int maxpar);
    virtual       ~TMinuit();
    virtual void   BuildArrays(int maxpar=15);
    virtual int  Command(const char* command);
    virtual int  DefineParameter(int parNo, const char* name, double initVal, double initErr, double lowerLimit,
                                 double upperLimit);
    virtual void   DeleteArrays();
    virtual int  Eval(int npar, double* grad, double& fval, double* par, int flag);
    virtual int  FixParameter(int parNo);
    //TMethodCall   *GetMethodCall() const {return fMethodCall;}
    //TObject       *GetObjectFit() const {return fObjectFit;}
    int          GetMaxIterations() const {
        return fMaxIterations;
    }
    virtual int  GetNumFixedPars() const;
    virtual int  GetNumFreePars() const;
    virtual int  GetNumPars() const;
    virtual int  GetParameter(int parNo, double& currentValue, double& currentError) const;
    //virtual TObject *GetPlot() const {return fPlot;}
    int          GetStatus() const {
        return fStatus;
    }
    virtual int  Migrad();
    virtual void   mnamin();
    virtual void   mnbins(double a1, double a2, int naa, double& bl, double& bh, int& nb, double& bwid);
    virtual void   mncalf(double* pvec, double& ycalf);
    virtual void   mncler();
    virtual void   mncntr(int ke1, int ke2, int& ierrf);
    virtual void   mncomd(const char* crdbin, int& icondn);
    virtual void   mncont(int ke1, int ke2, int nptu, double* xptu, double* yptu, int& ierrf);
    virtual void   mncrck(string crdbuf, int maxcwd, string& comand, int& lnc
                          ,  int mxp, double* plist, int& llist, int& ierr, int isyswr);
    virtual void   mncros(double& aopt, int& iercr);
    virtual void   mncuve();
    virtual void   mnderi();
    virtual void   mndxdi(double pint, int ipar, double& dxdi);
    virtual void   mneig(double* a, int ndima, int n, int mits, double* work, double precis, int& ifault);
    virtual void   mnemat(double* emat, int ndim);
    virtual void   mnerrs(int number, double& eplus, double& eminus, double& eparab, double& gcc);
    virtual void   mneval(double anext, double& fnext, int& ierev);
    virtual void   mnexcm(const char* comand, double* plist, int llist, int& ierflg) ;
    virtual void   mnexin(double* pint);
    virtual void   mnfixp(int iint, int& ierr);
    virtual void   mnfree(int k);
    virtual void   mngrad();
    virtual void   mnhelp(string comd);
    virtual void   mnhelp(const char* command="");
    virtual void   mnhess();
    virtual void   mnhes1();
    virtual void   mnimpr();
    virtual void   mninex(double* pint);
    virtual void   mninit(int i1, int i2, int i3);
    virtual void   mnlims();
    virtual void   mnline(double* start, double fstart, double* step, double slope, double toler);
    virtual void   mnmatu(int kode);
    virtual void   mnmigr();
    virtual void   mnmnos();
    virtual void   mnmnot(int ilax, int ilax2, double& val2pl, double& val2mi);
    virtual void   mnparm(int k, string cnamj, double uk, double wk, double a, double b, int& ierflg);
    virtual void   mnpars(string& crdbuf, int& icondn);
    virtual void   mnpfit(double* parx2p, double* pary2p, int npar2p, double* coef2p, double& sdev2p);
    virtual void   mnpint(double& pexti, int i, double& pinti);
    virtual void   mnplot(double* xpt, double* ypt, char* chpt, int nxypt, int npagwd, int npagln);
    virtual void   mnpout(int iuext, string& chnam, double& val, double& err, double& xlolim, double& xuplim,
                          int& iuint) const;
    virtual void   mnprin(int inkode, double fval);
    virtual void   mnpsdf();
    virtual void   mnrazz(double ynew, double* pnew, double* y, int& jh, int& jl);
    virtual void   mnrn15(double& val, int& inseed);
    virtual void   mnrset(int iopt);
    virtual void   mnsave();
    virtual void   mnscan();
    virtual void   mnseek();
    virtual void   mnset();
    virtual void   mnsimp();
    virtual void   mnstat(double& fmin, double& fedm, double& errdef, int& npari, int& nparx, int& istat);
    virtual void   mntiny(double epsp1, double& epsbak);
    bool         mnunpt(string& cfname);
    virtual void   mnvert(double* a, int l, int m, int n, int& ifail);
    virtual void   mnwarn(const char* copt, const char* corg, const char* cmes);
    virtual void   mnwerr();
    virtual int  Release(int parNo);
    virtual int  SetErrorDef(double up);
    virtual void   SetFCN(void (*fcn)(int&, double*, double& f, double*, int));

    //virtual void   SetGraphicsMode(bool mode=kTRUE) {fGraphicsMode = mode;}
    virtual void   SetMaxIterations(int maxiter=500) {
        fMaxIterations = maxiter;
    }
    //virtual void   SetObjectFit(TObject *obj) {fObjectFit=obj;}
    virtual int  SetPrintLevel(int printLevel=0);

#ifdef NEEDS_CAREFUL_CONVERSION
    virtual TObject* Contour(int npoints=10, int pa1=0, int pa2=1);
#endif
};

#endif

