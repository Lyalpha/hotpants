typedef struct
{
   int       x0,y0;       /* origin of stamp in region coords*/
   int       x,y;         /* center of stamp in region coords*/
   int       nx,ny;       /* size of stamp */
   int       *xss;        /* x location of test substamp centers */
   int       *yss;        /* y location of test substamp centers */
   int       nss;         /* number of detected substamps, 1 .. nss     */
   int       sscnt;       /* represents which nss to use,  0 .. nss - 1 */
   double    **vectors;   /* contains convolved image data */
   double    *krefArea;   /* contains kernel substamp data */
   double    **mat;       /* fitting matrices */
   double    *scprod;     /* kernel sum solution */
   double    sum;         /* sum of fabs, for sigma use in check_stamps */
   double    mean;
   double    median;
   double    mode;        /* sky estimate */
   double    sd;
   double    fwhm;
   double    lfwhm;
   double    chi2;        /* residual in kernel fitting */
   double    norm;        /* kernel sum */
   double    diff;        /* (norm - mean_ksum) * sqrt(sum) */
} stamp_struct;

/* GLOBAL VARS POSSIBLY SET ON COMMAND LINE */
char      *template, *image, *outim;

float     tUThresh, tUKThresh, tLThresh, tGain, tRdnoise, iUThresh, iUKThresh, iLThresh, iGain, iRdnoise;
char      *tNoiseIm, *iNoiseIm, *tMaskIm, *iMaskIm, *kernelImIn, *kernelImOut, *outMask;
float     tPedestal, iPedestal;
int       hwKernel;
float     kerFitThresh, scaleFitThresh, minFracGoodStamps;
float     kfSpreadMask1, kfSpreadMask2;
int       nRegX, nRegY;
char      *regFile;
char      *regKeyWord;
int       numRegKeyWord;
int       nStampY, nStampX, useFullSS;
int       nKSStamps, hwKSStamp;
char      *sstampFile;
int       findSSC;
int       kerOrder, bgOrder;
float     statSig, kerSigReject, kerFracMask;
char      *forceConvolve, *photNormalize, *figMerit;
int       rescaleOK;
float     fillVal, fillValNoise;
char      *noiseImage, *sigmaImage, *convImage;
int       inclNoiseImage, inclSigmaImage, inclConvImage, inclMaskImage, noClobber;
int       doKerInfo, outShort, outNShort;
float     outBzero, outBscale, outNiBzero, outNiBscale;
int       convolveVariance;
int       nThread;

/* GLOBAL VARS NOT SET ON COMMAND LINE */
int       ngauss, *deg_fixe;
float     *sigma_gauss;

//int       rPixX, rPixY;
int       nStamps, nCompKer, nC;

int       nComp, nBGVectors, nCompTotal;

int       fwKernel, fwStamp, fwKSStamp, kcStep;
int       cmpFile;
char      version[32];

/* verbose for debugging */
int        verbose;
/* cmp file stuff */
char       xyfilename[1000];
int        savexyflag;
float      *xcmp,*ycmp;
int        Ncmp;
