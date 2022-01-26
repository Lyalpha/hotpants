#include<stdio.h>
#include<string.h>
#include<math.h>
#include<stdlib.h>
#include<fitsio.h>
#include<time.h>
#include<unistd.h>
#include <pthread.h>

#include "defaults.h"
#include "globals.h"
#include "functions.h"

/* locks for the threads so we don't try to concurrently write to fits files */
pthread_mutex_t writelock = PTHREAD_MUTEX_INITIALIZER;

struct RegionData {
int regNum, nR;
int xMin, yMin, xMax, yMax;
long tNaxes[MAXDIM];
long iNaxes[MAXDIM];
long oNaxes[MAXDIM];
int tBitpix, tNaxis;
int iBitpix, iNaxis;
int oBitpix, oNaxis;
int rXMin, rYMin, rXMax, rYMax;
int kInfoNum;
};

void* processRegion(void *voidData) {

    struct RegionData *regData = voidData;
    int i = regData->regNum;
    int nR = regData->nR;
    int xMin = regData->xMin;
    int yMin = regData->yMin;
    int xMax = regData->xMax;
    int yMax = regData->yMax;
    long tNaxes[MAXDIM];
    long iNaxes[MAXDIM];
    long oNaxes[MAXDIM];
    for (int ii = 0; ii < MAXDIM; ii++) {
        tNaxes[ii] = regData->tNaxes[ii];
        iNaxes[ii] = regData->iNaxes[ii];
        oNaxes[ii] = regData->oNaxes[ii];
    }
    int tNaxis = regData->tNaxis;
    int iNaxis = regData->iNaxis;
    int rXMin = regData->rXMin;
    int rYMin = regData->rYMin;
    int rXMax = regData->rXMax;
    int rYMax = regData->rYMax;
    int kInfoNum = regData->kInfoNum;

    int j, k, l, m;
    int rXBMin, rYBMin, rXBMax, rYBMax;


    double sumKernel;                              /* photometric normalization */
    double meansigSubstamps, scatterSubstamps;      /* mean sigma and scatter of substamps, for fits header */
    double meansigSubstampsF, scatterSubstampsF;    /* mean sigma and scatter of substamps, using final noise */
    int NskippedSubstamps;                      /* Number of skipped substamps, for fits header */

    double diffrat;

    /* fits pointers */
    fitsfile *tPtr;
    fitsfile *iPtr;
    fitsfile *oPtr;
    fitsfile *ePtr;

    /* template image */
    double *tKerSol = NULL;
    float *tRData = NULL;
    stamp_struct *ctStamps = NULL;
    double tMerit = 0;

    /* comparison image */
    double *iKerSol = NULL;
    float *iRData = NULL;
    stamp_struct *ciStamps = NULL;
    double iMerit = 0;

    /* output image */
    float *oRData = NULL;

    /* extraneous image */
    float *eRData = NULL;

    /* extraneous output regions */
    int *misRData = NULL;  /* image    mask */
    int *mtsRData = NULL;  /* template mask */
    int *mRData = NULL;   /* bad input data mask */

    int xBufLo, xBufHi, yBufLo, yBufHi;     /* buffer */

    int sXMin, sYMin, sXMax, sYMax;         /* each stamp */
    int niS, ntS;                           /* stamp counter */
    int nS;                              /* region stamp counter */
    int convTmpl;                           /* which way the convolution goes */


    int buffer_size;                        /* the buffer_size to use, depending on number of regions */
    int       rPixX, rPixY;

    /* cfitsio */
    long pixMin[2], pixMax[2];
    long inc[2] = {1, 1};
    int fpixelOutX, fpixelOutY, lpixelOutX, lpixelOutY;
    int anynul;
    int status = 0;
    char hKeyword[1024], hInfo[2048];


    int nx2norm;
    double sum, mean, median, mode, sd, fwhm, lfwhm, x2norm;
    double summ, meanm, medianm, modem, sdm, fwhmm, lfwhmm;
    double nsum, nmean, nmedian, nmode, nsd, nfwhm, nlfwhm;
    double nsumm, nmeanm, nmedianm, nmodem, nsdm, nfwhmm, nlfwhmm;


    /* misc */
    float iSFrac, tSFrac;
    float regKerFitThresh;             /* Kernel fitting threshhold for this region */
    int flag;
    double inv1;
    int sBorder = hwKSStamp + hwKernel;
    float *temp, *temp2;
    int *indx;

    /* kernel info */
    double **check_mat = NULL, *check_vec = NULL;
    double *check_stack = NULL, **wxy = NULL;
    double *kernel_coeffs, *kernel, **kernel_vec;
    double *filter_x, *filter_y;

    /* malloc data structures */
    fprintf(stderr, "Allocating memory\n");
    if (!(temp = (float *) calloc((fwKSStamp + fwKernel) * fwKSStamp, sizeof(float))) ||
        !(indx = (int *) calloc((nCompTotal + 1 + 100), sizeof(int))) ||

        !(kernel = (double *) calloc(fwKernel * fwKernel, sizeof(double))) ||
        !(kernel_coeffs = (double *) calloc(nCompKer, sizeof(double))) ||
        !(kernel_vec = (double **) calloc(nCompKer, sizeof(double *))) ||
        !(filter_x = (double *) calloc(nCompKer * fwKernel, sizeof(double))) ||
        !(filter_y = (double *) calloc(nCompKer * fwKernel, sizeof(double))) ||

        !(check_mat = (double **) calloc(nC, sizeof(double *))) ||
        !(check_vec = (double *) calloc(nC, sizeof(double))) ||
        !(check_stack = (double *) calloc(nStamps, sizeof(double)))) {
        exit(1);
    }
    for (j = 0; j < nC; j++)
        if (!(check_mat[j] = (double *) calloc(nC, sizeof(double)))) {
        exit(1);
    }

    /* start processing  region */
    if (kernelImIn) {
        /* grab region info from kernelImIn */
        readKernel(kernelImIn, i, &tKerSol, &iKerSol, &rXMin, &rXMax, &rYMin, &rYMax,
                   &meansigSubstamps, &scatterSubstamps,
                   &meansigSubstampsF, &scatterSubstampsF,
                   &diffrat, &NskippedSubstamps);
    } else {
        meansigSubstamps = scatterSubstamps = 0.0;
        NskippedSubstamps = 0;
    }

    if (nR > 1) {
        /* buffered - add an extra half-stamp for merging regions */
        buffer_size = fwStamp / 2;
    } else {
        buffer_size = hwKernel;
    }
    rXBMin = imax(xMin, rXMin - buffer_size);
    rYBMin = imax(yMin, rYMin - buffer_size);
    rXBMax = imin(xMax, rXMax + buffer_size);
    rYBMax = imin(yMax, rYMax + buffer_size);

    /* size of the buffering, used when merging output image sections */
    xBufLo = rXMin - rXBMin;
    xBufHi = rXBMax - rXMax;
    yBufLo = rYMin - rYBMin;
    yBufHi = rYBMax - rYMax;

    /* size of buffered region, in pixels */
    rPixX = rXBMax - rXBMin + 1;
    rPixY = rYBMax - rYBMin + 1;

    /* detemine the limits of the input images */
    /* NOTE: cfitsio wants first pixel = 1, not 0 */
    pixMin[0] = rXBMin + 1;
    pixMin[1] = rYBMin + 1;
    pixMax[0] = rXBMax + 1;
    pixMax[1] = rYBMax + 1;

    /* determine the limits of the output images */
    fpixelOutX = rXBMin + xBufLo + 1;
    fpixelOutY = rYBMin + yBufLo + 1;
    lpixelOutX = fpixelOutX + (rPixX - xBufHi - xBufLo - 1);
    lpixelOutY = fpixelOutY + (rPixY - yBufHi - yBufLo - 1);


    fprintf(stderr, "[Region %d] Pixels            : %ld:%ld,%ld:%ld\n", i, pixMin[0], pixMax[0], pixMin[1],
            pixMax[1]);
    fprintf(stderr, "[Region %d] Vector Indices (buffered) : %d:%d,%d:%d\n", i, rXBMin, rXBMax, rYBMin, rYBMax);
    fprintf(stderr, "[Region %d] Vector Indices (good data): %d:%d,%d:%d\n", i, rXMin, rXMax, rYMin, rYMax);

    /* malloc standard input and output arrays */
    tRData = (float *) calloc(rPixX * rPixY, sizeof(float));
    iRData = (float *) calloc(rPixX * rPixY, sizeof(float));
    oRData = (float *) calloc(rPixX * rPixY, sizeof(float));
    eRData = (float *) calloc(rPixX * rPixY, sizeof(float));
    if (tRData == NULL || iRData == NULL || oRData == NULL || eRData == NULL) {
        fprintf(stderr, "[Region %d] Cannot Allocate Standard Data Arrays\n", i);
        exit(1);
    }

    /* set them to fillVal */
    fset(tRData, fillVal, rPixX, rPixY);
    fset(iRData, fillVal, rPixX, rPixY);

    /* set eRData and oRData to fillValNoise */
    fset(oRData, fillValNoise, rPixX, rPixY);
    fset(eRData, fillValNoise, rPixX, rPixY);

    /* malloc mask arrays */
    mRData = (int *) calloc(rPixX * rPixY, sizeof(int));
    misRData = (int *) calloc(rPixX * rPixY, sizeof(int));
    mtsRData = (int *) calloc(rPixX * rPixY, sizeof(int));
    if (mRData == NULL || misRData == NULL || mtsRData == NULL) {
        fprintf(stderr, "[Region %d] Cannot Allocate Mask Arrays\n", i);
        exit(1);
    }

    /* dont bother if info already from kernel image */
    if (!(kernelImIn)) {
        /* ciStamps contains stamp information about convolving image */
        /* ctStamps contains stamp information about convolving template */
        if (!(strncmp(forceConvolve, "i", 1) == 0)) {
            if (verbose >= 2) fprintf(stderr, "[Region %d] Allocating stamps...\n", i);
            if (!(ctStamps = (stamp_struct *) calloc(nStamps, sizeof(stamp_struct)))) {
                printf("[Region %d] Cannot Allocate Stamp List\n", i);
                exit(1);
            }
            if (allocateStamps(ctStamps, nStamps)) {
                fprintf(stderr, "[Region %d] Cannot Allocate Stamp Vector\n", i);
                exit(1);
            }
            tKerSol = (double *) calloc((nCompTotal + 1), sizeof(double));
        }

        if (!(strncmp(forceConvolve, "t", 1) == 0)) {
            if (!(ciStamps = (stamp_struct *) calloc(nStamps, sizeof(stamp_struct)))) {
                printf("[Region %d] Cannot Allocate Stamp List\n", i);
                exit(1);
            }
            if (allocateStamps(ciStamps, nStamps)) {
                fprintf(stderr, "[Region %d] Cannot Allocate Stamp Vector\n", i);
                exit(1);
            }
            iKerSol = (double *) calloc((nCompTotal + 1), sizeof(double));
        }
    }

    /* open up template */
    if (fits_open_file(&tPtr, template, 0, &status))
        printError(status);

    /* open up comparison image */
    if (fits_open_file(&iPtr, image, 0, &status))
        printError(status);

    /* read image region from fits files into input arrays */
    if (fits_read_subset_flt(tPtr, 0, tNaxis, tNaxes, pixMin, pixMax, inc, 0, tRData, &anynul, &status) ||
        fits_read_subset_flt(iPtr, 0, iNaxis, iNaxes, pixMin, pixMax, inc, 0, iRData, &anynul, &status)) {
        printError(status);
    }
    /* close unneeded files */
    if (fits_close_file(iPtr, &status) || fits_close_file(tPtr, &status)) {
        printError(status);
    }

    /* take off possible pedestal */
    if (tPedestal != 0. || iPedestal != 0.) {
        for (l = rPixX * rPixY; l--;) {
            tRData[l] -= tPedestal;
            iRData[l] -= iPedestal;
        }
    }

    /* we will use eRData and oRData here to make the APPROXIMATE
       noise image - approximate because you don't know which way to
       convolve yet */
    if (iNoiseIm) {
        if (fits_open_file(&ePtr, iNoiseIm, 0, &status))
            printError(status);
        if (fits_read_subset_flt(ePtr, 0, iNaxis, iNaxes, pixMin, pixMax, inc, 0, oRData, &anynul, &status) ||
            fits_close_file(ePtr, &status))
            printError(status);

        for (l = rPixX * rPixY; l--;)
            oRData[l] *= oRData[l];
    } else {
        makeNoiseImage4(oRData, iRData, 1. / iGain, iRdnoise / iGain, rPixX, rPixY);
    }

    if (tNoiseIm) {
        if (fits_open_file(&ePtr, tNoiseIm, 0, &status))
            printError(status);
        if (fits_read_subset_flt(ePtr, 0, tNaxis, tNaxes, pixMin, pixMax, inc, 0, eRData, &anynul, &status) ||
            fits_close_file(ePtr, &status))
            printError(status);

        for (l = rPixX * rPixY; l--;)
            eRData[l] *= eRData[l];
    } else {
        makeNoiseImage4(eRData, tRData, 1. / tGain, tRdnoise / tGain, rPixX, rPixY);
    }

    /* add noise portions, NOT SQRT HERE
       this is only used for getStampSig() */
    for (l = rPixX * rPixY; l--;)
        oRData[l] += eRData[l];

    /* NOTE : oRData is (temporarily) the approximate noise squared per pixel */
    /*        and eRData gets realloced later */

    /*
      mask out bad input pixels, spread by hwKernel * kfSpreadMask1
    */
    if (iMaskIm) {
        if (fits_open_file(&ePtr, iMaskIm, 0, &status))
            printError(status);
        if (fits_read_subset_int(ePtr, 0, tNaxis, tNaxes, pixMin, pixMax, inc, 0, misRData, &anynul, &status) ||
            fits_close_file(ePtr, &status))
            printError(status);

        /* register as read in from input mask */
        for (l = rPixX * rPixY; l--;) {
            misRData[l] |= FLAG_INPUT_MASK * (misRData[l] > 0);
            mRData[l] |= misRData[l];
        }
    }


    if (tMaskIm) {
        if (fits_open_file(&ePtr, tMaskIm, 0, &status))
            printError(status);
        if (fits_read_subset_int(ePtr, 0, tNaxis, tNaxes, pixMin, pixMax, inc, 0, mtsRData, &anynul, &status) ||
            fits_close_file(ePtr, &status))
            printError(status);

        /* register as read in from input mask */
        for (l = rPixX * rPixY; l--;) {
            mtsRData[l] |= FLAG_INPUT_MASK * (mtsRData[l] > 0);
            mRData[l] |= mtsRData[l];
        }
    }

    /* add own bits and spread */
    makeInputMask(tRData, iRData, mRData, rPixX, rPixY);

    /* make sure program can't find a valid stamp center around the edge of the image */
    /* mask with 0x100 & 0x400 */
    for (l = 0; l < rPixY; l++) {
        for (k = 0; k < sBorder; k++) {
            mRData[k + rPixX * l] |= (FLAG_T_BAD | FLAG_I_BAD);
        }
        for (k = rPixX - sBorder; k < rPixX; k++) {
            mRData[k + rPixX * l] |= (FLAG_T_BAD | FLAG_I_BAD);
        }
    }
    for (l = 0; l < sBorder; l++)
        for (k = sBorder; k < rPixX - sBorder; k++)
            mRData[k + rPixX * l] |= (FLAG_T_BAD | FLAG_I_BAD);
    for (l = rPixY - sBorder; l < rPixY; l++)
        for (k = sBorder; k < rPixX - sBorder; k++)
            mRData[k + rPixX * l] |= (FLAG_T_BAD | FLAG_I_BAD);

    status = 1;
    flag = 1;
    regKerFitThresh = kerFitThresh;
    if (!(kernelImIn)) {
        /* only do a single loop for now... */
        while ((status <= 2) && (flag)) {
            flag = 0;
            niS = 0;
            ntS = 0;

            for (l = 0; l < nStampY; l++) {
                for (k = 0; k < nStampX; k++) {

                    fprintf(stderr, "[Region %d] Build stamp  : t %4d i %4d (grid coord %2d %2d)\n", i, ntS, niS, k, l);
                    /* coordinates in the image, not the region */

                    /* NOTE : we keep float-valued 'rPixX / nStampX' in
                       here to not exclude a full stamps' width of space on
                       the right and upper sides of the image.  Allows the
                       stamps to be separated by a little bit but fills the
                       region more evenly. */

                    sXMin = rXBMin + k * rPixX / nStampX;
                    sYMin = rYBMin + l * rPixY / nStampY;
                    sXMax = imin(sXMin + fwStamp - 1, rXBMax);
                    sYMax = imin(sYMin + fwStamp - 1, rYBMax);

                    /* armin: reinitializing sscnt and nss? */
                    if (!(strncmp(forceConvolve, "i", 1) == 0)) {
                        ctStamps[ntS].sscnt = ctStamps[ntS].nss = 0;
                        ctStamps[ntS].chi2 = 0.0;
                    }
                    if (!(strncmp(forceConvolve, "t", 1) == 0)) {
                        ciStamps[niS].sscnt = ciStamps[niS].nss = 0;
                        ciStamps[niS].chi2 = 0.0;
                    }

                    if (sstampFile) {
                        if (verbose >= 2) fprintf(stderr, "[Region %d] Adding centers manually\n", i);
                        for (m = 0; m < Ncmp; m++) {
                            if ((xcmp[m] > sXMin + hwKernel + 1) && (xcmp[m] < sXMax - hwKernel - 1) &&
                                (ycmp[m] > sYMin + hwKernel + 1) && (ycmp[m] < sYMax - hwKernel - 1)) {

                                buildStamps(sXMin, sXMax, sYMin, sYMax, &niS, &ntS, 0, rXBMin, rYBMin,
                                            ciStamps, ctStamps, iRData, tRData,
                                            xcmp[m] - rXBMin, ycmp[m] - rYBMin, mRData, rPixX, rPixY);

                            }
                        }
                        if (findSSC) {
                            if (verbose >= 2)
                                fprintf(stderr, "[Region %d] Automatically finding additional centers\n", i);
                            /* and then choose enough stamps manually to fill the array */
                            buildStamps(sXMin, sXMax, sYMin, sYMax, &niS, &ntS, 1, rXBMin, rYBMin,
                                        ciStamps, ctStamps, iRData, tRData, 0, 0, mRData, rPixX, rPixY);
                        }
                    } else {
                        if (useFullSS) {
                            /* don't try and find centers, just use center of stamp */
                            buildStamps(sXMin, sXMax, sYMin, sYMax, &niS, &ntS, 0, rXBMin, rYBMin,
                                        ciStamps, ctStamps, iRData, tRData, 0, 0, mRData, rPixX, rPixY);
                        } else {
                            buildStamps(sXMin, sXMax, sYMin, sYMax, &niS, &ntS, 1, rXBMin, rYBMin,
                                        ciStamps, ctStamps, iRData, tRData, 0, 0, mRData, rPixX, rPixY);
                        }
                    }


                    if (!(strncmp(forceConvolve, "i", 1) == 0)) {
                        if (verbose >= 2)
                            fprintf(stderr, "[Region %d]     templ: %d substamps\n", i, ctStamps[ntS].nss);
                        if (ctStamps[ntS].nss > 0) ntS += 1;
                    }

                    if (!(strncmp(forceConvolve, "t", 1) == 0)) {
                        if (verbose >= 2)
                            fprintf(stderr, "[Region %d]     image: %d substamps\n", i, ciStamps[niS].nss);
                        if (ciStamps[niS].nss > 0) niS += 1;
                    }

                }
            }

            iSFrac = niS / (float) nStamps;
            tSFrac = ntS / (float) nStamps;

            if (strncmp(forceConvolve, "i", 1) == 0) {
                fprintf(stderr, "[Region %d] %d stamps built (%.2f%s)\n\n", i, niS, iSFrac, "%");
                if (iSFrac < minFracGoodStamps)
                    flag = 1;
            } else if (strncmp(forceConvolve, "t", 1) == 0) {
                fprintf(stderr, "[Region %d] %d stamps built (%.2f%s)\n\n", i, ntS, tSFrac, "%");
                if (tSFrac < minFracGoodStamps)
                    flag = 1;
            } else if ((iSFrac < minFracGoodStamps) || (tSFrac < minFracGoodStamps)) {
                fprintf(stderr, "[Region %d] %d and %d stamps built (%.2f%s, %.2f%s)\n\n", i, ntS, niS, tSFrac, "%",
                        iSFrac, "%");
                flag = 1;
            } else {
                fprintf(stderr, "[Region %d] %d and %d stamps built (%.2f%s, %.2f%s)\n\n", i, ntS, niS, tSFrac, "%",
                        iSFrac, "%");
                break;
            }

            /*
            success here requires that you have chosen the number of
            stamps such that you expect 1 object at good S/N in each
            one.  if too few stamps have been chosen, either

            A : asking for too many stamps

            B : too few objects at high S/N

            i think B is the more likely scenario, so we lower the
            fitting threshold.  do it only once now to be safe,
            since this is a new feature brought on by lack on sleep
            during Chilean obsreving runs...
            */

            if ((flag) && (status <= 1) && (scaleFitThresh < 1.)) {

                /* start over... */
                regKerFitThresh *= scaleFitThresh;

                fprintf(stderr, "[Region %d] Too few stamps were fit, scaling down fitting threshold to %.2f\n", i, regKerFitThresh);

                if (ctStamps) {
                    freeStampMem(ctStamps, nStamps);
                    free(ctStamps);
                    ctStamps = NULL;
                }
                if (ciStamps) {
                    freeStampMem(ciStamps, nStamps);
                    free(ciStamps);
                    ciStamps = NULL;
                }


                /* sheesh, reallocate stamps ... */

                if (!(strncmp(forceConvolve, "i", 1) == 0)) {
                    if (!(ctStamps = (stamp_struct *) calloc(nStamps, sizeof(stamp_struct)))) {
                        printf("[Region %d] Cannot Allocate Stamp List\n", i);
                        exit(1);
                    }
                    if (allocateStamps(ctStamps, nStamps)) {
                        fprintf(stderr, "[Region %d] Cannot Allocate Stamp Vector\n", i);
                        exit(1);
                    }
                }

                if (!(strncmp(forceConvolve, "t", 1) == 0)) {
                    if (!(ciStamps = (stamp_struct *) calloc(nStamps, sizeof(stamp_struct)))) {
                        printf("Cannot Allocate Stamp List\n");
                        exit(1);
                    }
                    if (allocateStamps(ciStamps, nStamps)) {
                        fprintf(stderr, "[Region %d] Cannot Allocate Stamp Vector\n", i);
                        exit(1);
                    }
                }

                /* allow getpsfcenter to use pixels again */
                for (l = rPixX * rPixY; l--;)
                    mRData[l] &= ~0xa00;

            }
            /* loop */
            status += 1;
        }
        status = 0;

        /* no stamps at all? */
        /* FIXME We are not in a region loop here - is it ok to just return here if these conditions met? */
        if ((niS == 0) && (ntS == 0)) {
            fprintf(stderr, "[Region %d] No stamps found!", i);
            return 0;
        }
        if (strncmp(forceConvolve, "i", 1) == 0)
            if (niS == 0) {
                fprintf(stderr, "[Region %d] No image stamps found!", i);
                return 0;
            }
        if (strncmp(forceConvolve, "t", 1) == 0)
            if (ntS == 0) {
                fprintf(stderr, "[Region %d] No template stamps found!", i);
                return 0;
            }

        /*

        NOTE ON MASKS : mRData represents mask such that
        tLThresh < tpix < tUThresh
        iLThresh < ipix < iUThresh
        Spread by hwKernel * kfSpreadMask1
        and if -inm, add to mRData
        */

        /* initialize kernel weight mask: kernel_vec */
        getKernelVec(filter_x, filter_y, kernel_vec);
        fprintf(stderr, "[Region %d] Filling Template sub-stamps\n", i);
        /* fit the kernel going each way unless told otherwise */
        if (!(strncmp(forceConvolve, "i", 1) == 0)) {
            for (k = 0; k < ntS; k++) {
                ctStamps[k].sscnt = 0;
                fillStamp(&ctStamps[k], tRData, iRData, mRData, rPixX, rPixY, filter_x, filter_y, temp);
            }
            if ((strncmp(forceConvolve, "b", 1) == 0)) {
                fprintf(stderr, "\n\n[Region %d] Trying to convolve the TEMPLATE to fit IMAGE\n", i);
                tMerit = check_stamps(ctStamps, ntS, iRData, oRData, mRData, rPixX, rPixY, check_mat, check_vec, check_stack, wxy, kernel_coeffs, kernel, kernel_vec, indx, temp);
                fprintf(stderr, "[Region %d]     Result : merit = %.3f\n", i, tMerit);
            } else
                tMerit = iMerit = 0;
        }

        fprintf(stderr, "[Region %d] Filling Image sub-stamps\n", i);
        if (!(strncmp(forceConvolve, "t", 1) == 0)) {

            for (k = 0; k < niS; k++) {
                ciStamps[k].sscnt = 0;
                fillStamp(&ciStamps[k], iRData, tRData, mRData, rPixX, rPixY, filter_x, filter_y, temp);
            }
            if ((strncmp(forceConvolve, "b", 1) == 0)) {
                fprintf(stderr, "\n\n[Region %d] Trying to convolve the IMAGE to fit TEMPLATE \n", i);
                iMerit = check_stamps(ciStamps, niS, tRData, oRData, mRData, rPixX, rPixY, check_mat, check_vec, check_stack, wxy, kernel_coeffs, kernel, kernel_vec, indx, temp);
                fprintf(stderr, "[Region %d]     Result : merit = %.3f\n", i, iMerit);
            } else
                iMerit = tMerit = 0;
        }

    } /* end of if not kernelImIn */
    else {
        /* initialize kernel weight mask: kernel_vec */
        status = 0;
        getKernelVec(filter_x, filter_y, kernel_vec);
    }

    /* decide which way to go */
    /* perhaps for good! */
    if ((strncmp(forceConvolve, "t", 1) == 0) ||
        ((tMerit < iMerit) && (!(strncmp(forceConvolve, "i", 1) == 0)))) {
        convTmpl = 1;
        if (ciStamps) free(ciStamps);
        ciStamps = NULL;

    } else {
        convTmpl = 0;
        if (ctStamps) free(ctStamps);
        ctStamps = NULL;
    }

    /* do it! */
    if (convTmpl) {
        fprintf(stderr, "\n\n[Region %d] Region %d:%d,%d:%d : Convolving TEMPLATE\n", i, rXMin, rXMax, rYMin, rYMax);

        freeStampMem(ciStamps, nStamps);
        /*allocateStamps(ciStamps, nStamps);*/
        nS = ntS;

        /* inside fitKernel, mRData is still used as the input bad pixel mask */
        /* and oRData is the input noise^2 */
        if (!(kernelImIn))
            fitKernel(ctStamps, iRData, tRData, oRData, tKerSol, &meansigSubstamps, &scatterSubstamps,
                      &NskippedSubstamps, mRData, rPixX, rPixY, nS, filter_x, filter_y, indx, temp);

        /* zero out oRData to accept output diff image */
        oRData = (float *) realloc(oRData, rPixX * rPixY * sizeof(float));
        fset(oRData, fillVal, rPixX, rPixY);

        /* re-do the output mask so that only convolved mask is spread */
        /* mtsRData is the mask which is naturally spread by the convolution */
        /* add in my own bits */
        for (l = rPixX * rPixY; l--;) {
            mtsRData[l] |= (FLAG_INPUT_ISBAD | FLAG_BAD_PIXVAL) * (tRData[l] == fillVal);
            mtsRData[l] |= (FLAG_INPUT_ISBAD | FLAG_SAT_PIXEL) * (tRData[l] >= tUThresh);
            mtsRData[l] |= (FLAG_INPUT_ISBAD | FLAG_LOW_PIXEL) * (tRData[l] <= tLThresh);
        }

        /* start mRData over as output mask */
        memset(mRData, 0, rPixX * rPixY * sizeof(int));

        /* make the template noise image to send to spatial_convolve */
        eRData = (float *) realloc(eRData, rPixX * rPixY * sizeof(float));
        fset(eRData, fillValNoise, rPixX, rPixY);

        if (tNoiseIm) {
            if (fits_open_file(&ePtr, tNoiseIm, 0, &status))
                printError(status);
            if (fits_read_subset_flt(ePtr, 0, tNaxis, tNaxes, pixMin, pixMax, inc, 0, eRData, &anynul, &status) ||
                fits_close_file(ePtr, &status))
                printError(status);

            for (l = rPixX * rPixY; l--;)
                eRData[l] *= eRData[l];
        } else {
            makeNoiseImage4(eRData, tRData, 1. / tGain, tRdnoise / tGain, rPixX, rPixY);
        }

        /* spatial_convolve effectively spreads the input mtsRData mask into global mRData output mask!  bitwise... */
        fprintf(stderr, "\n[Region %d] Convolving...\n", i);

        spatial_convolve(tRData, &eRData, rPixX, rPixY, tKerSol, oRData, mtsRData, mRData, rPixX, rPixY, kernel_coeffs, kernel, kernel_vec);

        /* correct for background */
        for (l = hwKernel; l < rPixY - hwKernel; l++)
            for (k = hwKernel; k < rPixX - hwKernel; k++)
                oRData[k + rPixX * l] += get_background(k, l, tKerSol, rPixX, rPixY);


        sumKernel = make_kernel(rXMin, rYMin, tKerSol, rPixX, rPixY, kernel_coeffs, kernel, kernel_vec);
        fprintf(stderr, "[Region %d] Sum Kernel at %d,%d: %f\n", i, rXMin, rYMin, sumKernel);
        sumKernel = make_kernel(rXMax, rYMax, tKerSol, rPixX, rPixY, kernel_coeffs, kernel, kernel_vec);
        fprintf(stderr, "[Region %d] Sum Kernel at %d,%d: %f\n", i, rXMax, rYMax, sumKernel);
        /* use middle of region to normalize image */
        sumKernel = make_kernel(rPixX / 2, rPixY / 2, tKerSol, rPixX, rPixY, kernel_coeffs, kernel, kernel_vec);
        fprintf(stderr, "[Region %d] Using Kernel Sum = %f\n\n", i, sumKernel);

        /* eRData now contains partial noise image */
        /* oRData now contains difference image */
        /* tRData is no longer needed */
        /* use to read in other noise image if necessary */
        fset(tRData, fillValNoise, rPixX, rPixY);
        if (iNoiseIm) {
            if (fits_open_file(&ePtr, iNoiseIm, 0, &status))
                printError(status);
            if (fits_read_subset_flt(ePtr, 0, iNaxis, iNaxes, pixMin, pixMax, inc, 0, tRData, &anynul, &status) ||
                fits_close_file(ePtr, &status))
                printError(status);

            for (l = rPixX * rPixY; l--;)
                tRData[l] *= tRData[l];
        } else {
            makeNoiseImage4(tRData, iRData, 1. / iGain, iRdnoise / iGain, rPixX, rPixY);
        }

        /* add noise portions to create final noise image */
        for (l = rPixX * rPixY; l--;)
            tRData[l] = sqrt(tRData[l] + eRData[l]);
        free(eRData);
        eRData = NULL;

        /* finally, add in NON-convolved input mask */
        /* do it AFTER spatial_convolve so that it is not spread */
        for (l = rPixX * rPixY; l--;) {
            mRData[l] |= (FLAG_OUTPUT_ISBAD | FLAG_INPUT_ISBAD | FLAG_BAD_PIXVAL) * (iRData[l] == fillVal);
            mRData[l] |= (FLAG_OUTPUT_ISBAD | FLAG_INPUT_ISBAD | FLAG_SAT_PIXEL) * (iRData[l] >= iUThresh);
            mRData[l] |= (FLAG_OUTPUT_ISBAD | FLAG_INPUT_ISBAD | FLAG_LOW_PIXEL) * (iRData[l] <= iLThresh);
            mRData[l] |= misRData[l];
            mRData[l] |= FLAG_OUTPUT_ISBAD * ((misRData[l] & FLAG_INPUT_ISBAD) > 0);
        }

        if (savexyflag) {
            savexy(ctStamps, ntS, rXBMin, rYBMin, i);
        }


        /*
        NOTE : the following arrays will store the following info
        oRData   : Convolved input template image AND difference image AND sigma image
        tRData   : Noise image
        iRData   : Input image
        */

    } else {
        fprintf(stderr, "\n\n[Region %d] Region %d:%d,%d:%d : Convolving IMAGE\n", i, rXMin, rXMax, rYMin, rYMax);

        freeStampMem(ctStamps, nStamps);
        /*allocateStamps(ctStamps, nStamps);*/
        nS = niS;

        /* inside fitKernel, mRData is still used as the input bad pixel mask */
        /* and oRData is the input noise^2 */
        if (!(kernelImIn))
            fitKernel(ciStamps, tRData, iRData, oRData, iKerSol, &meansigSubstamps, &scatterSubstamps,
                      &NskippedSubstamps, mRData, rPixX, rPixY, nS, filter_x, filter_y, indx, temp);

        /* zero out oRData to accept output diff image */
        oRData = (float *) realloc(oRData, rPixX * rPixY * sizeof(float));
        fset(oRData, fillVal, rPixX, rPixY);

        /* re-do the output mask so that only convolved mask is spread */
        /* misRData is the mask which is naturally spread by the convolution */
        /* add in my own bits */
        for (l = rPixX * rPixY; l--;) {
            misRData[l] |= (FLAG_INPUT_ISBAD | FLAG_BAD_PIXVAL) * (iRData[l] == fillVal);
            misRData[l] |= (FLAG_INPUT_ISBAD | FLAG_SAT_PIXEL) * (iRData[l] >= iUThresh);
            misRData[l] |= (FLAG_INPUT_ISBAD | FLAG_LOW_PIXEL) * (iRData[l] <= iLThresh);
        }

        /* start mRData over as output mask */
        memset(mRData, 0, rPixX * rPixY * sizeof(int));

        /* make the template noise image to send to spatial_convolve */
        eRData = (float *) realloc(eRData, rPixX * rPixY * sizeof(float));
        fset(eRData, fillValNoise, rPixX, rPixY);

        if (iNoiseIm) {
            if (fits_open_file(&ePtr, iNoiseIm, 0, &status))
                printError(status);
            if (fits_read_subset_flt(ePtr, 0, iNaxis, iNaxes, pixMin, pixMax, inc, 0, eRData, &anynul, &status) ||
                fits_close_file(ePtr, &status))
                printError(status);

            for (l = rPixX * rPixY; l--;)
                eRData[l] *= eRData[l];
        } else {
            makeNoiseImage4(eRData, iRData, 1. / iGain, iRdnoise / iGain, rPixX, rPixY);
        }

        /* spatial_convolve effectively spreads the input misRData mask into global mRData output mask!  bitwise... */
        fprintf(stderr, "\n[Region %d]  Convolving...\n", i);
        spatial_convolve(iRData, &eRData, rPixX, rPixY, iKerSol, oRData, misRData, mRData, rPixX, rPixY, kernel_coeffs, kernel, kernel_vec);

        /* correct for background */
        for (l = hwKernel; l < rPixY - hwKernel; l++)
            for (k = hwKernel; k < rPixX - hwKernel; k++)
                oRData[k + rPixX * l] += get_background(k, l, iKerSol, rPixX, rPixY);

        sumKernel = make_kernel(rXMin, rYMin, iKerSol, rPixX, rPixY, kernel_coeffs, kernel, kernel_vec);
        fprintf(stderr, "[Region %d]  Sum Kernel at %d,%d: %f\n", i, rXMin, rYMin, sumKernel);
        sumKernel = make_kernel(rXMax, rYMax, iKerSol, rPixX, rPixY, kernel_coeffs, kernel, kernel_vec);
        fprintf(stderr, "[Region %d]  Sum Kernel at %d,%d: %f\n", i, rXMax, rYMax, sumKernel);
        /* use middle of region to normalize image */
        sumKernel = make_kernel(rPixX / 2, rPixY / 2, iKerSol, rPixX, rPixY, kernel_coeffs, kernel, kernel_vec);
        fprintf(stderr, "[Region %d]  Using Kernel Sum = %f\n\n", i, sumKernel);


        /* eRData now contains partial noise image */
        /* oRData now contains difference image */
        /* iRData is no longer needed */
        /* use to read in other noise image if necessary */
        fset(iRData, fillValNoise, rPixX, rPixY);
        if (tNoiseIm) {
            if (fits_open_file(&ePtr, tNoiseIm, 0, &status))
                printError(status);
            if (fits_read_subset_flt(ePtr, 0, tNaxis, tNaxes, pixMin, pixMax, inc, 0, iRData, &anynul, &status) ||
                fits_close_file(ePtr, &status))
                printError(status);

            for (l = rPixX * rPixY; l--;)
                iRData[l] *= iRData[l];
        } else {
            makeNoiseImage4(iRData, tRData, 1. / tGain, tRdnoise / tGain, rPixX, rPixY);
        }

        /* add noise portions to create final noise image */
        for (l = rPixX * rPixY; l--;)
            iRData[l] = sqrt(iRData[l] + eRData[l]);
        free(eRData);
        eRData = NULL;

        /* finally, add in NON-convolved input mask */
        /* do it AFTER spatial_convolve so that it is not spread */
        for (l = rPixX * rPixY; l--;) {
            mRData[l] |= (FLAG_OUTPUT_ISBAD | FLAG_INPUT_ISBAD | FLAG_BAD_PIXVAL) * (tRData[l] == fillVal);
            mRData[l] |= (FLAG_OUTPUT_ISBAD | FLAG_INPUT_ISBAD | FLAG_SAT_PIXEL) * (tRData[l] >= tUThresh);
            mRData[l] |= (FLAG_OUTPUT_ISBAD | FLAG_INPUT_ISBAD | FLAG_LOW_PIXEL) * (tRData[l] <= tLThresh);
            mRData[l] |= mtsRData[l];
            mRData[l] |= FLAG_OUTPUT_ISBAD * ((mtsRData[l] & FLAG_INPUT_ISBAD) > 0);
        }

        if (savexyflag) {
            savexy(ciStamps, niS, pixMin[0], pixMin[1], i);
        }


        /*
        NOTE : the following arrays will store the following info
        oRData   : Convolved input image AND difference image AND sigma image
        tRData   : Input template image
        iRData   : Noise image
        */
    }

    /* mask edge of mask image */
    for (l = 0; l < rPixY; l++) {
        for (k = 0; k < hwKernel; k++)
            mRData[k + rPixX * l] |= FLAG_OUTPUT_ISBAD;
        for (k = rPixX - hwKernel; k < rPixX; k++)
            mRData[k + rPixX * l] |= FLAG_OUTPUT_ISBAD;
    }
    for (l = 0; l < hwKernel; l++)
        for (k = hwKernel; k < rPixX - hwKernel; k++)
            mRData[k + rPixX * l] |= FLAG_OUTPUT_ISBAD;
    for (l = rPixY - hwKernel; l < rPixY; l++)
        for (k = hwKernel; k < rPixX - hwKernel; k++)
            mRData[k + rPixX * l] |= FLAG_OUTPUT_ISBAD;

    fprintf(stderr, "[Region %d]  Creating and writing output images...\n", i);

    /*
      NOTE : this is split up the way it is so that the code is not
      such a memory hog.  Ideally we'd have arrays for the output
      data, convolved data, noise data, sigma data, etc...  Takes up
      a HUGE amount of memory though.  So the code gets a bit more
      unwieldy. Cost cutting is not pretty sometimes.  Hard lesson
      of life #17 learned at Lucent Technologies.
    */

    /*
      NOTE : we are NOT masking out the convolved image
    */
    inv1 = 1. / sumKernel;

    /* take a lock around all write operations to prevent multi-thread access to fits files */
    fprintf(stderr,"[Region %d] Acquiring writing lock\n", i);
    pthread_mutex_lock(&writelock);
    fprintf(stderr,"[Region %d] Acquired writing lock\n", i);

    /* open up output image */
    if (fits_open_file(&oPtr, outim, 1, &status))
        printError(status);

    if (inclConvImage || convImage) {
        /* renormalize by kernel if necessary */
        for (l = hwKernel; l < rPixY - hwKernel; l++) {
            for (k = hwKernel; k < rPixX - hwKernel; k++) {
                if ((strncmp(photNormalize, "u", 1) != 0) &&
                    ((convTmpl && strncmp(photNormalize, "t", 1) == 0) ||
                     (!convTmpl && strncmp(photNormalize, "i", 1) == 0)))
                    oRData[k + rPixX * l] *= inv1;
            }
        }

        if (inclConvImage) {
            /* if asked for, the convolved image is the second in the multi-fits layer */
            if (fits_movabs_hdu(oPtr, 2, NULL, &status))
                printError(status);
            if (outShort)
                if (fits_set_bscale(oPtr, outBscale, outBzero, &status))
                    printError(status);
            if (hp_fits_write_subset(oPtr, 0, 2, oNaxes, oRData, &status,
                                     outShort, outBzero, outBscale, fpixelOutX, fpixelOutY,
                                     lpixelOutX, lpixelOutY, xBufLo, yBufLo, mRData, rPixX, rPixY))
                printError(status);
        }

        if (convImage) {
            if (fits_open_file(&ePtr, convImage, 1, &status))
                printError(status);
            if (outShort)
                if (fits_set_bscale(ePtr, outBscale, outBzero, &status))
                    printError(status);

            if (hp_fits_write_subset(ePtr, 0, 2, oNaxes, oRData, &status,
                                     outShort, outBzero, outBscale, fpixelOutX, fpixelOutY,
                                     lpixelOutX, lpixelOutY, xBufLo, yBufLo, mRData, rPixX, rPixY) ||
                fits_close_file(ePtr, &status))
                printError(status);
        }

        /* undo above - ugly but speed vs. memory issue (vs. lack of cleverness) ... */
        for (l = hwKernel; l < rPixY - hwKernel; l++) {
            for (k = hwKernel; k < rPixX - hwKernel; k++) {
                if ((strncmp(photNormalize, "u", 1) != 0) &&
                    ((convTmpl && strncmp(photNormalize, "t", 1) == 0) ||
                     (!convTmpl && strncmp(photNormalize, "i", 1) == 0)))
                    oRData[k + rPixX * l] *= sumKernel;
            }
        }
    }

    /* next, reuse the convolved image array for the difference image array */
    if (convTmpl) {
        /* do subtraction, output is not yet masked */
        for (l = hwKernel; l < rPixY - hwKernel; l++) {
            for (k = hwKernel; k < rPixX - hwKernel; k++) {
                m = k + rPixX * l;
                /* do subtraction, in e- */
                oRData[m] -= iRData[m];
                /* renormalize if necessary */
                if ((strncmp(photNormalize, "u", 1) != 0) && (strncmp(photNormalize, "t", 1) == 0)) {
                    oRData[m] *= inv1;

                    /* NOTE TO AUTHOR - we rescale noise here, not in makeNoiseImage */
                    tRData[m] *= inv1;
                }
                oRData[m] *= -1.;
            }
        }

        /* get SSSIG and SSSCAT from final kernel solution! */
        if (!(kernelImIn)) {
            if (strncmp(figMerit, "v", 1) == 0) {
                temp2 = (float *) calloc(nS, sizeof(float));
                k = 0;
                for (l = 0; l < nS; l++) {
                    /* if was fit with a good legit substamp */
                    if (ctStamps[l].sscnt < ctStamps[l].nss) {
                        getFinalStampSig(&ctStamps[l], oRData, tRData, &sum, mRData, rPixX, rPixY);
                        temp2[k++] = sum;
                        /*fprintf(stderr, "SSSIG %d : %f\n", k, sum);*/
                    }
                }
                /* save the mean and scatter so that it can be saved in the fits header */
                sigma_clip(temp2, k, &meansigSubstampsF, &scatterSubstampsF, 10);
                fprintf(stderr, "[Region %d]   FINAL Mean sig: %6.3f stdev: %6.3f\n", i, meansigSubstampsF, scatterSubstampsF);
                free(temp2);
            }
        }
        freeStampMem(ctStamps, nStamps);

    } else {
        /* do subtraction, output is not yet masked */
        for (l = hwKernel; l < rPixY - hwKernel; l++) {
            for (k = hwKernel; k < rPixX - hwKernel; k++) {
                m = k + rPixX * l;
                /* do subtraction, in e- */
                oRData[m] -= tRData[m];
                /* renormalize if necessary */
                if ((strncmp(photNormalize, "u", 1) != 0) && (strncmp(photNormalize, "i", 1) == 0)) {
                    oRData[m] *= inv1;

                    /* NOTE TO AUTHOR - we rescale noise here, not in makeNoiseImage */
                    iRData[m] *= inv1;
                }
            }
        }

        /* get SSSIG and SSSCAT from final kernel solution! */
        if (!(kernelImIn)) {
            if (strncmp(figMerit, "v", 1) == 0) {
                temp2 = (float *) calloc(nS, sizeof(float));
                k = 0;
                for (l = 0; l < nS; l++) {
                    /* if was fit with a good legit substamp */
                    if (ciStamps[l].sscnt < ciStamps[l].nss) {
                        getFinalStampSig(&ciStamps[l], oRData, iRData, &sum, mRData, rPixX, rPixY);
                        temp2[k++] = sum;
                        /*fprintf(stderr, "SSSIG %d : %f\n", k, sum);*/
                    }
                }
                /* save the mean and scatter so that it can be saved in the fits header */
                sigma_clip(temp2, k, &meansigSubstampsF, &scatterSubstampsF, 10);
                fprintf(stderr, "[Region %d]    FINAL Mean sig: %6.3f stdev: %6.3f\n", i, meansigSubstampsF, scatterSubstampsF);
                free(temp2);
            }
        }
        freeStampMem(ciStamps, nStamps);

    }

    fprintf(stderr, "[Region %d]  Getting diffim stats for GOOD pixels : \n", i);
    getStampStats3(oRData, 0, 0, rPixX, rPixY,
                   &sum, &mean, &median,
                   &mode, &sd, &fwhm, &lfwhm, 0x0, 0xffff, 5, mRData, rPixX, rPixY);
    fprintf(stderr, "[Region %d]    Mean   : %.2f\n", i, mean);
    fprintf(stderr, "[Region %d]    Median : %.2f\n", i, median);
    fprintf(stderr, "[Region %d]    Mode   : %.2f\n", i, mode);
    fprintf(stderr, "[Region %d]    Stdev  : %.2f\n", i, sd);
    if (verbose >= 2) fprintf(stderr, "[Region %d]    FWHM   : %.2f\n", i, fwhm);
    if (verbose >= 2) fprintf(stderr, "[Region %d]    lFWHM  : %.2f\n", i, lfwhm);

    fprintf(stderr, "[Region %d]  Getting noiseim stats for GOOD pixels : \n", i);
    if (convTmpl) {
        getStampStats3(tRData, 0, 0, rPixX, rPixY,
                       &nsum, &nmean, &nmedian,
                       &nmode, &nsd, &nfwhm, &nlfwhm, 0x0, 0xffff, 5, mRData, rPixX, rPixY);
        getNoiseStats3(oRData, tRData, &x2norm, &nx2norm, 0x0, 0xffff, mRData, rPixX, rPixY);
    } else {
        getStampStats3(iRData, 0, 0, rPixX, rPixY,
                       &nsum, &nmean, &nmedian,
                       &nmode, &nsd, &nfwhm, &nlfwhm, 0x0, 0xffff, 5, mRData, rPixX, rPixY);
        getNoiseStats3(oRData, iRData, &x2norm, &nx2norm, 0x0, 0xffff, mRData, rPixX, rPixY);
    }

    fprintf(stderr, "[Region %d]    Mean   : %.2f\n", i, nmean);
    fprintf(stderr, "[Region %d]    Median : %.2f\n", i, nmedian);
    fprintf(stderr, "[Region %d]    Mode   : %.2f\n", i, nmode);
    fprintf(stderr, "[Region %d]    Stdev  : %.2f\n", i, nsd);
    if (verbose >= 2) fprintf(stderr, "[Region %d]    FWHM   : %.2f\n", i, nfwhm);
    if (verbose >= 2) fprintf(stderr, "[Region %d]    lFWHM  : %.2f\n", i, nlfwhm);

    /* find ratio for GOOD pixels only */
    if (!(kernelImIn))
        fprintf(stderr, "[Region %d]  Empirical / Expected Noise for GOOD pixels = %.2f\n", i, sd / nmean);
    fprintf(stderr, "[Region %d]  X2NORM = %.2f\n\n", i, x2norm);

    if (!(kernelImIn)) {
        diffrat = sd / nmean;
    }

    fprintf(stderr, "[Region %d] Getting diffim stats for OK pixels : \n", i);
    getStampStats3(oRData, 0, 0, rPixX, rPixY,
                   &summ, &meanm, &medianm,
                   &modem, &sdm, &fwhmm, &lfwhmm, 0xff, FLAG_OUTPUT_ISBAD, 5, mRData, rPixX, rPixY);
    fprintf(stderr, "[Region %d]    Mean   : %.2f\n", i, meanm);
    fprintf(stderr, "[Region %d]    Median : %.2f\n", i, medianm);
    fprintf(stderr, "[Region %d]    Mode   : %.2f\n", i, modem);
    fprintf(stderr, "[Region %d]    Stdev  : %.2f\n", i, sdm);
    if (verbose >= 2) fprintf(stderr, "[Region %d]    FWHM   : %.2f\n", i, fwhmm);
    if (verbose >= 2) fprintf(stderr, "[Region %d]    lFWHM  : %.2f\n", i, lfwhmm);

    fprintf(stderr, "[Region %d]  Getting noiseim stats for OK pixels : \n", i);
    if (convTmpl)
        getStampStats3(tRData, 0, 0, rPixX, rPixY,
                       &nsumm, &nmeanm, &nmedianm,
                       &nmodem, &nsdm, &nfwhmm, &nlfwhmm, 0xff, FLAG_OUTPUT_ISBAD, 5, mRData, rPixX, rPixY);
    else
        getStampStats3(iRData, 0, 0, rPixX, rPixY,
                       &nsumm, &nmeanm, &nmedianm,
                       &nmodem, &nsdm, &nfwhmm, &nlfwhmm, 0xff, FLAG_OUTPUT_ISBAD, 5, mRData, rPixX, rPixY);

    fprintf(stderr, "[Region %d]    Mean   : %.2f\n", i, nmeanm);
    fprintf(stderr, "[Region %d]    Median : %.2f\n", i, nmedianm);
    fprintf(stderr, "[Region %d]    Mode   : %.2f\n", i, nmodem);
    fprintf(stderr, "[Region %d]    Stdev  : %.2f\n", i, nsdm);
    if (verbose >= 2) fprintf(stderr, "[Region %d]    FWHM   : %.2f\n", i, nfwhmm);
    if (verbose >= 2) fprintf(stderr, "[Region %d]    lFWHM  : %.2f\n", i, nlfwhmm);

    /* find ratio for GOOD pixels only */
    if (!(kernelImIn))
        fprintf(stderr, "[Region %d]  Empirical / Expected Noise for OK pixels = %.2f\n\n", i, sdm / nmeanm);

    /* scale noise in OK pixels so that the ratios are equal! */
    if (rescaleOK) {
        if (!(kernelImIn))
            diffrat = (sdm / nmeanm) / diffrat;
        else
            /* its read in from the input image! */
            ;

        if (diffrat > 1) {
            fprintf(stderr, " Scale OK pixel noise by = %.2f\n", diffrat);
            if (convTmpl)
                for (l = rPixX * rPixY; l--;) {
                    if ((mRData[l] & 0xff) && (!(mRData[l] & FLAG_OUTPUT_ISBAD)))
                        tRData[l] *= diffrat;
                }
            else
                for (l = rPixX * rPixY; l--;) {
                    if ((mRData[l] & 0xff) && (!(mRData[l] & FLAG_OUTPUT_ISBAD)))
                        iRData[l] *= diffrat;
                }
        } else
            fprintf(stderr, " Leave OK pixel noise as-is\n");
    }

    /* BAD pixels are nothing but bad... */

    /* mask the output images */
    if (kfSpreadMask2 >= 0) {
        if (convTmpl) {
            for (l = hwKernel; l < rPixY - hwKernel; l++) {
                for (k = hwKernel; k < rPixX - hwKernel; k++) {
                    if (mRData[k + rPixX * l] & FLAG_OUTPUT_ISBAD) {
                        oRData[k + rPixX * l] = fillVal;
                        tRData[k + rPixX * l] = fillValNoise;
                    }
                }
            }
        } else {
            for (l = hwKernel; l < rPixY - hwKernel; l++) {
                for (k = hwKernel; k < rPixX - hwKernel; k++) {
                    if (mRData[k + rPixX * l] & FLAG_OUTPUT_ISBAD) {
                        oRData[k + rPixX * l] = fillVal;
                        iRData[k + rPixX * l] = fillValNoise;
                    }
                }
            }
        }
    }

    /* difference image is the first output layer */
    if (fits_movabs_hdu(oPtr, 1, NULL, &status) ||
        hp_fits_write_subset(oPtr, 0, 2, oNaxes, oRData, &status,
                             outShort, outBzero, outBscale, fpixelOutX, fpixelOutY,
                             lpixelOutX, lpixelOutY, xBufLo, yBufLo, mRData, rPixX, rPixY))
        printError(status);

    /* extraneous code, but the convolved image is the second layer if its asked for... */
    if (inclConvImage)
        if (fits_movrel_hdu(oPtr, 1, NULL, &status))
            printError(status);

    /* sigma image is next in the queue, oRData already masked */
    if (inclSigmaImage || sigmaImage) {
        for (l = hwKernel; l < rPixY - hwKernel; l++) {
            for (k = hwKernel; k < rPixX - hwKernel; k++) {
                if (oRData[k + rPixX * l] != fillVal) {
                    if (convTmpl) {
                        oRData[k + rPixX * l] /= tRData[k + rPixX * l];
                    } else {
                        oRData[k + rPixX * l] /= iRData[k + rPixX * l];
                    }
                }
            }
        }

        if (inclSigmaImage) {
            if (fits_movrel_hdu(oPtr, 1, NULL, &status) ||
                hp_fits_write_subset(oPtr, 0, 2, oNaxes, oRData, &status,
                                     outShort, outBzero, outBscale, fpixelOutX, fpixelOutY,
                                     lpixelOutX, lpixelOutY, xBufLo, yBufLo, mRData, rPixX, rPixY))
                printError(status);
        }
        if (sigmaImage) {
            if (fits_open_file(&ePtr, sigmaImage, 1, &status))
                printError(status);
            if (hp_fits_write_subset(ePtr, 0, 2, oNaxes, oRData, &status,
                                     outShort, outBzero, outBscale, fpixelOutX, fpixelOutY,
                                     lpixelOutX, lpixelOutY, xBufLo, yBufLo, mRData, rPixX, rPixY) ||
                fits_close_file(ePtr, &status))
                printError(status);
        }
    }

    /* noise image is next in the queue */
    sprintf(hKeyword, "NSCALO%02d", i);
    sprintf(hInfo, "%.4f", diffrat);

    if (inclNoiseImage) {
        if (fits_movrel_hdu(oPtr, 1, NULL, &status) ||
            fits_set_bscale(oPtr, outNiBscale, outNiBzero, &status))
            printError(status);

        if (convTmpl) {
            if (hp_fits_write_subset(oPtr, 0, 2, oNaxes, tRData, &status,
                                     outNShort, outNiBzero, outNiBscale, fpixelOutX, fpixelOutY,
                                     lpixelOutX, lpixelOutY, xBufLo, yBufLo, mRData, rPixX, rPixY))
                printError(status);
        } else {
            if (hp_fits_write_subset(oPtr, 0, 2, oNaxes, iRData, &status,
                                     outNShort, outNiBzero, outNiBscale, fpixelOutX, fpixelOutY,
                                     lpixelOutX, lpixelOutY, xBufLo, yBufLo, mRData, rPixX, rPixY))
                printError(status);
        }

        /* reset output stream */
        if (fits_write_key_str(oPtr, hKeyword, hInfo, "", &status) ||
            fits_set_bscale(oPtr, outBscale, outBzero, &status))
            printError(status);
    }
    if (noiseImage) {
        if (fits_open_file(&ePtr, noiseImage, 1, &status))
            printError(status);
        if (fits_set_bscale(ePtr, outNiBscale, outNiBzero, &status))
            printError(status);
        if (convTmpl) {
            if (hp_fits_write_subset(ePtr, 0, 2, oNaxes, tRData, &status,
                                     outNShort, outNiBzero, outNiBscale, fpixelOutX, fpixelOutY,
                                     lpixelOutX, lpixelOutY, xBufLo, yBufLo, mRData, rPixX, rPixY))
                printError(status);
        } else {
            if (hp_fits_write_subset(ePtr, 0, 2, oNaxes, iRData, &status,
                                     outNShort, outNiBzero, outNiBscale, fpixelOutX, fpixelOutY,
                                     lpixelOutX, lpixelOutY, xBufLo, yBufLo, mRData, rPixX, rPixY))
                printError(status);
        }
        if (fits_write_key_str(ePtr, hKeyword, hInfo, "", &status) ||
            fits_close_file(ePtr, &status))
            printError(status);
    }

    if (outMask) {
        if (fits_open_file(&ePtr, outMask, 1, &status))
            printError(status);
        if (hp_fits_write_subset_int(ePtr, 0, 2, oNaxes, mRData, &status,
                                     1, 32768, 1, fpixelOutX, fpixelOutY,
                                     lpixelOutX, lpixelOutY, xBufLo, yBufLo, mRData, rPixX, rPixY) ||
            fits_close_file(ePtr, &status))
            printError(status);
    }

    /* add fits header info */
    fits_movabs_hdu(oPtr, 1, NULL, &status);

    if (kernelImOut)
        if (fits_open_file(&ePtr, kernelImOut, 1, &status))
            printError(status);

    /* add comment to TITLE saying which image was actually convolved */
    sprintf(hKeyword, "REGION%02d", i);
    sprintf(hInfo, "[%d:%d,%d:%d]", rXMin + 1, rXMax + 1, rYMin + 1, rYMax + 1);
    fits_write_key_str(oPtr, hKeyword, hInfo, "", &status);
    if (kernelImOut)
        fits_write_key_str(ePtr, hKeyword, hInfo, "", &status);

    sprintf(hKeyword, "CONVOL%02d", i);
    sprintf(hInfo, "%s", (convTmpl ? "TEMPLATE" : "IMAGE"));
    fits_write_key_str(oPtr, hKeyword, hInfo, "", &status);
    if (kernelImOut)
        fits_write_key_str(ePtr, hKeyword, hInfo, "", &status);

    sprintf(hKeyword, "KSUM%02d", i);
    sprintf(hInfo, "%.4f", sumKernel);
    fits_write_key_str(oPtr, hKeyword, hInfo, "Kernel Sum", &status);
    if (kernelImOut)
        fits_write_key_str(ePtr, hKeyword, hInfo, "Kernel Sum", &status);

    sprintf(hKeyword, "SSSIG%02d", i);
    sprintf(hInfo, "%.4f", meansigSubstamps);
    fits_write_key_str(oPtr, hKeyword, hInfo, "Average Figure of Merit across Stamps", &status);
    if (kernelImOut)
        fits_write_key_str(ePtr, hKeyword, hInfo, "Average Figure of Merit across Stamps", &status);

    sprintf(hKeyword, "SSSCAT%02d", i);
    sprintf(hInfo, "%.4f", scatterSubstamps);
    fits_write_key_str(oPtr, hKeyword, hInfo, "Stdev in Figure of Merit", &status);
    if (kernelImOut)
        fits_write_key_str(ePtr, hKeyword, hInfo, "Stdev in Figure of Merit", &status);


    if (strncmp(figMerit, "v", 1) == 0) {
        sprintf(hKeyword, "FSIG%02d", i);
        sprintf(hInfo, "%.4f", meansigSubstampsF);
        fits_write_key_str(oPtr, hKeyword, hInfo, "Final SSSIG", &status);
        if (kernelImOut)
            fits_write_key_str(ePtr, hKeyword, hInfo, "Final SSSIG", &status);

        sprintf(hKeyword, "FSCAT%02d", i);
        sprintf(hInfo, "%.4f", scatterSubstampsF);
        fits_write_key_str(oPtr, hKeyword, hInfo, "Final SSSCAT", &status);
        if (kernelImOut)
            fits_write_key_str(ePtr, hKeyword, hInfo, "Final SSSCAT", &status);

    }

    sprintf(hKeyword, "X2NRM%02d", i);
    sprintf(hInfo, "%.4f", x2norm);
    fits_write_key_str(oPtr, hKeyword, hInfo, "1/N * SUM (diff/noise)^2", &status);
    if (kernelImOut)
        fits_write_key_str(ePtr, hKeyword, hInfo, "1/N * SUM (diff/noise)^2", &status);

    sprintf(hKeyword, "NX2NRM%02d", i);
    sprintf(hInfo, "%d", nx2norm);
    fits_write_key_str(oPtr, hKeyword, hInfo, "Number of pixels in X2NRM", &status);
    if (kernelImOut)
        fits_write_key_str(ePtr, hKeyword, hInfo, "Number of pixels in X2NRM", &status);

    sprintf(hKeyword, "DMEAN%02d", i);
    sprintf(hInfo, "%.4f", mean);
    fits_write_key_str(oPtr, hKeyword, hInfo, "Mean of diff image; good pixels", &status);
    if (kernelImOut)
        fits_write_key_str(ePtr, hKeyword, hInfo, "Mean of diff image; good pixels", &status);

    sprintf(hKeyword, "DSIGE%02d", i);
    sprintf(hInfo, "%.4f", sd);
    fits_write_key_str(oPtr, hKeyword, hInfo, "Stdev of diff image; good pixels", &status);
    if (kernelImOut)
        fits_write_key_str(ePtr, hKeyword, hInfo, "Stdev of diff image; good pixels", &status);

    sprintf(hKeyword, "DSIG%02d", i);
    sprintf(hInfo, "%.4f", nmean);
    fits_write_key_str(oPtr, hKeyword, hInfo, "Mean of noise image; good pixels", &status);
    if (kernelImOut)
        fits_write_key_str(ePtr, hKeyword, hInfo, "Mean of noise image; good pixels", &status);

    sprintf(hKeyword, "DMEANO%02d", i);
    sprintf(hInfo, "%.4f", meanm);
    fits_write_key_str(oPtr, hKeyword, hInfo, "Mean of diff image; OK pixels", &status);
    if (kernelImOut)
        fits_write_key_str(ePtr, hKeyword, hInfo, "Mean of diff image; OK pixels", &status);

    sprintf(hKeyword, "DSIGEO%02d", i);
    sprintf(hInfo, "%.4f", sdm);
    fits_write_key_str(oPtr, hKeyword, hInfo, "Stdev of diff image; OK pixels", &status);
    if (kernelImOut)
        fits_write_key_str(ePtr, hKeyword, hInfo, "Stdev of diff image; OK pixels", &status);

    sprintf(hKeyword, "DSIGO%02d", i);
    sprintf(hInfo, "%.4f", nmeanm);
    fits_write_key_str(oPtr, hKeyword, hInfo, "Mean of noise image; OK pixels", &status);
    if (kernelImOut)
        fits_write_key_str(ePtr, hKeyword, hInfo, "Mean of noise image; OK pixels", &status);

    if (rescaleOK) {
        sprintf(hKeyword, "NSCALO%02d", i);
        sprintf(hInfo, "%.4f", diffrat);
        fits_write_key_str(oPtr, hKeyword, hInfo, "Rescaling factor of OK noise pixels", &status);
        if (kernelImOut)
            fits_write_key_str(ePtr, hKeyword, hInfo, "Rescaling factor of OK noise pixels", &status);
    }

    if (kernelImOut)
        fits_close_file(ePtr, &status);

    if (doKerInfo) {
        if (convTmpl) {
            if (fits_movabs_hdu(oPtr, kInfoNum, NULL, &status) ||
                fits_write_col_dbl(oPtr, i + 1, 1, 1, (nCompTotal + 1), tKerSol, &status))
                printError(status);
        } else {
            if (fits_movabs_hdu(oPtr, kInfoNum, NULL, &status) ||
                fits_write_col_dbl(oPtr, i + 1, 1, 1, (nCompTotal + 1), iKerSol, &status))
                printError(status);
        }
    }
    if (kernelImOut) {
        if (convTmpl) {
            if (fits_open_file(&ePtr, kernelImOut, 1, &status) ||
                fits_movabs_hdu(ePtr, kInfoNum, NULL, &status) ||
                fits_write_col_dbl(ePtr, i + 1, 1, 1, (nCompTotal + 1), tKerSol, &status) ||
                fits_close_file(ePtr, &status))
                printError(status);
        } else {
            if (fits_open_file(&ePtr, kernelImOut, 1, &status) ||
                fits_movabs_hdu(ePtr, kInfoNum, NULL, &status) ||
                fits_write_col_dbl(ePtr, i + 1, 1, 1, (nCompTotal + 1), iKerSol, &status) ||
                fits_close_file(ePtr, &status))
                printError(status);
        }
    }

    /* close open fits files */
    if (fits_close_file(oPtr, &status)) {
        printError(status);
    }

    /* release the lock around all write operations to prevent multi-thread access to fits files */
    fprintf(stderr,"[Region %d] Releasing writing lock\n", i);
    pthread_mutex_unlock(&writelock);

    /* free things now we're done with them */
    free(tRData);
    tRData = NULL;
    free(iRData);
    iRData = NULL;
    free(oRData);
    oRData = NULL;
    free(eRData);
    eRData = NULL;
    free(mRData);
    mRData = NULL;
    free(misRData);
    misRData = NULL;
    free(mtsRData);
    mtsRData = NULL;

    if (ctStamps) {
        free(ctStamps);
        ctStamps = NULL;
    }
    if (ciStamps) {
        free(ciStamps);
        ciStamps = NULL;
    }
    if (iKerSol) {
        free(iKerSol);
        iKerSol = NULL;
    }
    if (tKerSol) {
        free(tKerSol);
        tKerSol = NULL;
    }

    if (temp) free(temp);
    if (indx) free(indx);
    if (kernel) free(kernel);
    if (kernel_coeffs) free(kernel_coeffs);
    for (j = 0; j < nCompKer; j++)
        if (kernel_vec[j]) free(kernel_vec[j]);
    if (kernel_vec) free(kernel_vec);
    if (filter_x) free(filter_x);
    if (filter_y) free(filter_y);
    if (check_vec) free(check_vec);
    for (j = 0; j < nC; j++)
        if (check_mat[j]) free(check_mat[j]);
    if (check_mat) free(check_mat);
    if (check_stack) free(check_stack);
    fprintf(stderr, "[Region %i] Finished\n\n", i);

}



int main(int argc, char *argv[]) {
    int i, j, k;                              /* generic indices */
    char scrStr[SCRLEN];                         /* scratch string */

    /*
      NOTE : The .fits images are 1-indexed, while these arrays contained within are
      0 indexed.  Thus anything you find below at x,y you will find displayed at
      position x+1,y+1.
    */

    /* template image */
    fitsfile *tPtr;
    int tBitpix, tNaxis;
    long tNaxes[MAXDIM];

    /* comparison image */
    fitsfile *iPtr;
    int iBitpix, iNaxis;
    long iNaxes[MAXDIM];

    /* output image */
    fitsfile *oPtr;
    int oBitpix, oNaxis;
    long oNaxes[MAXDIM];

    /* extraneous image; eRData now defined in globals.h */
    fitsfile *ePtr;

    /* extraneous output regions */

    int xMin, yMin, xMax, yMax;             /* whole image */

    int *rXMins, *rYMins, *rXMaxs, *rYMaxs; /* each region, good data */
    int rXMin, rYMin, rXMax, rYMax;         /* each region, good data */
    int nR;                                 /* region counter */

    /* cfitsio */
    int status = 0, kInfoNum;
    char **tform = NULL, **ttype = NULL, **tunit = NULL;
    char hKeyword[1024], hInfo[2048];

    /* misc */
    FILE *rFile;
    char *pstr;

    struct tm *tm;
    time_t thetime;

    /* SET VERSION */
    sprintf(version, "5.1.11");

    /* set global vars, and grab command line args */
    /*   return image names in fnames */
    vargs(argc, argv);

    /******/
    /* access input images, create output image shell */
    /******/

    /* open up, get template bitpix, # dimensions, image size */
    if (fits_open_file(&tPtr, template, 0, &status))
        printError(status);
    if (fits_get_img_param(tPtr, MAXDIM, &tBitpix, &tNaxis, tNaxes, &status))
        printError(status);

    /* if input noise image, open and check to make sure its the right size */
    /* use the ePtr fitsfile temporarily */
    if (tNoiseIm) {
        if (fits_open_file(&ePtr, tNoiseIm, 0, &status))
            printError(status);
        if (fits_get_img_param(ePtr, MAXDIM, &oBitpix, &oNaxis, oNaxes, &status))
            printError(status);

        if ((oNaxes[0] != tNaxes[0]) || (oNaxes[1] != tNaxes[1])) {
            fprintf(stderr, "WARNING : Input template noise array not same size as template, ignoring...\n");
            tNoiseIm = NULL;
        }

        if (fits_close_file(ePtr, &status))
            printError(status);
    }
    /* open up, get comparison image bitpix, # dimensions, image size */
    if (fits_open_file(&iPtr, image, 0, &status))
        printError(status);
    if (fits_get_img_param(iPtr, MAXDIM, &iBitpix, &iNaxis, iNaxes, &status))
        printError(status);

    /* if input noise image, open and check to make sure its the right size */
    /* use the ePtr fitsfile temporarily */
    if (iNoiseIm) {
        if (fits_open_file(&ePtr, iNoiseIm, 0, &status))
            printError(status);
        if (fits_get_img_param(ePtr, MAXDIM, &oBitpix, &oNaxis, oNaxes, &status))
            printError(status);

        if ((oNaxes[0] != iNaxes[0]) || (oNaxes[1] != iNaxes[1])) {
            fprintf(stderr, "WARNING : Input image noise array not same size as image, ignoring...\n");
            iNoiseIm = NULL;
        }
        if (fits_close_file(ePtr, &status))
            printError(status);
    }
    if (tNoiseIm) {
        if (fits_open_file(&ePtr, tNoiseIm, 0, &status))
            printError(status);
        if (fits_get_img_param(ePtr, MAXDIM, &oBitpix, &oNaxis, oNaxes, &status))
            printError(status);

        if ((oNaxes[0] != iNaxes[0]) || (oNaxes[1] != iNaxes[1])) {
            fprintf(stderr, "WARNING : Input template noise array not same size as image, ignoring...\n");
            tNoiseIm = NULL;
        }
        if (fits_close_file(ePtr, &status))
            printError(status);
    }

    /* same with mask images */
    if (iMaskIm) {
        if (fits_open_file(&ePtr, iMaskIm, 0, &status))
            printError(status);
        if (fits_get_img_param(ePtr, MAXDIM, &oBitpix, &oNaxis, oNaxes, &status))
            printError(status);

        if ((oNaxes[0] != iNaxes[0]) || (oNaxes[1] != iNaxes[1])) {
            fprintf(stderr, "WARNING : Input image mask array not same size as image, ignoring...\n");
            iMaskIm = NULL;
        }
        if (fits_close_file(ePtr, &status))
            printError(status);
    }
    if (tMaskIm) {
        if (fits_open_file(&ePtr, tMaskIm, 0, &status))
            printError(status);
        if (fits_get_img_param(ePtr, MAXDIM, &oBitpix, &oNaxis, oNaxes, &status))
            printError(status);

        if ((oNaxes[0] != iNaxes[0]) || (oNaxes[1] != iNaxes[1])) {
            fprintf(stderr, "WARNING : Input template mask array not same size as image, ignoring...\n");
            tMaskIm = NULL;
        }
        if (fits_close_file(ePtr, &status))
            printError(status);
    }

    /* let em know whats going on... */
    fprintf(stderr, "Doing : %s -\n", image);
    fprintf(stderr, "        %s =\n", template);
    fprintf(stderr, "        %s\n", outim);
    fprintf(stderr, "   Good templ data : %.1f -> %.1f\n", tLThresh, tUThresh);
    fprintf(stderr, "   Good image data : %.1f -> %.1f\n", iLThresh, iUThresh);

    /* ADU pedestal? */
    tUThresh -= tPedestal;
    tUKThresh -= tPedestal;
    tLThresh -= tPedestal;
    iUThresh -= iPedestal;
    iUKThresh -= iPedestal;
    iLThresh -= iPedestal;

    /* overwrite output image? */
    if (!noClobber) { sprintf(scrStr, "!%s", outim); }
    else { sprintf(scrStr, "%s", outim); }

    oNaxes[0] = imax(tNaxes[0], iNaxes[0]);
    oNaxes[1] = imax(tNaxes[1], iNaxes[1]);
    oNaxis = 2;

    if (outShort)
        oBitpix = SHORT_IMG;
    else
        oBitpix = FLOAT_IMG;

    /* create and open new empty output FITS file */
    if (fits_create_file(&oPtr, scrStr, &status) ||
        fits_create_img(oPtr, oBitpix, oNaxis, oNaxes, &status))
        printError(status);

    /* copy over the fits keywords */
    if (hp_fits_copy_header(iPtr, oPtr, &status))
        printError(status);

    /* add info; format taken from SWARP */
    fits_write_key_str(oPtr, "COMMENT", "", "", &status);
    fits_write_key_str(oPtr, "SOFTNAME", "HOTPanTS", "The software that differenced this image", &status);
    fits_write_key_str(oPtr, "SOFTVERS", version, "Version", &status);
    fits_write_key_str(oPtr, "SOFTAUTH", "Andrew Becker <becker@astro.washington.edu>", "Maintainer", &status);
    fits_write_key_str(oPtr, "SOFTINST", "University of Washington", "Institute", &status);
    fits_write_key_str(oPtr, "COMMENT", "", "", &status);

    /* user name */
    if (!(pstr = getenv("USERNAME")))       /* Cygwin,... */
        if ((pstr = getenv("LOGNAME")))       /* Linux,... */
            fits_write_key_str(oPtr, "AUTHOR", pstr, "Who ran the software", &status);
    if (!pstr)
        fits_write_key_str(oPtr, "AUTHOR", "unknown", "Who ran the software", &status);

    /* host name */
    if (!gethostname(hInfo, 80))
        fits_write_key_str(oPtr, "ORIGIN", hInfo, "Where it was done", &status);

    /* what was it fired up */
    thetime = time(NULL);
    tm = gmtime(&thetime);
    sprintf(hInfo, "%04d-%02d-%02dT%02d:%02d:%02d",
            tm->tm_year + 1900, tm->tm_mon + 1, tm->tm_mday,
            tm->tm_hour, tm->tm_min, tm->tm_sec);
    fits_write_key_str(oPtr, "DATE", hInfo, "When it was started (GMT)", &status);

    fits_write_key_str(oPtr, "COMMENT", "", "", &status);
    /* don't let a failure here kill us */
    status = 0;

    /* change bscale/bzero AFTER transfering input image header */
    if (fits_update_key_flt(oPtr, "BZERO", outBzero, -5, "", &status) ||
        fits_update_key_flt(oPtr, "BSCALE", outBscale, -5, "", &status))
        printError(status);

    /* make new image extensions for any extraneous output data */
    if (inclConvImage) {
        if (fits_insert_img(oPtr, oBitpix, oNaxis, oNaxes, &status) ||
            fits_update_key(oPtr, TSTRING, "OBJECT", "Convolved Image", "", &status))
            printError(status);
    }
    if (convImage) {
        if (!noClobber) { sprintf(scrStr, "!%s", convImage); }
        else { sprintf(scrStr, "%s", convImage); }
        if (fits_create_file(&ePtr, scrStr, &status) ||
            fits_create_img(ePtr, oBitpix, oNaxis, oNaxes, &status) ||
            hp_fits_copy_header(iPtr, ePtr, &status) ||
            fits_update_key(ePtr, TSTRING, "OBJECT", "Convolved Image", "", &status) ||
            fits_close_file(ePtr, &status))
            printError(status);
    }

    if (inclSigmaImage) {
        if (fits_insert_img(oPtr, oBitpix, 2, oNaxes, &status) ||
            fits_update_key(oPtr, TSTRING, "OBJECT", "Noise-scaled (sigma) Difference Image", "", &status))
            printError(status);
    }
    if (sigmaImage) {
        if (!noClobber) { sprintf(scrStr, "!%s", sigmaImage); }
        else { sprintf(scrStr, "%s", sigmaImage); }
        if (fits_create_file(&ePtr, scrStr, &status) ||
            fits_create_img(ePtr, oBitpix, oNaxis, oNaxes, &status) ||
            hp_fits_copy_header(iPtr, ePtr, &status) ||
            fits_update_key(ePtr, TSTRING, "OBJECT", "Noise-scaled (sigma) Difference Image", "", &status) ||
            fits_close_file(ePtr, &status))
            printError(status);
    }

    if (inclNoiseImage) {
        if (outNShort) {
            if (fits_insert_img(oPtr, SHORT_IMG, 2, oNaxes, &status))
                printError(status);
        } else {
            if (fits_insert_img(oPtr, FLOAT_IMG, 2, oNaxes, &status))
                printError(status);
        }
        if (fits_update_key(oPtr, TSTRING, "OBJECT", "HOTPanTS Noise Image", "", &status))
            printError(status);

        /* manually set bscale and bzero for noise image layer */
        if (outNShort) {
            if (fits_update_key_flt(oPtr, "BZERO", outNiBzero, -5, "", &status) ||
                fits_update_key_flt(oPtr, "BSCALE", outNiBscale, -5, "", &status))
                printError(status);
        }
    }
    if (noiseImage) {
        if (!noClobber) { sprintf(scrStr, "!%s", noiseImage); }
        else { sprintf(scrStr, "%s", noiseImage); }

        if (fits_create_file(&ePtr, scrStr, &status))
            printError(status);

        if (outNShort) {
            if (fits_create_img(ePtr, SHORT_IMG, oNaxis, oNaxes, &status))
                printError(status);
        } else {
            if (fits_create_img(ePtr, FLOAT_IMG, oNaxis, oNaxes, &status))
                printError(status);
        }

        if (hp_fits_copy_header(iPtr, ePtr, &status) ||
            fits_update_key(ePtr, TSTRING, "OBJECT", "HOTPanTS Noise Image", "", &status) ||
            fits_update_key_flt(ePtr, "GAIN", 1., -1, "No gain in noise image", &status) ||
            fits_update_key_flt(ePtr, "RDNOISE", 0., -1, "No rdnoise in noise image", &status) ||
            fits_write_key_flt(ePtr, "MASKVAL", fillValNoise, -6, "Value of Masked Pixels", &status))
            printError(status);

        if (outNShort)
            if (fits_update_key_flt(ePtr, "BZERO", outNiBzero, -5, "", &status) ||
                fits_update_key_flt(ePtr, "BSCALE", outNiBscale, -5, "", &status))
                printError(status);

        if (fits_close_file(ePtr, &status))
            printError(status);
    }

    if (outMask) {
        if (!noClobber) { sprintf(scrStr, "!%s", outMask); }
        else { sprintf(scrStr, "%s", outMask); }
        if (fits_create_file(&ePtr, scrStr, &status) ||
            fits_create_img(ePtr, SHORT_IMG, oNaxis, oNaxes, &status) ||
            hp_fits_copy_header(iPtr, ePtr, &status) ||
            fits_update_key(ePtr, TSTRING, "OBJECT", "Hotpants output Mask Image", "", &status) ||
            fits_update_key_flt(ePtr, "BZERO", 32768, -5, "", &status) ||
            fits_update_key_flt(ePtr, "BSCALE", 1, -5, "", &status) ||
            fits_close_file(ePtr, &status))
            printError(status);
    }

    if (kernelImIn) {
        /* possibly override defaults with info in kernel image */
        getKernelInfo(kernelImIn);
        fprintf(stderr, "   received kernel info\n");
    }

    /* determine the size of the data structures */
    nCompKer = 0;
    for (i = 0; i < ngauss; i++)
        nCompKer += ((deg_fixe[i] + 1) * (deg_fixe[i] + 2)) / 2;

    nComp = ((kerOrder + 1) * (kerOrder + 2)) / 2;
    nC = nCompKer + 2;
    nBGVectors = ((bgOrder + 1) * (bgOrder + 2)) / 2;
    nCompTotal = nCompKer * nComp + nBGVectors;

    fwKernel = hwKernel * 2 + 1;         /* kernel size */

    if (useFullSS) {
        fwKSStamp = fwKernel;     /* substamp size */
        fwStamp = fwKernel;
        nStampX = (int) (imin(tNaxes[0], iNaxes[0]) / nRegX / fwStamp);
        nStampY = (int) (imin(tNaxes[1], iNaxes[1]) / nRegY / fwStamp);
        fprintf(stderr, "Using maximial number of stamps : %d x %d\n", nStampX, nStampY);
    } else {
        fwKSStamp = hwKSStamp * 2 + 1;     /* substamp size */

        /* estimate of fwStamp for all these mallocs... */
        /* smallest dimension, minus kernel width */
        fwStamp = imin(imin(tNaxes[0], iNaxes[0]) / nRegX / nStampX,
                       imin(tNaxes[1], iNaxes[1]) / nRegY / nStampY);
        fwStamp -= fwKernel;
        fwStamp -= fwStamp % 2 == 0 ? 1 : 0; /* hmmm, an odd shape... */

        /* insanity checking */
        if (fwStamp < fwKSStamp) {
            fwStamp = fwKSStamp + fwKernel;
            fwStamp -= fwStamp % 2 == 0 ? 1 : 0;

            nStampX = imin(tNaxes[0], iNaxes[0]) / nRegX / fwStamp;
            nStampY = imin(tNaxes[1], iNaxes[1]) / nRegY / fwStamp;

            fprintf(stderr, "WARNING : too many stamps requested\n");
            fprintf(stderr, "          using nsx = %d, nsy = %d\n", nStampX, nStampY);
        }

    }

    kcStep = kcStep ? kcStep : fwKernel; /* size of step in spatial_convolve */

    nStamps = nStampX * nStampY;


    /******/
    /* determine pixel limits of regions, and of stamps in region */
    /******/

    xMin = 0;
    yMin = 0;
    xMax = -1 + imin(tNaxes[0], iNaxes[0]);
    yMax = -1 + imin(tNaxes[1], iNaxes[1]);

    nR = 0;
    if (regFile) {
        rFile = fopen(regFile, "r");
        while (fscanf(rFile, "%s\n", scrStr) != EOF)
            nR++;
        rewind(rFile);

        if (!(rXMins = (int *) calloc(nR, sizeof(int))) ||
            !(rXMaxs = (int *) calloc(nR, sizeof(int))) ||
            !(rYMins = (int *) calloc(nR, sizeof(int))) ||
            !(rYMaxs = (int *) calloc(nR, sizeof(int))))
            exit(1);

        nR = 0;
        while (fscanf(rFile, "%d:%d,%d:%d", &rXMin, &rXMax, &rYMin, &rYMax) != EOF) {
            /* range of good data for the region */
            rXMins[nR] = imax(rXMin, xMin);
            rYMins[nR] = imax(rYMin, yMin);
            rXMaxs[nR] = imin(rXMax, xMax);
            rYMaxs[nR] = imin(rYMax, yMax);
            nR++;
        }
        fclose(rFile);
    } else if (regKeyWord) {
        if (!(rXMins = (int *) calloc(numRegKeyWord, sizeof(int))) ||
            !(rXMaxs = (int *) calloc(numRegKeyWord, sizeof(int))) ||
            !(rYMins = (int *) calloc(numRegKeyWord, sizeof(int))) ||
            !(rYMaxs = (int *) calloc(numRegKeyWord, sizeof(int))))
            exit(1);

        for (nR = 0; nR < numRegKeyWord; nR++) {
            sprintf(hKeyword, "%s%d", regKeyWord, nR);

            if (fits_read_key_str(tPtr, hKeyword, hInfo, scrStr, &status)) {
                rXMins[nR] = xMin;
                rXMaxs[nR] = xMax;
                rYMins[nR] = yMin;
                rYMaxs[nR] = yMax;

                status = 0;
            } else
                sscanf(hInfo, "[%d:%d,%d:%d]", &rXMins[nR], &rXMaxs[nR], &rYMins[nR], &rYMaxs[nR]);

            if (fits_read_key_str(iPtr, hKeyword, hInfo, scrStr, &status)) {
                rXMin = xMin;
                rXMax = xMax;
                rYMin = yMin;
                rYMax = yMax;

                status = 0;
            } else
                sscanf(hInfo, "[%d:%d,%d:%d]", &rXMin, &rXMax, &rYMin, &rYMax);

            rXMins[nR] = imax(rXMin, rXMins[nR]);
            rYMins[nR] = imax(rYMin, rYMins[nR]);
            rXMaxs[nR] = imin(rXMax, rXMaxs[nR]);
            rYMaxs[nR] = imin(rYMax, rYMaxs[nR]);

            rXMins[nR] = imax(xMin, rXMins[nR]);
            rYMins[nR] = imax(yMin, rYMins[nR]);
            rXMaxs[nR] = imin(xMax, rXMaxs[nR]);
            rYMaxs[nR] = imin(yMax, rYMaxs[nR]);

            /* fprintf(stderr, "%d %d %d %d\n", rXMins[nR], rYMins[nR], rXMaxs[nR], rYMaxs[nR]); */

        }
    } else {
        if (!(rXMins = (int *) calloc(nRegX * nRegY, sizeof(int))) ||
            !(rXMaxs = (int *) calloc(nRegX * nRegY, sizeof(int))) ||
            !(rYMins = (int *) calloc(nRegX * nRegY, sizeof(int))) ||
            !(rYMaxs = (int *) calloc(nRegX * nRegY, sizeof(int))))
            exit(1);

        /* in cfitsio, data are striped along x dimen, thus all loops
           are y outer, x inner.  or at least they should be... */

        for (j = 0; j < nRegY; j++) {
            for (i = 0; i < nRegX; i++) {
                /* range of good data for the region */
                rXMins[nR] = xMin + i * xMax / nRegX;
                rYMins[nR] = yMin + j * yMax / nRegY;
                rXMaxs[nR] = imin((i + 1) * xMax / nRegX, xMax);
                rYMaxs[nR] = imin((j + 1) * yMax / nRegY, yMax);
                nR++;
            }
        }
    }

    /* now that we know the number of regions, if we want a fits
       binary table for the kernel info, create it here */

    if ((doKerInfo) || (kernelImOut)) {
        tform = (char **) calloc(nR, sizeof(char *));
        ttype = (char **) calloc(nR, sizeof(char *));
        tunit = (char **) calloc(nR, sizeof(char *));
        for (k = 0; k < nR; k++) {
            tform[k] = (char *) malloc(10 * sizeof(char));
            ttype[k] = (char *) malloc(10 * sizeof(char));
            tunit[k] = (char *) malloc(sizeof(char));
            strcpy(tform[k], "D10.8");
            sprintf(ttype[k], "Region%d", k);
            strcpy(tunit[k], "");
        }

        if (doKerInfo) {
            /* lets try it as a layer at the **bottom** of the output difference image */
            if (fits_insert_btbl(oPtr, (nCompTotal + 1), nR, ttype, tform, tunit, "Convolution Kernel Information", 0L,
                                 &status) ||
                fits_update_key(oPtr, TSTRING, "OBJECT", "Convolution Kernel Information", "", &status) ||
                fits_get_hdu_num(oPtr, &kInfoNum))
                printError(status);
        }
        if (kernelImOut) {
            if (!noClobber) { sprintf(scrStr, "!%s", kernelImOut); }
            else { sprintf(scrStr, "%s", kernelImOut); }
            if (fits_create_file(&ePtr, scrStr, &status) ||
                fits_insert_btbl(ePtr, (nCompTotal + 1), nR, ttype, tform, tunit, "Convolution Kernel Information", 0L,
                                 &status) ||
                fits_update_key(ePtr, TSTRING, "OBJECT", "Convolution Kernel Information", "", &status) ||
                fits_get_hdu_num(ePtr, &kInfoNum) ||
                fits_close_file(ePtr, &status))
                printError(status);
        }
    }

    xcmp = ycmp = 0;
    if (sstampFile) {
        /* downsize for speed */
        /* armin: this screws things up! */
        /*fwStamp = fwKSStamp + fwKernel;*/

        /* armin: call function to load x,y from file */
        loadxyfile(sstampFile, cmpFile);
    }

    /* limit number of threads to number of regions */
    if (nThread > nR) {
        fprintf(stderr, "Limiting number of threads to number of regions (%d)\n", nR);
        nThread = nR;
    }

    /* set up threads */
//    FIXME allow for nThread < nR. For now add this:
    nThread = nR;

    pthread_t thread[nThread];

    struct RegionData regData[nThread];
//    int regPerThread=(nR+nThread-1)/nThread;
    for (i=0; i<nThread; i++) {
//        regData[i].regNumStart=i*regPerThread;
//        regData[i].regNumStop=(i+1)*regPerThread;
        regData[i].regNum=i;
        regData[i].nR=nR;
        regData[i].xMin=xMin;
        regData[i].yMin=yMin;
        regData[i].xMax=xMax;
        regData[i].yMax=yMax;
        for (int ii = 0; ii < MAXDIM; ii++) {
            regData[i].tNaxes[ii] = tNaxes[ii];
            regData[i].iNaxes[ii] = iNaxes[ii];
            regData[i].oNaxes[ii] = oNaxes[ii];
        }
        regData[i].tNaxis=tNaxis;
        regData[i].iNaxis=iNaxis;
        regData[i].rXMin = rXMins[i];
        regData[i].rYMin = rYMins[i];
        regData[i].rXMax = rXMaxs[i];
        regData[i].rYMax = rYMaxs[i];
        regData[i].kInfoNum=kInfoNum;
        /* TODO What if nThread != nR - need to pass regPerThread to each thread? */
    }

    /* ensure last thread does not try to process a non-existant region */
//    regData[nThread - 1].stop=nR;
    /* start the threads and wait for them to complete */

    for (i=0; i<nThread; i++) {
        pthread_create(&thread[i], NULL, processRegion, &regData[i]);
    }

    for (i=0; i<nThread; i++) {
        pthread_join(thread[i], NULL);
    }


    /* add fits header info */
    fits_movabs_hdu(oPtr, 1, NULL, &status);

    /*fits_write_key_str(oPtr, "HPVSN", version, "HOTPANTS Version", &status);*/

    sprintf(hKeyword, "DIFFCMD");
    sprintf(hInfo, "%s", argv[0]);
    for (i = 1; i < argc; i++)
        sprintf(hInfo, "%s %s", hInfo, argv[i]);
    fits_write_key_longstr(oPtr, hKeyword, hInfo, "", &status);

    fits_write_key_longstr(oPtr, "TARGET", image, "HOTPanTS : Input Image", &status);
    fits_write_key_longstr(oPtr, "TEMPLATE", template, "HOTPanTS : Reference Image", &status);
    fits_write_key_longstr(oPtr, "DIFFIM", outim, "HOTPanTS : Difference Image", &status);

    fits_write_key_str(oPtr, "PHOTNORM", photNormalize, "Direction of photometric normalization", &status);
    fits_write_key_lng(oPtr, "NREGION", nRegX * nRegY, "Number of independent regions", &status);
    fits_write_key_flt(oPtr, "MASKVAL", fillVal, -6, "Value of Masked Pixels", &status);

    if (doKerInfo) {
        /* in primary HDU */
        fits_write_key_log(oPtr, "KERINFO", 1, "", &status);

        /* to the convolution kernel table */
        fits_movabs_hdu(oPtr, kInfoNum, NULL, &status);

        sprintf(hKeyword, "NGAUSS");
        fits_write_key_lng(oPtr, hKeyword, ngauss, "Number of Gaussian Basis Functions", &status);
        for (k = 0; k < ngauss; k++) {
            sprintf(hKeyword, "DGAUSS%d", k + 1);
            fits_write_key_lng(oPtr, hKeyword, deg_fixe[k], "Polynomial Degree", &status);

            sprintf(hKeyword, "SGAUSS%d", k + 1);
            fits_write_key_flt(oPtr, hKeyword, 1. / sqrt(2. * sigma_gauss[k]), 8, "Gaussian Sigma", &status);
        }

        fits_write_key_lng(oPtr, "FWKERN", fwKernel, "Kernel Width in Pixels", &status);
        fits_write_key_lng(oPtr, "CKORDER", kerOrder, "Spatial Convolution Order", &status);
        fits_write_key_lng(oPtr, "BGORDER", bgOrder, "Background Order", &status);
    } else {
        if (!(kernelImOut))
            fits_write_key_log(oPtr, "KERINFO", 0, "", &status);
    }

    if (kernelImOut) {
        /* in primary HDU */
        if (fits_open_file(&ePtr, kernelImOut, 1, &status))
            printError(status);

        /* in diffim file */
        fits_write_key_log(oPtr, "KERINFO", 0, "", &status);
        fits_write_key_str(oPtr, "KERFILE", kernelImOut, "", &status);
        /* in kernel */
        fits_write_key_log(ePtr, "KERINFO", 1, "", &status);

        /* to the convolution kernel table */
        fits_movabs_hdu(ePtr, kInfoNum, NULL, &status);

        fits_write_key_str(ePtr, "PHOTNORM", photNormalize, "Direction of photometric normalization", &status);

        sprintf(hKeyword, "NGAUSS");
        fits_write_key_lng(ePtr, hKeyword, ngauss, "Number of Gaussian Basis Functions", &status);
        fits_write_key_lng(oPtr, hKeyword, ngauss, "Number of Gaussian Basis Functions", &status);
        for (k = 0; k < ngauss; k++) {
            sprintf(hKeyword, "DGAUSS%d", k + 1);
            fits_write_key_lng(ePtr, hKeyword, deg_fixe[k], "Polynomial Degree", &status);
            fits_write_key_lng(oPtr, hKeyword, deg_fixe[k], "Polynomial Degree", &status);

            sprintf(hKeyword, "SGAUSS%d", k + 1);
            fits_write_key_flt(ePtr, hKeyword, 1. / sqrt(2. * sigma_gauss[k]), 8, "Gaussian Sigma", &status);
            fits_write_key_flt(oPtr, hKeyword, 1. / sqrt(2. * sigma_gauss[k]), 8, "Gaussian Sigma", &status);
        }

        fits_write_key_lng(ePtr, "FWKERN", fwKernel, "Kernel Width in Pixels", &status);
        fits_write_key_lng(ePtr, "CKORDER", kerOrder, "Spatial Convolution Order", &status);
        fits_write_key_lng(ePtr, "BGORDER", bgOrder, "Background Order", &status);
        fits_write_key_lng(oPtr, "FWKERN", fwKernel, "Kernel Width in Pixels", &status);
        fits_write_key_lng(oPtr, "CKORDER", kerOrder, "Spatial Convolution Order", &status);
        fits_write_key_lng(oPtr, "BGORDER", bgOrder, "Background Order", &status);
        fits_close_file(ePtr, &status);
    } else {
        fits_write_key_log(oPtr, "KERINFO", 0, "", &status);
    }

    /* close em down */
    if (fits_close_file(iPtr, &status) || fits_close_file(tPtr, &status) || fits_close_file(oPtr, &status)) {
        printError(status);
    }

    /* free anything left */
    free(rXMins);
    free(rXMaxs);
    free(rYMins);
    free(rYMaxs);

    free(deg_fixe);
    free(sigma_gauss);

    if (xcmp) free(xcmp);
    if (ycmp) free(ycmp);


    if ((doKerInfo) || (kernelImOut)) {
        for (k = 0; k < nR; k++) {
            if (tform[k]) free(tform[k]);
            if (ttype[k]) free(ttype[k]);
            if (tunit[k]) free(tunit[k]);
        }
        if (tform) free(tform);
        if (ttype) free(ttype);
        if (tunit) free(tunit);
    }

    fprintf(stderr, "SUCCESS\n");

    /* freedom itself! */
    return 0;
}
