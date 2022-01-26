#include <fitsio.h>

/* Alard.c */
void getKernelVec(double *, double *, double **);

int fillStamp(stamp_struct *, float *, float *, int *, int, int, double *, double *, float *);

double *kernel_vector(int, int, int, int, int *, double *, double *, double **);

void xy_conv_stamp(stamp_struct *, const float *, int, int, int, const double *, const double *, float *);

void fitKernel(stamp_struct *, float *, float *, float *, double *, double *, double *, int *, int *, int, int, int, double *, double *, int *, float *);

void build_matrix0(stamp_struct *);

void build_scprod0(stamp_struct *, const float *, int);

double check_stamps(stamp_struct *, int, float *, float *, int *, int, int, double **, double *, double *, double **, double *, double *, double **, int *, float *);

void build_matrix(stamp_struct *, int, double **, int, int, double **);

void build_scprod(stamp_struct *, int, const float *, double *, int, double **);

void getStampSig(stamp_struct *, double *, const float *, double *, double *, double *, int *, int, int, float *);

void getFinalStampSig(stamp_struct *, const float *, const float *, double *, const int *, int, int);

char check_again(stamp_struct *, double *, float *, float *, float *, double *, double *, int *, int *, int, int, int, double *, double *, float *);

void spatial_convolve(const float *, float **, int, int, double *, float *, const int *, int *, int, int, double *, double *, double **);

double make_kernel(int, int, const double *, int, int, double *, double *, double **);

double get_background(int, int, const double *, int, int);

void make_model(stamp_struct *, const double *, float *, int, int);

int ludcmp(double **, int, int *, double *);

void lubksb(double **, int, const int *, double *);

/* Functions.c */
int allocateStamps(stamp_struct *, int);

void buildStamps(int, int, int, int, int *, int *, int, int, int,
                 stamp_struct *, stamp_struct *, float *, float *,
                 float, float, int *, int, int);

void cutStamp(float *, float *, int, int, int, int, int, stamp_struct *);

int cutSStamp(stamp_struct *, float *, int *, int, int);

double checkPsfCenter(float *, int, int, int, int, int, int, double, float, float,
                      int, int, int, int, int *, int, int);

int getPsfCenters(stamp_struct *, float *, int, int, double, int, int, int *, int, int);

int getStampStats3(float *, int, int, int, int, double *, double *, double *, double *, double *, double *, double *, int,
               int, int, int *, int, int);

void getNoiseStats3(float *, float *, double *, int *, int, int, int *, int, int);

int sigma_clip(float *, int, double *, double *, int);

void freeStampMem(stamp_struct *, int);

void makeNoiseImage4(float *, float *, float, float, int, int);

void getKernelInfo(char *);

void readKernel(char *, int, double **, double **, int *, int *, int *, int *, double *, double *, double *, double *,
                double *, int *);

void fits_get_kernel_btbl(fitsfile *, double **, int);

void spreadMask(int *, int, int, int);

void makeInputMask(float *, float *, int *, int, int);

int hp_fits_copy_header(fitsfile *, fitsfile *, int *);

void hp_fits_correct_data(float *, int, float, float, int, int *);

void hp_fits_correct_data_int(int *, int, float, float, int, int *);

int hp_fits_write_subset(fitsfile *, long, long, long *,
                         float *, int *,
                         int, float, float,
                         int, int, int, int, int, int, int *, int, int);

int hp_fits_write_subset_int(fitsfile *, long, long, long *,
                             int *, int *,
                             int, float, float,
                             int, int, int, int, int, int, int *, int, int);

void fset(float *, double, int, int);

void dfset(double *, double, int, int);

void printError(int);

double ran1(int *);

void quick_sort(double *, int *, int);

int imin(int, int);

int imax(int, int);

/* armin */
void savexy(stamp_struct *, int, long, long, int);

void loadxyfile(char *, int);

/* Vargs.c */
void vargs(int, char *[]);
