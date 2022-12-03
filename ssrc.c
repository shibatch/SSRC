#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <time.h>

#include "dither.h"

#include "sleefdft.h"

#define VERSION "1.33"

#ifndef HIGHPREC
typedef float REAL;
double AA=140;
double DF=2000;
#define M 30
#else
typedef double REAL;
double AA=150;
double DF=200;
#define M 30
#endif

#define MAXNCH 10

#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795028842
#endif

#define RANDBUFLEN 65536

#define RINT(x) ((x) >= 0 ? ((int)((x) + 0.5)) : ((int)((x) - 0.5)))

double fact[M+1];

static inline void ignoreReturnValue(int x) {}

int quiet = 0;
int lastshowed2;
time_t starttime, lastshowed;

struct SSRCDither *shaper[MAXNCH];

int do_shaping(double s, double *peak, int dtype, int ch) {
  int ret;
  SSRCDither_quantizeDouble(shaper[ch], &ret, &s, 1, 1.0);
  return ret;
}

double alpha(double a)
{
  if (a <= 21) return 0;
  if (a <= 50) return 0.5842*pow(a-21,0.4)+0.07886*(a-21);
  return 0.1102*(a-8.7);
}

double izero(double x) {
  double ret = 1;
  int m;

  for(m=M;m>=1;m--) {
    double t = pow(x/2, m)/fact[m];
    ret += t*t;
  }

  return ret;
}

double win(double n,int len,double alp,double iza)
{
  return izero(alp*sqrt(1-4*n*n/((len-1.0)*(len-1.0))))/iza;
}

double sinc(double x)
{
  return x == 0 ? 1 : sin(x)/x;
}

double hn_lpf(int n,double lpf,double fs)
{
  double t = 1/fs;
  double omega = 2*M_PI*lpf;
  return 2*lpf*t*sinc(n*omega*t);
}

void usage(void)
{
  printf("http://shibatch.sourceforge.net/\n\n");
  printf("usage: ssrc [<options>] <source wav file> <destination wav file>\n");
  printf("options : --rate <sampling rate>     output sample rate\n");
  printf("          --att <attenuation(dB)>    attenuate signal\n");
  printf("          --bits <number of bits>    output quantization bit length\n");
  printf("          --tmpfile <file name>      specify temporal file\n");
  printf("          --twopass                  two pass processing to avoid clipping\n");
  printf("          --normalize                normalize the wave file\n");
  printf("          --quiet                    nothing displayed except error\n");
  printf("          --dither [<type>]          dither options\n");
  printf("                                       0    : ATH-based noise shaping, low intensity\n");
  printf("                                       2    : ATH-based noise shaping, mid intensity\n");
  printf("                                       6    : ATH-based noise shaping, high intensity\n");
  printf("                                       99   : No noise shaping\n");
  printf("                                       help : Show all available noise shapers\n");
  printf("          --pdf <type> [<amp>]       select p.d.f. of noise\n");
  printf("                                       0 : Rectangular\n");
  printf("                                       1 : Triangular\n");
  //printf("                                       2 : Gaussian\n");
  printf("                                       3 : Two-level (experimental)\n");
  printf("          --profile <type>           specify profile\n");
  printf("                                       short  : shorter filter length\n");
  printf("                                       normal : the default quality\n");
  printf("                                       long   : longer filter length\n");
}

void fmterr(int x)
{
  fprintf(stderr,"unknown error %d\n",x);
  exit(-1);
}

void setstarttime(void)
{
  starttime = time(NULL);
  lastshowed = 0;
  lastshowed2 = -1;
}

void showprogress(double p)
{
  int eta,pc;
  time_t t;
  if (quiet) return;

  t = time(NULL)-starttime;
  if (p == 0) eta = 0; else eta = t*(1-p)/p;

  pc = (int)(p*100);

  if (pc != lastshowed2 || t != lastshowed) {
    printf(" %3d%% processed",pc);
    lastshowed2 = pc;
  }
  if(t != lastshowed) {
    printf(", ETA =%4dsec",eta);
    lastshowed = t;
  }
  printf("\r");
  fflush(stdout);
}

int64_t gcd(int64_t x, int64_t y)
{
    int64_t t;

    while (y != 0) {
        t = x % y;  x = y;  y = t;
    }
    return x;
}

double upsample(FILE *fpi,FILE *fpo,int nch,int bps,int dbps,int64_t sfrq,int64_t dfrq,double gain,unsigned int chunklen,int twopass,int dither)
{
  int64_t frqgcd,osf;
  int64_t fs1, fs2, n2;
  int64_t n1,n1x,n1y,n2b;

  REAL **stage1,*stage2;
  int filter2len;
  int *f1order,*f1inc;
  struct SleefDFT *dftf = NULL, *dftb = NULL;
  uint8_t *rawinbuf,*rawoutbuf;
  REAL *inbuf,*outbuf;
  REAL **buf1,**buf2;
  double peak=0;
  int spcount = 0;
  int i,j;

  filter2len = 1; /* stage 2 filter length */

  /* Make stage 1 filter */

  {
    double aa = AA; /* stop band attenuation(dB) */
    double lpf,delta,d,df,alp,iza;
    double guard = 2;

    frqgcd = gcd(sfrq,dfrq);

    fs1 = (int64_t)sfrq / frqgcd * dfrq;

    if (fs1/dfrq == 1) osf = 1;
    else if (fs1/dfrq % 2 == 0) osf = 2;
    else if (fs1/dfrq % 3 == 0) osf = 3;
    else {
      fprintf(stderr,"Resampling from %dHz to %dHz is not supported.\n",(int)sfrq, (int)dfrq);
      fprintf(stderr,"%d/gcd(%d,%d)=%d must be divided by 2 or 3.\n",(int)sfrq,(int)sfrq,(int)dfrq,(int)(fs1/dfrq));
      exit(-1);
    }

    df = (dfrq*osf/2 - sfrq/2) * 2 / guard;
    lpf = sfrq/2 + (dfrq*osf/2 - sfrq/2)/guard;

    delta = pow(10,-aa/20);
    if (aa <= 21) d = 0.9222; else d = (aa-7.95)/14.36;

    n1 = fs1/df*d+1;
    if (n1 % 2 == 0) n1++;

    alp = alpha(aa);
    iza = izero(alp);

    n1y = fs1/sfrq;
    n1x = n1/n1y+1;

    f1order = calloc(n1y*osf,sizeof(int));
    for(i=0;i<n1y*osf;i++) {
      f1order[i] = fs1/sfrq-(i*(fs1/(dfrq*osf)))%(fs1/sfrq);
      if (f1order[i] == fs1/sfrq) f1order[i] = 0;
    }

    f1inc = calloc(n1y*osf,sizeof(int));
    for(i=0;i<n1y*osf;i++) {
      f1inc[i] = f1order[i] < fs1/(dfrq*osf) ? nch : 0;
    }

    stage1 = malloc(n1y*sizeof(REAL *));
    stage1[0] = calloc(n1x*n1y,sizeof(REAL));

    for(i=1;i<n1y;i++) {
      stage1[i] = &(stage1[0][n1x*i]);
      for(j=0;j<n1x;j++) stage1[i][j] = 0;
    }

    for(i=-(n1/2);i<=n1/2;i++)
      {
	stage1[(i+n1/2)%n1y][(i+n1/2)/n1y] = win(i,n1,alp,iza)*hn_lpf(i,lpf,fs1)*fs1/sfrq;
      }
  }

  /* Make stage 2 filter */

  {
    double aa = AA; /* stop band attenuation(dB) */
    double lpf,delta,d,df,alp,iza;
    int ipsize,wsize;

    delta = pow(10,-aa/20);
    if (aa <= 21) d = 0.9222; else d = (aa-7.95)/14.36;

    fs2 = (int64_t)dfrq*osf;

    for(i=1;;i = i * 2)
      {
	n2 = filter2len * i;
	if (n2 % 2 == 0) n2--;
	df = (fs2*d)/(n2-1);
	lpf = sfrq/2;
	if (df < DF) break;
      }

    alp = alpha(aa);
    iza = izero(alp);

    for(n2b=1;n2b<n2;n2b*=2);
    n2b *= 2;

    stage2 = SleefDFT_malloc(n2b * sizeof(REAL));

    for(i=-(n2/2);i<=n2/2;i++) {
      stage2[i+n2/2] = win(i,n2,alp,iza)*hn_lpf(i,lpf,fs2)/n2b*2;
    }

    dftf = SleefDFT_init(SLEEF_MODE_REAL | SLEEF_MODE_ALT | SLEEF_MODE_FORWARD , n2b);
    dftb = SleefDFT_init(SLEEF_MODE_REAL | SLEEF_MODE_ALT | SLEEF_MODE_BACKWARD, n2b);

    SleefDFT_execute(dftf, stage2, stage2);
  }

  /* Apply filters */

  setstarttime();

  {
    int n2b2 = n2b/2;
    int rp;            // Position of the next sample read from inbuf at fs1
    int ds;            // Number of samples disposed next at sfrq
    int nsmplwrt1 = 0; // Number of samples that will be sent to the stage2 filter, that is
                       // calculated from the actual number of samples read from the file
    int nsmplwrt2 = 0; // Number of samples that will be sent to the stage2 filter, that is
                       // calculated from the actual number of samples read from the file
    int s1p;           // Residual of the number of samples output from the stage1 filter divided by (n1y * osf)
    int init,ending;
    int64_t sumread,sumwrite;
    int osc;
    REAL *ip,*ip_backup;
    int s1p_backup,osc_backup;
    int k,ch,p;
    int inbuflen;
    int delay = 0;

    buf1 = malloc(nch*sizeof(REAL *));
    for(i=0;i<nch;i++)
      {
	buf1[i] = calloc((n2b2/osf+1),sizeof(REAL));
	for(j=0;j<(n2b2/osf+1);j++) buf1[i][j] = 0;
      }

    buf2 = malloc(sizeof(REAL *)*nch);
    for(i=0;i<nch;i++) buf2[i] = SleefDFT_malloc(n2b * sizeof(REAL));

    rawinbuf  = calloc(nch*(n2b2+n1x+2), bps);
    rawoutbuf = calloc(nch*(n2b2/osf+1), dbps);

    inbuf  = calloc(nch*(n2b2+n1x+2),sizeof(REAL));
    outbuf = calloc(nch*(n2b2/osf+1),sizeof(REAL));

    s1p = 0;
    rp  = 0;
    ds  = 0;
    osc = 0;

    init = 1;
    ending = 0;
    inbuflen = n1/2/(fs1/sfrq)+1;
    delay = (double)n2/2/(fs2/dfrq);

    sumread = sumwrite = 0;

    for(;;)
      {
	int nsmplread,toberead,toberead2;

	toberead2 = toberead = ceil((double)n2b2*sfrq/(dfrq*osf))+1+n1x-inbuflen;
	if (toberead+sumread > chunklen) {
	  toberead = chunklen-sumread;
	}

	nsmplread = fread(rawinbuf,1,bps*nch*toberead,fpi);
	nsmplread /= bps*nch;

	switch(bps)
	  {
	  case 1:
	    for(i = 0; i < nsmplread * nch; i++)
	      inbuf[nch * inbuflen + i] =
		(1 / (REAL)0x7f) * ((REAL)((uint8_t *)rawinbuf)[i]-128);
	    break;

	  case 2:
#ifndef BIGENDIAN
	    for(i=0;i<nsmplread*nch;i++)
	      inbuf[nch*inbuflen+i] = (1/(REAL)0x7fff)*(REAL)((int16_t *)rawinbuf)[i];
#else
	    for(i=0;i<nsmplread*nch;i++) {
	      inbuf[nch*inbuflen+i] = (1/(REAL)0x7fff)*
		(((int)rawinbuf[i*2]) |
		 (((int)((int8_t *)rawinbuf)[i*2+1]) << 8));
	    }
#endif
	    break;

	  case 3:
	    for(i=0;i<nsmplread*nch;i++) {
	      inbuf[nch*inbuflen+i] = (1/(REAL)0x7fffff)*
		((((int)rawinbuf[i*3  ]) << 0 ) |
		 (((int)rawinbuf[i*3+1]) << 8 ) |
		 (((int)((int8_t *)rawinbuf)[i*3+2]) << 16));
	    }
	    break;

	  case 4:
	    for(i=0;i<nsmplread*nch;i++) {
	      inbuf[nch*inbuflen+i] = (1/(REAL)0x7fffffff)*
		((((int)rawinbuf[i*4  ]) << 0 ) |
		 (((int)rawinbuf[i*4+1]) << 8 ) |
		 (((int)rawinbuf[i*4+2]) << 16) |
		 (((int)((int8_t *)rawinbuf)[i*4+3]) << 24));
	    }
	    break;
	  }

	assert(i <= nch*(n2b2+n1x+2));
	for(;i<nch*toberead2;i++) inbuf[nch*inbuflen+i] = 0;

	inbuflen += toberead2;

	sumread += nsmplread;

	ending = feof(fpi) || sumread >= chunklen;

	nsmplwrt1 = n2b2;

	// apply stage 1 filter

	ip = &inbuf[((sfrq*(rp-1)+fs1)/fs1)*nch];

	s1p_backup = s1p;
	ip_backup  = ip;
	osc_backup = osc;

	for(ch=0;ch<nch;ch++)
	  {
	    REAL *op = &outbuf[ch];
	    int fdo = fs1/(dfrq*osf),no = n1y*osf;

	    s1p = s1p_backup; ip = ip_backup+ch;

	    switch(n1x)
	      {
	      case 7:
		for(p=0;p<nsmplwrt1;p++)
		  {
		    int s1o = f1order[s1p];

		    buf2[ch][p] =
		      stage1[s1o][0] * *(ip+0*nch)+
		      stage1[s1o][1] * *(ip+1*nch)+
		      stage1[s1o][2] * *(ip+2*nch)+
		      stage1[s1o][3] * *(ip+3*nch)+
		      stage1[s1o][4] * *(ip+4*nch)+
		      stage1[s1o][5] * *(ip+5*nch)+
		      stage1[s1o][6] * *(ip+6*nch);

		    ip += f1inc[s1p];

		    s1p++;
		    if (s1p == no) s1p = 0;
		  }
		break;

	      case 9:
		for(p=0;p<nsmplwrt1;p++)
		  {
		    int s1o = f1order[s1p];

		    buf2[ch][p] =
		      stage1[s1o][0] * *(ip+0*nch)+
		      stage1[s1o][1] * *(ip+1*nch)+
		      stage1[s1o][2] * *(ip+2*nch)+
		      stage1[s1o][3] * *(ip+3*nch)+
		      stage1[s1o][4] * *(ip+4*nch)+
		      stage1[s1o][5] * *(ip+5*nch)+
		      stage1[s1o][6] * *(ip+6*nch)+
		      stage1[s1o][7] * *(ip+7*nch)+
		      stage1[s1o][8] * *(ip+8*nch);

		    ip += f1inc[s1p];

		    s1p++;
		    if (s1p == no) s1p = 0;
		  }
		break;

	      default:
		for(p=0;p<nsmplwrt1;p++)
		  {
		    REAL tmp = 0;
		    REAL *ip2=ip;

		    int s1o = f1order[s1p];

		    for(i=0;i<n1x;i++)
		      {
			tmp += stage1[s1o][i] * *ip2;
			ip2 += nch;
		      }
		    buf2[ch][p] = tmp;

		    ip += f1inc[s1p];

		    s1p++;
		    if (s1p == no) s1p = 0;
		  }
		break;
	      }

	    osc = osc_backup;

	    // apply stage 2 filter

	    for(p=nsmplwrt1;p<n2b;p++) buf2[ch][p] = 0;

	    SleefDFT_execute(dftf, buf2[ch], buf2[ch]);

	    buf2[ch][0] = stage2[0]*buf2[ch][0];
	    buf2[ch][1] = stage2[1]*buf2[ch][1]; 

	    for(i=1;i<n2b/2;i++)
	      {
		REAL re,im;

		re = stage2[i*2  ]*buf2[ch][i*2] - stage2[i*2+1]*buf2[ch][i*2+1];
		im = stage2[i*2+1]*buf2[ch][i*2] + stage2[i*2  ]*buf2[ch][i*2+1];

		buf2[ch][i*2  ] = re;
		buf2[ch][i*2+1] = im;
	      }

	    SleefDFT_execute(dftb, buf2[ch], buf2[ch]);

	    for(i=osc,j=0;i<n2b2;i+=osf,j++)
	      {
		REAL f = (buf1[ch][j] + buf2[ch][i]);
		op[j*nch] = f;
	      }

	    nsmplwrt2 = j;

	    osc = i - n2b2;

	    for(j=0;i<n2b;i+=osf,j++)
	      buf1[ch][j] = buf2[ch][i];
	  }

	rp += nsmplwrt1 * (sfrq / frqgcd) / osf;

	if (twopass) {
	  for(i=0;i<nsmplwrt2*nch;i++)
	    {
	      REAL f = outbuf[i] > 0 ? outbuf[i] : -outbuf[i];
	      peak = peak < f ? f : peak;
	      ((REAL *)rawoutbuf)[i] = outbuf[i];
	    }
	} else {
	  switch(dbps)
	    {
	    case 1:
	      {
		REAL gain2 = gain * (REAL)0x7f;
		ch = 0;

		for(i=0;i<nsmplwrt2*nch;i++)
		  {
		    int s;

		    if (dither != -1) {
		      s = do_shaping(outbuf[i]*gain2,&peak,dither,ch);
		    } else {
		      s = RINT(outbuf[i]*gain2);

		      if (s < -0x80) {
			double d = (double)s/-0x80;
			peak = peak < d ? d : peak;
			s = -0x80;
		      }
		      if (0x7f <  s) {
			double d = (double)s/ 0x7f;
			peak = peak < d ? d : peak;
			s =  0x7f;
		      }
		    }

		    ((uint8_t *)rawoutbuf)[i] = s + 0x80;

		    ch++;
		    if (ch == nch) ch = 0;
		  }
	      }
	    break;

	    case 2:
	      {
		REAL gain2 = gain * (REAL)0x7fff;
		ch = 0;

		for(i=0;i<nsmplwrt2*nch;i++)
		  {
		    int s;

		    if (dither != -1) {
		      s = do_shaping(outbuf[i]*gain2,&peak,dither,ch);
		    } else {
		      s = RINT(outbuf[i]*gain2);

		      if (s < -0x8000) {
			double d = (double)s/-0x8000;
			peak = peak < d ? d : peak;
			s = -0x8000;
		      }
		      if (0x7fff <  s) {
			double d = (double)s/ 0x7fff;
			peak = peak < d ? d : peak;
			s =  0x7fff;
		      }
		    }

#ifndef BIGENDIAN
		    ((int16_t *)rawoutbuf)[i] = s;
#else
		    ((int8_t *)rawoutbuf)[i*2  ] = s & 255; s >>= 8;
		    ((int8_t *)rawoutbuf)[i*2+1] = s & 255;
#endif
		    ch++;
		    if (ch == nch) ch = 0;
		  }
	      }
	    break;

	    case 3:
	      {
		REAL gain2 = gain * (REAL)0x7fffff;
		ch = 0;

		for(i=0;i<nsmplwrt2*nch;i++)
		  {
		    int s;

		    if (dither != -1) {
		      s = do_shaping(outbuf[i]*gain2,&peak,dither,ch);
		    } else {
		      s = RINT(outbuf[i]*gain2);

		      if (s < -0x800000) {
			double d = (double)s/-0x800000;
			peak = peak < d ? d : peak;
			s = -0x800000;
		      }
		      if (0x7fffff <  s) {
			double d = (double)s/ 0x7fffff;
			peak = peak < d ? d : peak;
			s =  0x7fffff;
		      }
		    }

		    ((int8_t *)rawoutbuf)[i*3  ] = s & 255; s >>= 8;
		    ((int8_t *)rawoutbuf)[i*3+1] = s & 255; s >>= 8;
		    ((int8_t *)rawoutbuf)[i*3+2] = s & 255;

		    ch++;
		    if (ch == nch) ch = 0;
		  }
	      }
	    break;

	    }
	}

	if (!init) {
	  if (ending) {
	    if ((double)sumread*dfrq/sfrq+2 > sumwrite+nsmplwrt2) {
	      if (dbps*nch*nsmplwrt2 != fwrite(rawoutbuf,1,dbps*nch*nsmplwrt2,fpo)) {
		fprintf(stderr,"ERROR: -41\n");
		abort();
	      }
	      sumwrite += nsmplwrt2;
	    } else {
	      if (dbps*nch*(floor((double)sumread*dfrq/sfrq)+2-sumwrite) != 
		  fwrite(rawoutbuf,1,dbps*nch*(floor((double)sumread*dfrq/sfrq)+2-sumwrite),fpo)) {
		fprintf(stderr,"fwrite error(2).\n");
		abort();
	      }
	      break;
	    }
	  } else {
	    if (dbps*nch*nsmplwrt2 != fwrite(rawoutbuf,1,dbps*nch*nsmplwrt2,fpo)) {
	      fprintf(stderr,"fwrite error(3).\n");
	      abort();
	    }
	    sumwrite += nsmplwrt2;
	  }
	} else {
	  int pos,len;
	  if (nsmplwrt2 < delay) {
	    delay -= nsmplwrt2;
	  } else {
	    if (ending) {
	      if ((double)sumread*dfrq/sfrq+2 > sumwrite+nsmplwrt2-delay) {
		if (dbps*nch*(nsmplwrt2-delay) != fwrite(rawoutbuf+dbps*nch*delay,1,dbps*nch*(nsmplwrt2-delay),fpo)) {
		  fprintf(stderr,"fwrite error(4).\n");
		  abort();
		}
		sumwrite += nsmplwrt2-delay;
	      } else {
		if (dbps*nch*(floor((double)sumread*dfrq/sfrq)+2-sumwrite-delay) !=
		    fwrite(rawoutbuf+dbps*nch*delay,1,dbps*nch*(floor((double)sumread*dfrq/sfrq)+2-sumwrite-delay),fpo)) {
		  fprintf(stderr,"fwrite error(5).\n");
		  abort();
		}
		break;
	      }
	    } else {
	      if (dbps*nch*(nsmplwrt2-delay) != fwrite(rawoutbuf+dbps*nch*delay,1,dbps*nch*(nsmplwrt2-delay),fpo)) {
		fprintf(stderr,"fwrite error(6).\n");
		abort();
	      }
	      sumwrite += nsmplwrt2-delay;
	      init = 0;
	    }
	  }
	}

	{
	  int ds = (rp-1)/(fs1/sfrq);

	  assert(inbuflen >= ds);

	  memmove(inbuf,inbuf+nch*ds,sizeof(REAL)*nch*(inbuflen-ds));
	  inbuflen -= ds;
	  rp -= ds*(fs1/sfrq);
	}

	if ((spcount++ & 7) == 7) showprogress((double)sumread / chunklen);
      }
  }

  showprogress(1);

  free(f1order);
  free(f1inc);
  free(stage1[0]);
  free(stage1);
  SleefDFT_free(stage2);
  SleefDFT_dispose(dftf);
  SleefDFT_dispose(dftb);
  for(i=0;i<nch;i++) free(buf1[i]);
  free(buf1);
  for(i=0;i<nch;i++) SleefDFT_free(buf2[i]);
  free(buf2);
  free(inbuf);
  free(outbuf);
  free(rawinbuf);
  free(rawoutbuf);

  if (dither != -1) {
    static int64_t dlmin[] = { 0, -0x80, -0x8000, -0x800000, -0x80000000ULL };
    static int64_t dlmax[] = { 0,  0x7f,  0x7fff,  0x7fffff,  0x7fffffffULL };
    double min = 0, max = 0;
    for(i=0;i<nch;i++) {
      double p[2];
      SSRCDither_getPeaks(shaper[i], p);
      if (p[1] < min) min = p[1];
      if (p[0] > max) max = p[0];
    }
    if (min < dlmin[dbps]) peak = min / dlmin[dbps];
    if (max > dlmax[dbps]) peak = max / dlmax[dbps];
  }

  return peak;
}

double downsample(FILE *fpi,FILE *fpo,int nch,int bps,int dbps,int64_t sfrq,int64_t dfrq,double gain,unsigned int chunklen,int twopass,int dither)
{
  int64_t frqgcd,osf;
  int64_t fs1, fs2;
  int64_t n2,n2x,n2y,n1,n1b;

  REAL *stage1,**stage2;
  int filter1len;
  int *f2order,*f2inc;
  struct SleefDFT *dftf = NULL, *dftb = NULL;
  uint8_t *rawinbuf,*rawoutbuf;
  REAL *inbuf,*outbuf;
  REAL **buf1,**buf2;
  int i,j;
  int spcount = 0;
  double peak=0;

  filter1len = 1; /* stage 1 filter length */

  /* Make stage 1 filter */

  {
    double aa = AA; /* stop band attenuation(dB) */
    double lpf,delta,d,df,alp,iza;
    int ipsize,wsize;

    frqgcd = gcd(sfrq,dfrq);

    if (dfrq/frqgcd == 1) osf = 1;
    else if (dfrq/frqgcd % 2 == 0) osf = 2;
    else if (dfrq/frqgcd % 3 == 0) osf = 3;
    else {
      fprintf(stderr,"Resampling from %dHz to %dHz is not supported.\n",(int)sfrq,(int)dfrq);
      fprintf(stderr,"%d/gcd(%d,%d)=%d must be divided by 2 or 3.\n",(int)dfrq,(int)sfrq,(int)dfrq,(int)(dfrq/frqgcd));
      exit(-1);
    }

    fs1 = (int64_t)sfrq*osf;

    delta = pow(10,-aa/20);
    if (aa <= 21) d = 0.9222; else d = (aa-7.95)/14.36;

    n1 = filter1len;
    for(i=1;;i = i * 2)
      {
	n1 = filter1len * i;
	if (n1 % 2 == 0) n1--;
	df = (fs1*d)/(n1-1);
	lpf = (dfrq-df)/2;
	if (df < DF) break;
      }

    alp = alpha(aa);
    iza = izero(alp);

    for(n1b=1;n1b<n1;n1b*=2);
    n1b *= 2;

    stage1 = SleefDFT_malloc(n1b * sizeof(REAL));

    for(i=0;i<n1b;i++) stage1[i] = 0;

    for(i=-(n1/2);i<=n1/2;i++) {
      stage1[i+n1/2] = win(i,n1,alp,iza)*hn_lpf(i,lpf,fs1)*fs1/sfrq/n1b*2;
    }

    dftf = SleefDFT_init(SLEEF_MODE_REAL | SLEEF_MODE_ALT | SLEEF_MODE_FORWARD , n1b);
    dftb = SleefDFT_init(SLEEF_MODE_REAL | SLEEF_MODE_ALT | SLEEF_MODE_BACKWARD, n1b);

    SleefDFT_execute(dftf, stage1, stage1);
  }

  /* Make stage 2 filter */

  if (osf == 1) {
    fs2 = (int64_t)sfrq/frqgcd*dfrq;
    n2 = 1;
    n2y = n2x = 1;
    f2order = calloc(n2y,sizeof(int));
    f2order[0] = 0;
    f2inc = calloc(n2y,sizeof(int));
    f2inc[0] = sfrq/dfrq;
    stage2 = malloc(sizeof(REAL *)*n2y);
    stage2[0] = calloc(n2x*n2y,sizeof(REAL));
    stage2[0][0] = 1;
  } else {
    double aa = AA; /* stop band attenuation(dB) */
    double lpf,delta,d,df,alp,iza;
    double guard = 2;

    fs2 = (int64_t)sfrq / frqgcd * dfrq;

    df = (fs1/2 - sfrq/2) * 2 / guard;
    lpf = sfrq/2 + (fs1/2 - sfrq/2)/guard;

    delta = pow(10,-aa/20);
    if (aa <= 21) d = 0.9222; else d = (aa-7.95)/14.36;

    n2 = (int64_t)fs2/df*d+1;
    if (n2 % 2 == 0) n2++;

    alp = alpha(aa);
    iza = izero(alp);

    n2y = fs2/fs1; // non-zero samples are at intervals of n2y at fs2
    n2x = n2/n2y+1;

    f2order = calloc(n2y,sizeof(int));
    for(i=0;i<n2y;i++) {
      f2order[i] = fs2/fs1-(i*(fs2/dfrq))%(fs2/fs1);
      if (f2order[i] == fs2/fs1) f2order[i] = 0;
    }

    f2inc = calloc(n2y,sizeof(int));
    for(i=0;i<n2y;i++) {
      f2inc[i] = (fs2/dfrq-f2order[i])/(fs2/fs1)+1;
      if (f2order[i+1==n2y ? 0 : i+1] == 0) f2inc[i]--;
    }

    stage2 = malloc(sizeof(REAL *)*n2y);
    stage2[0] = calloc(n2x*n2y,sizeof(REAL));

    for(i=1;i<n2y;i++) {
      stage2[i] = &(stage2[0][n2x*i]);
      for(j=0;j<n2x;j++) stage2[i][j] = 0;
    }

    for(i=-(n2/2);i<=n2/2;i++)
      {
	stage2[(i+n2/2)%n2y][(i+n2/2)/n2y] = win(i,n2,alp,iza)*hn_lpf(i,lpf,fs2)*fs2/fs1;
      }
  }

  /* Apply filters */

  setstarttime();

  {
    int n1b2 = n1b/2;
    int rp;        // inbufのfs1での次に読むサンプルの場所を保持
    int rps;       // rpを(fs1/sfrq=osf)で割った余り
    int rp2;       // buf2のfs2での次に読むサンプルの場所を保持
    int ds;        // 次にdisposeするsfrqでのサンプル数
    int nsmplwrt1 = 0; // 実際にファイルからinbufに読み込まれた値から計算した
                   // stage2 filterに渡されるサンプル数
    int nsmplwrt2 = 0; // 実際にファイルからinbufに読み込まれた値から計算した
                   // stage2 filterに渡されるサンプル数
    int s2p;       // stage1 filterから出力されたサンプルの数をn1y*osfで割った余り
    int init,ending;
    int osc;
    REAL *bp,*bp_backup; // rp2から計算される．buf2の次に読むサンプルの位置
    int rps_backup,s2p_backup,osc_backup;
    int k,ch,p;
    int inbuflen=0;
    int64_t sumread,sumwrite;
    int delay = 0;
    REAL *op;

    //    |....B....|....C....|   buf1      n1b2+n1b2
    //|.A.|....D....|             buf2  n2x+n1b2
    //
    // まずinbufからBにosf倍サンプリングしながらコピー
    // Cはクリア
    // BCにstage 1 filterをかける
    // DにBを足す
    // ADにstage 2 filterをかける
    // Dの後ろをAに移動
    // CをDにコピー

    buf1 = malloc(sizeof(REAL *)*nch);
    for(i=0;i<nch;i++)
      buf1[i] = SleefDFT_malloc(sizeof(REAL) * n1b);

    buf2 = malloc(sizeof(REAL *)*nch);
    for(i=0;i<nch;i++) {
      buf2[i] = calloc(n2x+1+n1b2,sizeof(REAL));
      for(j=0;j<n2x+1+n1b2;j++) buf2[i][j] = 0;
    }

    rawinbuf  = calloc(nch*(n1b2/osf+osf+1),bps);
    rawoutbuf = calloc(((double)n1b2*sfrq/dfrq+1),dbps*nch);
    inbuf = calloc(nch*(n1b2/osf+osf+1),sizeof(REAL));
    outbuf = calloc(nch*((double)n1b2*sfrq/dfrq+1),sizeof(REAL));

    op = outbuf;

    s2p = 0;
    rp  = 0;
    rps = 0;
    ds  = 0;
    osc = 0;
    rp2 = 0;

    init = 1;
    ending = 0;
    delay = (double)n1/2/((double)fs1/dfrq)+(double)n2/2/((double)fs2/dfrq);

    sumread = sumwrite = 0;

    for(;;)
      {
	int nsmplread;
	int toberead;

	toberead = (n1b2-rps-1)/osf+1;
	if (toberead+sumread > chunklen) {
	  toberead = chunklen-sumread;
	}

	nsmplread = fread(rawinbuf,1,bps*nch*toberead,fpi);
	nsmplread /= bps*nch;

	switch(bps)
	  {
	  case 1:
	    for(i = 0; i < nsmplread * nch; i++)
	      inbuf[nch * inbuflen + i] =
		(1 / (REAL)0x7f) * ((REAL)((uint8_t *)rawinbuf)[i]-128);
	    break;

	  case 2:
#ifndef BIGENDIAN
	    for(i=0;i<nsmplread*nch;i++)
	      inbuf[nch*inbuflen+i] = (1/(REAL)0x7fff)*(REAL)((int16_t *)rawinbuf)[i];
#else
	    for(i=0;i<nsmplread*nch;i++) {
	      inbuf[nch*inbuflen+i] = (1/(REAL)0x7fff)*
		(((int)rawinbuf[i*2]) |
		 (((int)((int8_t *)rawinbuf)[i*2+1]) << 8));
	    }
#endif
	    break;

	  case 3:
	    for(i=0;i<nsmplread*nch;i++) {
	      inbuf[nch*inbuflen+i] = (1/(REAL)0x7fffff)*
		((((int)rawinbuf[i*3  ]) << 0 ) |
		 (((int)rawinbuf[i*3+1]) << 8 ) |
		 (((int)((int8_t *)rawinbuf)[i*3+2]) << 16));
	    }
	    break;

	  case 4:
	    for(i=0;i<nsmplread*nch;i++) {
	      inbuf[nch*inbuflen+i] = (1/(REAL)0x7fffffff)*
		((((int)rawinbuf[i*4  ]) << 0 ) |
		 (((int)rawinbuf[i*4+1]) << 8 ) |
		 (((int)rawinbuf[i*4+2]) << 16) |
		 (((int)((int8_t *)rawinbuf)[i*4+3]) << 24));
	    }
	    break;
	  }
	assert(i <= nch*(n1b2/osf+osf+1));

	for(;i<nch*toberead;i++) inbuf[i] = 0;

	sumread += nsmplread;

	ending = feof(fpi) || sumread >= chunklen;

	rps_backup = rps;
	s2p_backup = s2p;

	for(ch=0;ch<nch;ch++)
	  {
	    rps = rps_backup;

	    for(k=0;k<rps;k++) buf1[ch][k] = 0;

	    for(i=rps,j=0;i<n1b2;i+=osf,j++)
	      {
		assert(j < ((n1b2-rps-1)/osf+1));

		buf1[ch][i] = inbuf[j*nch+ch];

		for(k=i+1;k<i+osf;k++) buf1[ch][k] = 0;
	      }

	    assert(j == ((n1b2-rps-1)/osf+1));

	    for(k=n1b2;k<n1b;k++) buf1[ch][k] = 0;

	    rps = i - n1b2;
	    rp += j;

	    SleefDFT_execute(dftf, buf1[ch], buf1[ch]);

	    buf1[ch][0] = stage1[0]*buf1[ch][0];
	    buf1[ch][1] = stage1[1]*buf1[ch][1]; 

	    for(i=1;i<n1b2;i++)
	      {
		REAL re,im;

		re = stage1[i*2  ]*buf1[ch][i*2] - stage1[i*2+1]*buf1[ch][i*2+1];
		im = stage1[i*2+1]*buf1[ch][i*2] + stage1[i*2  ]*buf1[ch][i*2+1];

		buf1[ch][i*2  ] = re;
		buf1[ch][i*2+1] = im;
	      }

	    SleefDFT_execute(dftb, buf1[ch], buf1[ch]);

	    for(i=0;i<n1b2;i++) {
	      buf2[ch][n2x+1+i] += buf1[ch][i];
	    }

	    {
	      int t1 = rp2/(fs2/fs1);
	      if (rp2%(fs2/fs1) != 0) t1++;

	      bp = &(buf2[ch][t1]);
	    }

	    s2p = s2p_backup;

	    for(p=0;bp-buf2[ch]<n1b2+1;p++)
	      {
		REAL tmp = 0;
		REAL *bp2;
		int s;
		int s2o;

		bp2 = bp;
		s2o = f2order[s2p];
		bp += f2inc[s2p];
		s2p++;

		if (s2p == n2y) s2p = 0;

		assert((bp2-&(buf2[ch][0]))*(fs2/fs1)-(rp2+p*(fs2/dfrq)) == s2o);

		for(i=0;i<n2x;i++)
		  tmp += stage2[s2o][i] * *bp2++;

		op[p*nch+ch] = tmp;
	      }

	    nsmplwrt2 = p;
	  }

	rp2 += nsmplwrt2 * (fs2 / dfrq);

	if (twopass) {
	  for(i=0;i<nsmplwrt2*nch;i++)
	    {
	      REAL f = outbuf[i] > 0 ? outbuf[i] : -outbuf[i];
	      peak = peak < f ? f : peak;
	      ((REAL *)rawoutbuf)[i] = outbuf[i];
	    }
	} else {
	  switch(dbps)
	    {
	    case 1:
	      {
		REAL gain2 = gain * (REAL)0x7f;
		ch = 0;

		for(i=0;i<nsmplwrt2*nch;i++)
		  {
		    int s;

		    if (dither != -1) {
		      s = do_shaping(outbuf[i]*gain2,&peak,dither,ch);
		    } else {
		      s = RINT(outbuf[i]*gain2);

		      if (s < -0x80) {
			double d = (double)s/-0x80;
			peak = peak < d ? d : peak;
			s = -0x80;
		      }
		      if (0x7f <  s) {
			double d = (double)s/ 0x7f;
			peak = peak < d ? d : peak;
			s =  0x7f;
		      }
		    }

		    ((uint8_t *)rawoutbuf)[i] = s + 0x80;

		    ch++;
		    if (ch == nch) ch = 0;
		  }
	      }
	    break;

	    case 2:
	      {
		REAL gain2 = gain*(REAL)0x7fff;
		ch = 0;

		for(i=0;i<nsmplwrt2*nch;i++)
		  {
		    int s;

		    if (dither != -1) {
		      s = do_shaping(outbuf[i]*gain2,&peak,dither,ch);
		    } else {
		      s = RINT(outbuf[i]*gain2);

		      if (s < -0x8000) {
			double d = (double)s/-0x8000;
			peak = peak < d ? d : peak;
			s = -0x8000;
		      }
		      if (0x7fff <  s) {
			double d = (double)s/ 0x7fff;
			peak = peak < d ? d : peak;
			s =  0x7fff;
		      }
		    }

#ifndef BIGENDIAN
		    ((int16_t *)rawoutbuf)[i] = s;
#else
		    ((int8_t *)rawoutbuf)[i*2  ] = s & 255; s >>= 8;
		    ((int8_t *)rawoutbuf)[i*2+1] = s & 255;
#endif

		    ch++;
		    if (ch == nch) ch = 0;
		  }
	      }
	    break;

	    case 3:
	      {
		REAL gain2 = gain * (REAL)0x7fffff;
		ch = 0;

		for(i=0;i<nsmplwrt2*nch;i++)
		  {
		    int s;

		    if (dither != -1) {
		      s = do_shaping(outbuf[i]*gain2,&peak,dither,ch);
		    } else {
		      s = RINT(outbuf[i]*gain2);

		      if (s < -0x800000) {
			double d = (double)s/-0x800000;
			peak = peak < d ? d : peak;
			s = -0x800000;
		      }
		      if (0x7fffff <  s) {
			double d = (double)s/ 0x7fffff;
			peak = peak < d ? d : peak;
			s =  0x7fffff;
		      }
		    }

		    ((int8_t *)rawoutbuf)[i*3  ] = s & 255; s >>= 8;
		    ((int8_t *)rawoutbuf)[i*3+1] = s & 255; s >>= 8;
		    ((int8_t *)rawoutbuf)[i*3+2] = s & 255;

		    ch++;
		    if (ch == nch) ch = 0;
		  }
	      }
	    break;
	    }
	}

	if (!init) {
	  if (ending) {
	    if ((double)sumread*dfrq/sfrq+2 > sumwrite+nsmplwrt2) {
	      if (dbps*nch*nsmplwrt2 != fwrite(rawoutbuf,1,dbps*nch*nsmplwrt2,fpo)) {
	      }
	      sumwrite += nsmplwrt2;
	    } else {
	      fwrite(rawoutbuf,1,dbps*nch*(floor((double)sumread*dfrq/sfrq)+2-sumwrite),fpo);
	      break;
	    }
	  } else {
	    fwrite(rawoutbuf,1,dbps*nch*nsmplwrt2,fpo);
	    sumwrite += nsmplwrt2;
	  }
	} else {
	  int pos,len;
	  if (nsmplwrt2 < delay) {
	    delay -= nsmplwrt2;
	  } else {
	    if (ending) {
	      if ((double)sumread*dfrq/sfrq+2 > sumwrite+nsmplwrt2-delay) {
		fwrite(rawoutbuf+dbps*nch*delay,1,dbps*nch*(nsmplwrt2-delay),fpo);
		sumwrite += nsmplwrt2-delay;
	      } else {
		fwrite(rawoutbuf+dbps*nch*delay,1,dbps*nch*(floor((double)sumread*dfrq/sfrq)+2-sumwrite-delay),fpo);
		break;
	      }
	    } else {
	      fwrite(rawoutbuf+dbps*nch*delay,1,dbps*nch*(nsmplwrt2-delay),fpo);
	      sumwrite += nsmplwrt2-delay;
	      init = 0;
	    }
	  }
	}

	{
	  int ds = (rp2-1)/(fs2/fs1);

	  if (ds > n1b2) ds = n1b2;

	  for(ch=0;ch<nch;ch++)
	    memmove(buf2[ch],buf2[ch]+ds,sizeof(REAL)*(n2x+1+n1b2-ds));

	  rp2 -= ds*(fs2/fs1);
	}

	for(ch=0;ch<nch;ch++)
	  memcpy(buf2[ch]+n2x+1,buf1[ch]+n1b2,sizeof(REAL)*n1b2);

	if ((spcount++ & 7) == 7) showprogress((double)sumread / chunklen);
      }
  }

  showprogress(1);

  SleefDFT_free(stage1);
  SleefDFT_dispose(dftf);
  SleefDFT_dispose(dftb);
  free(f2order);
  free(f2inc);
  free(stage2[0]);
  free(stage2);
  for(i=0;i<nch;i++) SleefDFT_free(buf1[i]);
  free(buf1);
  for(i=0;i<nch;i++) free(buf2[i]);
  free(buf2);
  free(inbuf);
  free(outbuf);
  free(rawinbuf);
  free(rawoutbuf);

  if (dither != -1) {
    static int64_t dlmin[] = { 0, -0x80, -0x8000, -0x800000, -0x80000000ULL };
    static int64_t dlmax[] = { 0,  0x7f,  0x7fff,  0x7fffff,  0x7fffffffULL };
    double min = 0, max = 0;
    for(i=0;i<nch;i++) {
      double p[2];
      SSRCDither_getPeaks(shaper[i], p);
      if (p[1] < min) min = p[1];
      if (p[0] > max) max = p[0];
    }
    if (min < dlmin[dbps]) peak = min / dlmin[dbps];
    if (max > dlmax[dbps]) peak = max / dlmax[dbps];
  }

  return peak;
}

double no_src(FILE *fpi,FILE *fpo,int nch,int bps,int dbps,double gain,int chunklen,int twopass,int dither)
{
  double peak=0;
  int ch=0,sumread=0;

  setstarttime();

  while(sumread < chunklen*nch)
    {
      REAL f;
      int s;
      uint8_t buf[4];

      switch(bps) {
      case 1:
	ignoreReturnValue(fread(buf,1,1,fpi));
	f = (1 / (REAL)0x7f) * ((REAL)((uint8_t *)buf)[0]-128);
	break;
      case 2:
	ignoreReturnValue(fread(buf,2,1,fpi));
#if 0
	f = (1 / (REAL)0x7fff) * ((REAL)(((int16_t *)buf)[0]);
#else
	f = (1 / (REAL)0x7fff) * (((int)buf[0]) | (((int)(((int8_t *)buf)[1])) << 8));
#endif
	break;
      case 3:
	ignoreReturnValue(fread(buf,3,1,fpi));
	f = (1 / (REAL)0x7fffff) * 
	  ((((int)buf[0]) << 0 ) |
	   (((int)buf[1]) << 8 ) |
	   (((int)((int8_t *)buf)[2]) << 16));
	break;
      case 4:
	ignoreReturnValue(fread(buf,4,1,fpi));
	f = (1 / (REAL)0x7fffffff) * 
	  ((((int)buf[0]) << 0 ) |
	   (((int)buf[1]) << 8 ) |
	   (((int)buf[2]) << 16) |
	   (((int)((int8_t *)buf)[3]) << 24));
	break;
      };

      if (feof(fpi)) break;
      f *= gain;

      if (!twopass) {
	switch(dbps) {
	case 1:
	  f *= 0x7f;
	  s = dither != -1 ? do_shaping(f,&peak,dither,ch) : RINT(f);
	  buf[0] = s + 128;
	  fwrite(buf,sizeof(int8_t),1,fpo);
	  break;
	case 2:
	  f *= 0x7fff;
	  s = dither != -1 ? do_shaping(f,&peak,dither,ch) : RINT(f);
	  buf[0] = s & 255; s >>= 8;
	  buf[1] = s & 255;
	  fwrite(buf,sizeof(int8_t),2,fpo);
	  break;
	case 3:
	  f *= 0x7fffff;
	  s = dither != -1 ? do_shaping(f,&peak,dither,ch) : RINT(f);
	  buf[0] = s & 255; s >>= 8;
	  buf[1] = s & 255; s >>= 8;
	  buf[2] = s & 255;
	  fwrite(buf,sizeof(int8_t),3,fpo);
	  break;
	};
      } else {
	REAL p = f > 0 ? f : -f;
	peak = peak < p ? p : peak;
	fwrite(&f,sizeof(REAL),1,fpo);
      }

      ch++;
      if (ch == nch) ch = 0;
      sumread++;

      if ((sumread & 0x3ffff) == 0) showprogress((double)sumread / (chunklen*nch));
    }

  showprogress(1);

  if (dither != -1) {
    static int64_t dlmin[] = { 0, -0x80, -0x8000, -0x800000, -0x80000000ULL };
    static int64_t dlmax[] = { 0,  0x7f,  0x7fff,  0x7fffff,  0x7fffffffULL };
    double min = 0, max = 0;
    int i;
    for(i=0;i<nch;i++) {
      double p[2];
      SSRCDither_getPeaks(shaper[i], p);
      if (p[1] < min) min = p[1];
      if (p[0] > max) max = p[0];
    }
    if (min < dlmin[dbps]) peak = min / dlmin[dbps];
    if (max > dlmax[dbps]) peak = max / dlmax[dbps];
  }

  return peak;
}

int extract_int(uint8_t *buf)
{
#ifndef BIGENDIAN
  return *(int *)buf;
#else
  return ((int)buf[0]) | (((int)buf[1]) << 8) | 
    (((int)buf[2]) << 16) | (((int)((int8_t *)buf)[3]) << 24);
#endif
}

unsigned int extract_uint(uint8_t *buf)
{
#ifndef BIGENDIAN
  return *(int *)buf;
#else
  return ((unsigned int)buf[0]) | (((unsigned int)buf[1]) << 8) | 
    (((unsigned int)buf[2]) << 16) | (((unsigned int)((int8_t *)buf)[3]) << 24);
#endif
}

int16_t extract_short(uint8_t *buf)
{
#ifndef BIGENDIAN
  return *(int16_t *)buf;
#else
  return ((int16_t)buf[0]) | (((int16_t)((int8_t *)buf)[1]) << 8);
#endif
}

void bury_int(uint8_t *buf,int i)
{
#ifndef BIGENDIAN
  *(int *)buf = i;
#else
  buf[0] = i & 0xff; i >>= 8;
  buf[1] = i & 0xff; i >>= 8;
  buf[2] = i & 0xff; i >>= 8;
  buf[3] = i & 0xff;
#endif
}

void bury_uint(uint8_t *buf,unsigned int i)
{
#ifndef BIGENDIAN
  *(int *)buf = i;
#else
  buf[0] = i & 0xff; i >>= 8;
  buf[1] = i & 0xff; i >>= 8;
  buf[2] = i & 0xff; i >>= 8;
  buf[3] = i & 0xff;
#endif
}

void bury_short(uint8_t *buf,int16_t s)
{
#ifndef BIGENDIAN
  *(int16_t *)buf = s;
#else
  buf[0] = s & 0xff; s >>= 8;
  buf[1] = s & 0xff;
#endif
}

int fread_int(FILE *fp)
{
#ifndef BIGENDIAN
  int ret;
  ignoreReturnValue(fread(&ret,4,1,fp));

  return ret;
#else
  uint8_t buf[4];
  ignoreReturnValue(fread(&buf,1,4,fp));
  return extract_int(buf);
#endif
}

unsigned int fread_uint(FILE *fp)
{
#ifndef BIGENDIAN
  unsigned int ret;
  ignoreReturnValue(fread(&ret,4,1,fp));

  return ret;
#else
  uint8_t buf[4];
  ignoreReturnValue(fread(&buf,1,4,fp));
  return extract_uint(buf);
#endif
}

int fread_short(FILE *fp)
{
#ifndef BIGENDIAN
  int16_t ret;
  ignoreReturnValue(fread(&ret,2,1,fp));
  return ret;
#else
  uint8_t buf[2];
  ignoreReturnValue(fread(&buf,1,2,fp));
  return extract_short(buf);
#endif
}

void fwrite_int(FILE *fp,int i)
{
#ifndef BIGENDIAN
  fwrite(&i,4,1,fp);
#else
  uint8_t buf[4];
  bury_int(buf,i);
  fwrite(&buf,1,4,fp);
#endif
}

void fwrite_uint(FILE *fp, uint32_t i)
{
#ifndef BIGENDIAN
  fwrite(&i,4,1,fp);
#else
  uint8_t buf[4];
  bury_uint(buf,i);
  fwrite(&buf,1,4,fp);
#endif
}

void fwrite_short(FILE *fp, int16_t s)
{
#ifndef BIGENDIAN
  fwrite(&s,2,1,fp);
#else
  uint8_t buf[4];
  bury_int(buf,s);
  fwrite(&buf,1,2,fp);
#endif
}

int main(int argc, char **argv)
{
  char *sfn,*dfn,*tmpfn=NULL;
  FILE *fpi = NULL,*fpo = NULL,*fpt = NULL;
  int twopass,normalize,dither,pdf,samp=0;
  int nch,bps;
  unsigned int length;
  int sfrq,dfrq,dbps;
  double att,peak,noiseamp;
  int i, j, k;

  // parse command line options

  dfrq = -1;
  att = 0;
  dbps = -1;
  twopass = 0;
  normalize = 0;
  dither = -1;
  pdf = 1;
  noiseamp = 1;

  for(i=1;i<argc;i++)
    {
      if (argv[i][0] != '-') break;

      if (strcmp(argv[i],"--rate") == 0) {
	dfrq = atoi(argv[++i]);
	continue;
      }

      if (strcmp(argv[i],"--att") == 0) {
	att = atof(argv[++i]);
	continue;
      }

      if (strcmp(argv[i],"--bits") == 0) {
	dbps = atoi(argv[++i]);
	if (dbps != 8 && dbps != 16 && dbps != 24) {
	  fprintf(stderr,"Error: Only 8bit, 16bit and 24bit PCM are supported.\n");
	  exit(-1);
	}
	dbps /= 8;
	continue;
      }

      if (strcmp(argv[i],"--twopass") == 0) {
	twopass = 1;
	continue;
      }

      if (strcmp(argv[i],"--normalize") == 0) {
	twopass = 1;
	normalize = 1;
	continue;
      }

      if (strcmp(argv[i],"--dither") == 0) {
	if (strcmp(argv[i+1], "help") == 0) {
	  printf("\nAll available noise shaper IDs:\n");
	  const int *fslist = SSRCDither_getAllSupportedFS();
	  for(j=0;fslist[j] != -1;j++) {
	    int ns = SSRCDither_getNumAvailableShapers(fslist[j]);
	    int *nsids = calloc(ns, sizeof(int));
	    SSRCDither_getAvailableShaperIDs(nsids, ns, fslist[j]);
	    for(k=0;k<ns;k++) {
	      printf("fs = %d, ID = %d : %s\n", fslist[j], nsids[k], SSRCDither_getNameForShaperID(fslist[j], nsids[k]));
	    }
	  }
	  exit(-1);
	}
	char *endptr;
	dither = strtol(argv[i+1],&endptr,10);
	if (*endptr == '\0') {
	  i++;
	} else {
	  fprintf(stderr,"Error: Unrecognized dither type %s\n", argv[i+1]);
	  exit(-1);
	}
	continue;
      }

      if (strcmp(argv[i],"--pdf") == 0) {
	char *endptr;
	pdf = strtol(argv[i+1],&endptr,10);
	if (*endptr == '\0') {
	  if (pdf != 0 && pdf != 1 && pdf != 3) {
	    fprintf(stderr,"unrecognized p.d.f. type : %s\n",argv[i+1]);
	    exit(-1);
	  }
	  i++;
	} else {
	  fprintf(stderr,"unrecognized p.d.f. type : %s\n",argv[i+1]);
	  exit(-1);
	}

	noiseamp = strtod(argv[i+1],&endptr);
	if (*endptr == '\0') {
	  i++;
	} else {
	  static double presets[] = {0.7, 0.9, 0, 0.5};
	  noiseamp = presets[pdf];
	}

	continue;
      }

      if (strcmp(argv[i],"--quiet") == 0) {
	quiet = 1;
	continue;
      }

      if (strcmp(argv[i],"--tmpfile") == 0) {
	tmpfn = argv[++i];
	continue;
      }

      if (strcmp(argv[i],"--profile") == 0) {
	if (strcmp(argv[i+1],"short") == 0) {
	  AA = 96;
	  DF = 2000;
	} else if (strcmp(argv[i+1],"normal") == 0) {
	  AA = 140;
	  DF = 2000;
	} else if (strcmp(argv[i+1],"long") == 0) {
	  AA = 150;
	  DF = 200;
	} else {
	  fprintf(stderr,"unrecognized profile : %s\n",argv[i+1]);
	  exit(-1);
	}
	i++;
	continue;
      }

      fprintf(stderr,"unrecognized option : %s\n",argv[i]);
      exit(-1);
    }

#ifndef HIGHPREC
  if (!quiet) printf("Shibatch sampling rate converter version " VERSION "(single precision)\n\n");
#else
  if (!quiet) printf("Shibatch sampling rate converter version " VERSION "(double precision)\n\n");
#endif

  if (argc-i != 2) {usage(); exit(-1);}

  sfn = argv[i];
  dfn = argv[i+1];

  fpi = fopen(sfn,"rb");

  if (!fpi) {fprintf(stderr,"cannot open input file.\n"); exit(-1);}

  /* read wav header */

  {
    uint8_t ibuf[576*2*2];
    int16_t word;
    int fmtheaderlen;
    uint16_t formattag;

    if (getc(fpi) != 'R') fmterr(1);
    if (getc(fpi) != 'I') fmterr(1);
    if (getc(fpi) != 'F') fmterr(1);
    if (getc(fpi) != 'F') fmterr(1);

    ignoreReturnValue(fread(ibuf,1,4,fpi));

    if (getc(fpi) != 'W') fmterr(2);
    if (getc(fpi) != 'A') fmterr(2);
    if (getc(fpi) != 'V') fmterr(2);
    if (getc(fpi) != 'E') fmterr(2);
    if (getc(fpi) != 'f') fmterr(2);
    if (getc(fpi) != 'm') fmterr(2);
    if (getc(fpi) != 't') fmterr(2);
    if (getc(fpi) != ' ') fmterr(2);

    fmtheaderlen = fread_int(fpi);
    if (fmtheaderlen < 16 || fmtheaderlen > sizeof ibuf) fmterr(3);
    ignoreReturnValue(fread(ibuf,fmtheaderlen,1,fpi));

    formattag = extract_short(&ibuf[0]);
    if (formattag != 1 && !(formattag == 0xFFFE && extract_short(&ibuf[24]) == 1)) {
      fprintf(stderr,"Error: Only PCM is supported.\n");
      exit(-1);
    }
    nch = extract_short(&ibuf[2]);
    sfrq = extract_int(&ibuf[4]);
    bps = extract_int(&ibuf[8]);
    if ((int)bps % sfrq*nch != 0) fmterr(4);

    bps /= sfrq*nch;

    for(;;)
      {
	int8_t buf[4];
	buf[0] = getc(fpi);
	buf[1] = getc(fpi);
	buf[2] = getc(fpi);
	buf[3] = getc(fpi);
	length = fread_uint(fpi);
	if (buf[0] == 'd' && buf[1] == 'a' && buf[2] == 't' && buf[3] == 'a') break;
	if (feof(fpi)) break;
	fseek(fpi, length, SEEK_CUR);
      }
    if (feof(fpi)) {
      fprintf(stderr,"Couldn't find data chunk\n");
      exit(-1);
    }
  }

  if (bps != 1 && bps != 2 && bps != 3 && bps != 4) {
    fprintf(stderr,"Error : Only 8bit, 16bit, 24bit and 32bit PCM are supported.\n");
    exit(-1);
  }

  if (dbps == -1) {
    if (bps != 1) dbps = bps;
    else dbps = 2;
    if (dbps == 4) dbps = 3;
  }

  if (dfrq == -1) dfrq = sfrq;

  if (dither != -1 && SSRCDither_getNameForShaperID(dfrq, dither) == NULL) {
    fprintf(stderr, "Noise shaper %d is not available.\n", dither);
    dither = -1;
  }

  if (sfrq < dfrq && DF < sfrq / 24.0) DF = sfrq/24.0;
  if (dfrq < sfrq && DF < dfrq / 24.0) DF = dfrq/24.0;

  if (!quiet) {
    const char *ptype[] = {
      "rectangular", "triangular", "gaussian", "two-level"
    };
    printf("frequency : %d -> %d\n",sfrq,dfrq);
    printf("attenuation : %gdB\n",att);
    printf("bits per sample : %d -> %d\n",bps*8,dbps*8);
    printf("nchannels : %d\n",nch);
    printf("length : %d bytes, %g secs\n",length,(double)length/bps/nch/sfrq);
    if (dither == -1) {
      printf("dither type : none\n");
    } else {
      printf("dither type : %s, %s p.d.f, amp = %g\n", SSRCDither_getNameForShaperID(dfrq, dither), ptype[pdf], noiseamp);
    }
    printf("\n");
  }

  if (twopass) {
    if (tmpfn) {
      fpt = fopen(tmpfn,"w+b");
    } else {
      fpt = tmpfile();
    }
    if (!fpt) {fprintf(stderr,"cannot open temporary file.\n"); exit(-1);}
  }

  fpo = fopen(dfn,"wb");

  if (!fpo) {fprintf(stderr,"cannot open output file.\n"); exit(-1);}

  /* generate wav header */

  {
    int16_t word;
    int dword;

    fwrite("RIFF",4,1,fpo);
    dword = 0;
    fwrite_int(fpo,dword);

    fwrite("WAVEfmt ",8,1,fpo);
    dword = 16;
    fwrite_int(fpo,dword);
    word = 1;
    fwrite_short(fpo,word); /* format category, PCM */
    word = nch;
    fwrite_short(fpo,word); /* channels */
    dword = dfrq;
    fwrite_int(fpo,dword); /* sampling rate */
    dword = dfrq*nch*dbps;
    fwrite_int(fpo,dword); /* bytes per sec */
    word = dbps*nch;
    fwrite_short(fpo,word); /* block alignment */
    word = dbps*8;
    fwrite_short(fpo,word); /* bits per sample */

    fwrite("data",4,1,fpo);
    dword = 0;
    fwrite_int(fpo,dword);
  }

  {
    int i,j;

    for(i=0;i<=M;i++)
      {
	fact[i] = 1;
	for(j=1;j<=i;j++) fact[i] *= j;
      }
  }

  if (dither != -1) {
    int min,max;
    if (dbps == 1) {min = -0x80; max = 0x7f;}
    if (dbps == 2) {min = -0x8000; max = 0x7fff;}
    if (dbps == 3) {min = -0x800000; max = 0x7fffff;}
    if (dbps == 4) {min = -0x80000000; max = 0x7fffffff;}

    for(i=0;i<nch;i++) {
      shaper[i] = SSRCDither_init(dfrq, min, max, dither, pdf, noiseamp, 0);
      assert(shaper[i] != NULL);
    }
  }

  if (twopass) {
    REAL gain = 1.0;
    int ch=0;
    unsigned int fptlen, sumread;

    if (!quiet) printf("Pass 1\n");

    if (normalize) {
      if (sfrq < dfrq) peak = upsample(fpi,fpt,nch,bps,sizeof(REAL),sfrq,dfrq,1,length/bps/nch,twopass,dither);
      else if (sfrq > dfrq) peak = downsample(fpi,fpt,nch,bps,sizeof(REAL),sfrq,dfrq,1,length/bps/nch,twopass,dither);
      else peak = no_src(fpi,fpt,nch,bps,sizeof(REAL),1,length/bps/nch,twopass,dither);
    } else {
      if (sfrq < dfrq) peak = upsample(fpi,fpt,nch,bps,sizeof(REAL),sfrq,dfrq,pow(10,-att/20),length/bps/nch,twopass,dither);
      else if (sfrq > dfrq) peak = downsample(fpi,fpt,nch,bps,sizeof(REAL),sfrq,dfrq,pow(10,-att/20),length/bps/nch,twopass,dither);
      else peak = no_src(fpi,fpt,nch,bps,sizeof(REAL),pow(10,-att/20),length/bps/nch,twopass,dither);
    }

    if (!quiet) printf("\npeak : %gdB\n",20*log10(peak));

    if (!normalize) {
      if (peak < pow(10,-att/20)) peak = 1;
      else peak *= pow(10,att/20);
    } else peak *= pow(10,att/20);

    if (!quiet) printf("\nPass 2\n");
    
    if (dither != -1) {
      switch(dbps)
	{
	case 1:
	  gain = (normalize || peak >= (0x7f-samp)/(double)0x7f) ? 1/peak*(0x7f-samp) : 1/peak*0x7f;
	  break;
	case 2:
	  gain = (normalize || peak >= (0x7fff-samp)/(double)0x7fff) ? 1/peak*(0x7fff-samp) : 1/peak*0x7fff;
	  break;
	case 3:
	  gain = (normalize || peak >= (0x7fffff-samp)/(double)0x7fffff) ? 1/peak*(0x7fffff-samp) : 1/peak*0x7fffff;
	  break;
	}
    } else {
      switch(dbps)
	{
	case 1:
	  gain = 1/peak * 0x7f;
	  break;
	case 2:
	  gain = 1/peak * 0x7fff;
	  break;
	case 3:
	  gain = 1/peak * 0x7fffff;
	  break;
	}
    }

    setstarttime();

    fptlen = ftell(fpt) / sizeof(REAL);
    sumread = 0;

    fseek(fpt,0,SEEK_SET);
    for(;;)
      {
	REAL f;
	int s;

	if (fread(&f,sizeof(REAL),1,fpt) == 0) break;
	f *= gain;
	sumread++;

	switch(dbps) {
	case 1:
	  {
	    uint8_t buf[1];
	    s = dither != -1 ? do_shaping(f,&peak,dither,ch) : RINT(f);

	    buf[0] = s + 128;

	    fwrite(buf,sizeof(int8_t),1,fpo);
	  }
	  break;
	case 2:
	  {
	    int8_t buf[2];
	    s = dither != -1 ? do_shaping(f,&peak,dither,ch) : RINT(f);

	    buf[0] = s & 255; s >>= 8;
	    buf[1] = s & 255;

	    fwrite(buf,sizeof(int8_t),2,fpo);
	  }
	  break;
	case 3:
	  {
	    int8_t buf[3];
	    s = dither != -1 ? do_shaping(f,&peak,dither,ch) : RINT(f);

	    buf[0] = s & 255; s >>= 8;
	    buf[1] = s & 255; s >>= 8;
	    buf[2] = s & 255;

	    fwrite(buf,sizeof(int8_t),3,fpo);
	  }
	  break;
	}

	ch++;
	if (ch == nch) ch = 0;

	if ((sumread & 0x3ffff) == 0) showprogress((double)sumread / fptlen);
      }
    showprogress(1);
    if (!quiet) printf("\n");
    fclose(fpt);
    if (tmpfn != NULL) {
      if (remove(tmpfn))
	fprintf(stderr,"Failed to remove %s\n",tmpfn);
    }
  } else {
    if (sfrq < dfrq) peak = upsample(fpi,fpo,nch,bps,dbps,sfrq,dfrq,pow(10,-att/20),length/bps/nch,twopass,dither);
    else if (sfrq > dfrq) peak = downsample(fpi,fpo,nch,bps,dbps,sfrq,dfrq,pow(10,-att/20),length/bps/nch,twopass,dither);
    else peak = no_src(fpi,fpo,nch,bps,dbps,pow(10,-att/20),length/bps/nch,twopass,dither);
    if (!quiet) printf("\n");
  }

  if (dither != -1) {
    for(i=0;i<nch;i++) SSRCDither_dispose(shaper[i]);
  }

  if (!twopass && peak > 1) {
    if (!quiet) printf("clipping detected : %gdB\n",20*log10(peak));
  }

  {
    int16_t word;
    int dword;
    int len;

    fseek(fpo,0,SEEK_END);
    len = ftell(fpo);

    fseek(fpo,4,SEEK_SET);
    dword = len-8;
    fwrite_int(fpo,dword);

    fseek(fpo,40,SEEK_SET);
    dword = len-44;
    fwrite_int(fpo,dword);
  }

  exit(0);
}
