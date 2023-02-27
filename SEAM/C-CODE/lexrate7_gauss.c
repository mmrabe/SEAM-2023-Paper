/*-----------------------------------------------------------
 Iteration function for lexical processing
 (Version 3.0, 08-JUL-09)
 (Version 4.0, 19-FEB-11)
 (Version 5.0, 15-JAN-17)
 (Version 5.1, 28-AUG-18)
 -------------------------------------------------------------
 parameters:
 NW            number of words in current sentence
 a[1..NW]      vector of lexical activities
 as[1..NW]     lexical processing state (1/2/3)
 delta         parameter of lexical processing rates
 ppf           global inhibition
 k             current fixation position
 sacc          saccade in execution (0/1)
 t             time  
 */

#define sq(A) ((A)*(A))
#define cu(A) ((A)*(A)*(A))
#define realmin 1.0e-10
#define SWIFT_VARIANT_LEXRATE "gauss"

void lexrate(double *a,int *as,double *view,int *len,double aord,double alpha,double delta0,double delta1,
             double asym,double *pre,double proc,double ppf,double theta,double decay,double eta,
             int NW,double kpos,int k,int sacc,double *procrate,int fitting)

{
    double     dynfac, deltar, deltal, prate;
    double     ecc, leftact, hl, hr, normproc;
    int        j, l;
    
	/* dynamical processing span */ 
    dynfac = 1.0 - a[k];
	deltar = delta0 + delta1*dynfac;
	deltal = delta0*asym;
	
    /* normalization constant for the asymmetric Gaussian processing span */
    normproc = sqrt(3.1415926536/2.0)*(deltal+deltar);
    hl = 2.0*sq(deltal);
    hr = 2.0*sq(deltar);
    
    /* calculating lexical processing rates */
    for ( j=1; j<=NW; j++ ) {
        procrate[j] = 0.0;
        for ( l=1; l<=len[j]; l++ )  {
            ecc = view[j] - len[j]/2.0 + 1.0*l - 0.5 - kpos;  // modified by Ralf on Dec. 22, 2016
            if ( ecc<0.0 )
                prate = exp(-sq(ecc)/hl)/normproc;
            else
                prate = exp(-sq(ecc)/hr)/normproc;
            procrate[j] += prate;
        }
        procrate[j] *= exp(-eta*log(len[j]));
        
        /* computing activation to the left */
        //for ( l=1, leftact=1.0; l<j; l++ )  leftact *= exp(-ppf*a[l]/aord);
        for ( l=1, leftact=1.0; l<j; l++ )  leftact *= exp(-ppf*a[l]);
        procrate[j] *= leftact;
        
        /* predictability effect */
        procrate[j] *= exp(theta*sq(pre[j]));
    }
    
    /* adding decay if as[.]==2 */
    for ( j=1; j<=NW; j++ ) {
        if ( as[j]==2 )  {
            procrate[j] *= proc;
            if ( procrate[j]<decay )  procrate[j] = decay;
        }
        if ( sacc==1 && as[j]==1 )  procrate[j] = 0.0;
    }
}

#undef sq
#undef realmin
