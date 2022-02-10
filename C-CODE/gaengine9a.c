/*------------------------------------------
 Function GAENGINE_FITTING.C
 ===================
 Simulation of the SWIFT model along experimental data
 --------------------------------------------*/

// Set dummy variable dimensions
// Dummy variables are additional columns in corpus and fixseq files
// They can be integers (idum), doubles (ddum) or strings (cdum)
// The number of dummy variables is defined in the following constants
// The order of variables in corpus/fixseq files is always:
// corpus: sentence length | word length | word freq | word predictability | idum[1..n] | ddum[1..n] | cdum[1..n]
// fixseq: sentence | fixated word | fixated letter | fixation duration | saccade duration | idum[1..n] | ddum[1..n] | cdum[1..n]
#define N_FIXSEQ_IDUM 0
#define N_FIXSEQ_DDUM 0
#define N_FIXSEQ_CDUM 0
#define N_CORPUS_IDUM 0
#define N_CORPUS_DDUM 0
#define N_CORPUS_CDUM 2
#define N_LOGLIKS 3

// What is the version of this algorithm?
#define SWIFT_VERSION_MAJOR 14
#define SWIFT_VERSION_MINOR 1

#define SWIFT_VARIANT "actr1-lrp"


// Numerical recipes includes
//#include "./NumericalRecipes/nrutil.c"
#include "logsumexp.c"


#include "gausslike.c"
#include "epallike.c"


#if defined(LEXRATE_GAUSS)
#include "lexrate7_gauss.c" // gaussian span
#elif defined(LEXRATE_INVPARAB)
#include "lexrate8_invparab.c" // gaussian span
#else
#warning You have not specified a lexrate module. Using LEXRATE_INVPARAB by default.
#include "lexrate8_invparab.c" // inverse parabolic span
#endif

#if defined(EXECSACC_GAUSS)
#include "execsacc_gauss.c"
#elif defined(EXECSACC_INVGAUSS)
#include "execsacc_invgauss.c"
#elif defined(EXECSACC_INVGAMMA)
#include "execsacc_invgamma.c"
#elif defined(EXECSACC_GAMMA)
#include "execsacc_gamma.c"
#else
#warning You have not specified a saccade execution module. Using EXECSACC_INVGAMMA by default.
#include "execsacc_gamma.c"
#endif

#include "selectar3a.c" // distance-dependent target selection
//#include "selectar3b.c" // distance-independent target selection
//#include "selectar3d.c" // distance-independent target selection with improved fallback

#include "logliktimer4.c"



// This is a shorthand version to require that specific parameters have been set
#define require(params, par) require_parameter(params, swift_parameter_ids.par)

