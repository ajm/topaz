#ifndef H_SSW_WRAPPER
#define H_SSW_WRAPPER

#include "seq.h"
#include "hit.h"
#include "options.h"

#include "ssw.h"


s_profile* get_ssw_profile(seq_t* s1, options_t* opt);
void align_using_ssw_profile(s_profile* profile, seq_t* s2, hit_t* hit, options_t* opt);
void align_using_ssw(seq_t* s1, seq_t* s2, hit_t* hit, options_t* opt);

#endif

