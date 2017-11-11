#ifndef H_GENERIC
#define H_GENERIC


#define PROGRAM "TOPAZ"

#define MIN(_a, _b) (((_a) < (_b)) ? (_a) : (_b))
#define MAX(_a, _b) (((_a) > (_b)) ? (_a) : (_b))

#define LCP_SIZE (16)
#define LCP_MAX ((1 << LCP_SIZE) - 1)

#define INDEX_COMMAND 0
#define SEARCH_COMMAND 1


#endif

