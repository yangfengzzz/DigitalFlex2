#ifndef __Version_h__
#define __Version_h__

#define STRINGIZE_HELPER(x) #x
#define STRINGIZE(x) STRINGIZE_HELPER(x)
#define WARNING(desc) message(__FILE__ "(" STRINGIZE(__LINE__) ") : Warning: " #desc)

#define GIT_SHA1 "619a911372dca5f48c562d39c13b9f83eda6e10a"
#define GIT_REFSPEC "refs/heads/master"
#define GIT_LOCAL_STATUS "DIRTY"

#define SPLISHSPLASH_VERSION "2.11.2"

#ifdef DL_OUTPUT
#pragma WARNING(Local changes not committed.)
#endif

#endif
