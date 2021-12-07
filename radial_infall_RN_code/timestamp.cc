// timestamp.c exists solely to contain the VERSTRING and to be compiled
// every time any of the other source files changes
//
// File created on Jul 27, 2006
// Last Change: Thu Jul 27 11:04 2006
// Author: Roland Haas
// 

#include "header.h"

char VERSTRING[] = "$VER: " __FILE__ " compiled on " __DATE__ " at " __TIME__;
