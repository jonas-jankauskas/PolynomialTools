/** @file		debug.h
	@author		Jonas Jankauskas
	@date		November 22, 2019
	@version	1.0
	@note		Flint version 2.5.2
	@brief		Simple C debuging system
	@details	Header file with macros that print various messages (usually, the content of various FLINT variables) when the program is compiled with -DDEBUG=lvl option, where lvl indicates the level of detail (0, 1, 2, ... etc). By default, these macros are disabled to keep the code as fast as possible. The code was inspired by answers in the discussion given in https://stackoverflow.com/questions/1644868/define-macro-for-debug-printing-in-c
*/

#ifndef DOXYGEN_SHOULD_SKIP_THIS

	#ifndef DEBUG_H
		
		#define DEBUG_H
		
		#include <stdarg.h>
	
		#include "flint/flint.h"
		#include "flint/fmpz.h"
		#include "flint/fmpq.h"
		#include "flint/fmpq_poly.h"
	
		#ifdef DEBUG
			#define DEBUG_PRINT_AT($lvl, fmt, ...) do { if (DEBUG>=$lvl) flint_printf("# %s:%d:%s():\n" fmt, __FILE__, __LINE__, __func__, ##__VA_ARGS__); } while (0)
			#define DEBUG_MSG_AT($lvl, fmt, ...) do { if (DEBUG>=$lvl) flint_printf(fmt, ##__VA_ARGS__); } while (0)
			#define DEBUG_FMPZ_AT($lvl, $var) do { if (DEBUG>=$lvl) fmpz_print($var); } while (0)
			#define DEBUG_FMPQ_AT($lvl, $var) do { if (DEBUG>=$lvl) fmpq_print($var); } while (0)
			#define DEBUG_FMPQ_POLY_AT($lvl, $var) do { if (DEBUG>=$lvl) fmpq_poly_print($var); } while (0)
			#define DEBUG_T($lvl, $n, $T, $sig) do { if (DEBUG>=$lvl) print_T_data($n, $T, $sig); } while (0)
			#define DEBUG_ENTER_AT($lvl) do { if (DEBUG>=$lvl) flint_printf("# %s:%d: %s(): entering...>\n", __FILE__, __LINE__, __FUNCTION__);  } while (0)
			#define DEBUG_LEAVE_AT($lvl) do { if (DEBUG>=$lvl) flint_printf("# %s:%d: %s(): ...leaving <\n",  __FILE__, __LINE__, __FUNCTION__);  } while (0)
		#else
			#define DEBUG_PRINT_AT($lvl, fmt, ...)
			#define DEBUG_MSG_AT($lvl, fmt, ...)
			#define DEBUG_FMPZ_AT($lvl, $var)
			#define DEBUG_FMPQ_AT($lvl, $var)
			#define DEBUG_FMPQ_POLY_AT($lvl, $var)
			#define DEBUG_T($lvl, $n, $T, $sig)
			#define DEBUG_ENTER_AT($lvl)
			#define DEBUG_LEAVE_AT($lvl)
		#endif

#endif

#endif /* DOXYGEN_SHOULD_SKIP_THIS */
	