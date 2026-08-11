#ifndef UTILS
#define UTILS
#include <cstdarg>
#include <cstdio>
#include <cstdlib>

#define vrna_alloc(S)       calloc(1, (S))
#define vrna_realloc(p, S)  realloc(p, S)

#define ANSI_COLOR_BRIGHT     "\x1b[1m"
#define ANSI_COLOR_UNDERLINE  "\x1b[4m"
#define ANSI_COLOR_RED        "\x1b[31m"
#define ANSI_COLOR_GREEN      "\x1b[32m"
#define ANSI_COLOR_YELLOW     "\x1b[33m"
#define ANSI_COLOR_BLUE       "\x1b[34m"
#define ANSI_COLOR_MAGENTA    "\x1b[35m"
#define ANSI_COLOR_CYAN       "\x1b[36m"
#define ANSI_COLOR_RED_B      "\x1b[1;31m"
#define ANSI_COLOR_GREEN_B    "\x1b[1;32m"
#define ANSI_COLOR_YELLOW_B   "\x1b[1;33m"
#define ANSI_COLOR_BLUE_B     "\x1b[1;34m"
#define ANSI_COLOR_MAGENTA_B  "\x1b[1;35m"
#define ANSI_COLOR_CYAN_B     "\x1b[1;36m"
#define ANSI_COLOR_RESET      "\x1b[0m"

/**
 *  @brief  Log option to turn on call tracing
 *
 *  When this option is set via vrna_log_options_set()
 *  the internal logging system will include a call
 *  trace to the log output, i.e. the source code file
 *  and line numbers will be included in the log message.
 *
 *  @see  vrna_log_options_set(), vrna_log_options(), vrna_log_reset(),
 *        #VRNA_LOG_OPTION_QUIET, #VRNA_LOG_OPTION_TRACE_TIME,
 *        #VRNA_LOG_OPTION_DEFAULT
 */
#define VRNA_LOG_OPTION_TRACE_CALL  2U
/**
 *  @brief  Log option to turn on time stamp
 *
 *  When this option is set via vrna_log_options_set()
 *  the internal logging system will include a time stamp
 *  to the log output, i.e. the time when the log message
 *  was issued will be included in the log message.
 *
 *  @see  vrna_log_options_set(), vrna_log_options(), vrna_log_reset(),
 *        #VRNA_LOG_OPTION_QUIET, #VRNA_LOG_OPTION_TRACE_CALL,
 *        #VRNA_LOG_OPTION_DEFAULT
 */
#define VRNA_LOG_OPTION_TRACE_TIME  4U
/**
 *  @brief  Default log level
 *
 *  @see  vrna_log_level_set(), vrna_log_reset(),
 *        vrna_log_level(), #vrna_log_levels_e
 */
#define VRNA_LOG_LEVEL_DEFAULT      VRNA_LOG_LEVEL_ERROR

/**
 *  @brief  Log option representing the default options
 *
 *  When this option is set via vrna_log_options_set()
 *  the default options will be set.
 *
 *  @see  vrna_log_options_set(), vrna_log_options(), vrna_log_reset(),
 *        #VRNA_LOG_OPTION_QUIET, #VRNA_LOG_OPTION_TRACE_CALL,
 *        #VRNA_LOG_OPTION_TRACE_TIME
 */
#define VRNA_LOG_OPTION_DEFAULT     0U

/**
 *  @brief  Log option to turn off internal logging
 *
 *  When this option is set via vrna_log_options_set()
 *  the internal logging system will be deactivated and
 *  only user-defined callbacks will be seeing any logs.
 *
 *  @see  vrna_log_options_set(), vrna_log_options(), vrna_log_reset(),
 *        #VRNA_LOG_OPTION_TRACE_CALL, #VRNA_LOG_OPTION_TRACE_TIME,
 *        #VRNA_LOG_OPTION_DEFAULT
 */
#define VRNA_LOG_OPTION_QUIET       1U

/**
 *  @brief  The log levels
 */
typedef enum {
  VRNA_LOG_LEVEL_UNKNOWN  = -1,   /**< Unknown log level */
  VRNA_LOG_LEVEL_DEBUG    = 10,   /**< Debug log level */
  VRNA_LOG_LEVEL_INFO     = 20,   /**< Info log level */
  VRNA_LOG_LEVEL_WARNING  = 30,   /**< Warning log level */
  VRNA_LOG_LEVEL_ERROR    = 40,   /**< Error log level */
  VRNA_LOG_LEVEL_CRITICAL = 50,   /**< Critical log level */
  VRNA_LOG_LEVEL_SILENT   = 999   /**< Silent log level */
} vrna_log_levels_e;

typedef void (*vrna_logdata_free_f)(void *data);

/**
 *  @brief  The lock function prototype that may be passed to the logging system
 *
 *  @see  vrna_log_lock_set()
 *
 *  @param  lock      A parameter indicating whether to lock (lock != 0) or unlock (lock == 0)
 *  @param  lock_data An arbitrary user-defined data pointer for the user-defined locking system
 */
typedef void (*vrna_log_lock_f)(int   lock,
                                void  *lock_data);

/**
 *  @brief  A log event
 */
typedef struct vrna_log_event_s {
  const char        *format_string; /**< The printf-like format string containing the log information */
  va_list           params;         /**< The parameters for the printf-like format string */
  vrna_log_levels_e level;          /**< The log level */
  int               line_number;    /**< The source code line number that issued the log */
  const char        *file_name;     /**< The source code file that issued the log */
} vrna_log_event_t;


inline void vrna_message_warning(const char *format, ...){
    va_list args;

    va_start(args, format);
    fprintf(stderr, "WARNING: ");
    vfprintf(stderr, format, args);
    fprintf(stderr, "\n");
    va_end(args);
}

inline void vrna_message_error(const char *format,...){
    va_list args;

    va_start(args, format);
    fprintf(stderr, "ERROR: ");
    vfprintf(stderr, format, args);
    fprintf(stderr, "\n");
    va_end(args);
}

/**
 *  @brief  Issue a warning log message
 *
 *  This macro expects a printf-like format string followed by a variable list of
 *  arguments for the format string and passes this content to the log system.
 *
 *  @see  #vrna_log_debug, #vrna_log_info, #vrna_log_error, #vrna_log_critical,
 *        vrna_log(), vrna_log_level_set(), vrna_log_options_set(), vrna_log_fp_set()
 */
#define vrna_log_warning(...) \
        do { \
          vrna_log(VRNA_LOG_LEVEL_WARNING, __FILE__, __LINE__, __VA_ARGS__); \
        } while (0)

/**
 *  @brief  Issue an error log message
 *
 *  This macro expects a printf-like format string followed by a variable list of
 *  arguments for the format string and passes this content to the log system.
 *
 *  @see  #vrna_log_debug, #vrna_log_info, #vrna_log_warning, #vrna_log_critical,
 *        vrna_log(), vrna_log_level_set(), vrna_log_options_set(), vrna_log_fp_set()
 */
#define vrna_log_error(...) \
        do { \
          vrna_log(VRNA_LOG_LEVEL_ERROR, __FILE__, __LINE__, __VA_ARGS__); \
        } while (0)

/**
 *  @brief  Issue a log message
 *
 *  This is the low-level log message function. Usually, you don't want to call
 *  it directly but rather call one of the following high-level macros instead:
 *
 *  - #vrna_log_debug
 *  - #vrna_log_info
 *  - #vrna_log_warning
 *  - #vrna_log_error
 *  - #vrna_log_critical
 *
 *  @see  #vrna_log_debug, #vrna_log_info, #vrna_log_warning, #vrna_log_error,
 *        #vrna_log_critical, vrna_log_level_set(), vrna_log_options_set(),
 *        vrna_log_fp_set()
 *
 *  @param  level         The log level
 *  @param  file_name     The source code file name of the file that issued the log
 *  @param  line_number   The source code line number that issued the log
 *  @param  format_string The printf-like format string containing the log message
 *  @param  ...           The variable argument list for the printf-like @p format_string
 */
void
vrna_log(vrna_log_levels_e  level,
         const char         *file_name,
         int                line_number,
         const char         *format_string,
         ...);

/**
 *  @brief  The log callback function prototype
 *
 *  @see    vrna_log_cb_add(), vrna_log_cb_num(), vrna_log_cb_remove()
 *
 *  @param  event     The log event
 *  @param  log_data  An arbitrary user-defined data pointer for the user-define log message receiver
 */
typedef void (*vrna_log_cb_f)(vrna_log_event_t  *event,
                              void              *log_data);
#endif

/**
 *  @brief get a random number from [0..1]
 *
 *  @see  vrna_int_urn(), vrna_init_rand()
 *  @note Usually implemented by calling @e erand48().
 *  @return   A random number in range [0..1]
 */
double vrna_urn(void);