#include "ViennaRNA/utils.hh"
#include "ViennaRNA/array.hh"
#include <ctime>
#include <unistd.h>

/**
 *  @brief Define an array
 */
#define vrna_array(Type) Type *

typedef struct {
  vrna_log_cb_f       cb;
  void                *cb_data;
  vrna_logdata_free_f data_release;
  vrna_log_levels_e   level;
} logger_callback;

static struct {
  FILE              *default_file;
  vrna_log_levels_e default_level;
  unsigned int      options;
  vrna_log_lock_f   lock;
  void              *lock_data;
  vrna_array(logger_callback) callbacks;
#if VRNA_WITH_PTHREADS
  pthread_mutex_t   mtx;            /* semaphore to prevent concurrent access */
#endif
} logger = {
  .default_file   = NULL,
  .default_level  = VRNA_LOG_LEVEL_DEFAULT,
  .options        = VRNA_LOG_OPTION_DEFAULT,
  .lock           = NULL,
  .lock_data      = NULL,
  .callbacks      = NULL,
#if VRNA_WITH_PTHREADS
  .mtx            = PTHREAD_MUTEX_INITIALIZER
#endif
};

/*
 #################################
 # PRIVATE FUNCTION DECLARATIONS #
 #################################
 */
static void log_v(vrna_log_event_t *event);

/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
static const char* get_log_level_string(vrna_log_levels_e level)
{
  switch (level) {
    case VRNA_LOG_LEVEL_DEBUG:
      return "[DEBUG]";
    case VRNA_LOG_LEVEL_INFO:
      return "[INFO]";
    case VRNA_LOG_LEVEL_WARNING:
      return "[WARNING]";
    case VRNA_LOG_LEVEL_ERROR:
      return "[ERROR]";
    case VRNA_LOG_LEVEL_CRITICAL:
      return "[FATAL]";
    default:
      return "[UNKNOWN]";
  }
}
static const char* get_log_level_color(vrna_log_levels_e level)
{
  switch (level) {
    case VRNA_LOG_LEVEL_DEBUG:
      return ANSI_COLOR_CYAN_B;
    case VRNA_LOG_LEVEL_INFO:
      return ANSI_COLOR_BLUE_B;
    case VRNA_LOG_LEVEL_WARNING:
      return ANSI_COLOR_YELLOW_B;
    case VRNA_LOG_LEVEL_ERROR:
      return ANSI_COLOR_RED_B;
    case VRNA_LOG_LEVEL_CRITICAL:
      return ANSI_COLOR_MAGENTA_B;
    default:
      return ANSI_COLOR_RESET;
  }
}

static void
lock(void);


static void
unlock(void);

void
vrna_log(vrna_log_levels_e  level,
         const char         *file_name,
         int                line_number,
         const char         *format_string,
         ...)
{
  vrna_log_event_t event{};
  event.format_string = format_string;
  event.level         = level;
  event.line_number   = line_number;
  event.file_name     = file_name;

  va_start(event.params, format_string);
  log_v(&event);
  va_end(event.params);
}

/*
 #################################
 # STATIC helper functions below #
 #################################
 */
static void log_default(vrna_log_event_t *event){
  if (!logger.default_file)
    logger.default_file = stderr;

  /* print time unless turned off explicitely */
  if (logger.options & VRNA_LOG_OPTION_TRACE_TIME) {
    char    timebuf[64];
    time_t  t = time(NULL);
    timebuf[strftime(timebuf, sizeof(timebuf), "%H:%M:%S", localtime(&t))] = '\0';
    fprintf(logger.default_file, "%s ", timebuf);
  }

  /* print log level */
  if (isatty(fileno(logger.default_file))) {
    fprintf(logger.default_file,
            "%s%-9s" ANSI_COLOR_RESET " ",
            get_log_level_color(event->level),
            get_log_level_string(event->level));
  }

  /* print file name / line number trace unless turned off explicitely */
  if (logger.options & VRNA_LOG_OPTION_TRACE_CALL) {
    if (isatty(fileno(logger.default_file))) {
      fprintf(logger.default_file,
              "\x1b[90m%s:%d:" ANSI_COLOR_RESET " ",
              event->file_name,
              event->line_number);
    }
  }

  /* print actual message */
  vfprintf(logger.default_file, event->format_string, event->params);
  fprintf(logger.default_file, "\n");
  fflush(logger.default_file);
}

static void
log_v(vrna_log_event_t *event)
{
  lock();

  /* initialize the logger, if not done already */
  if (logger.callbacks == NULL) {
    /* initialize callbacks if not done so far */
    vrna_array_init(logger.callbacks);
  }

  /* process log output for default implementation */
  if (!(logger.options & VRNA_LOG_OPTION_QUIET)) {
    /* print log if not in quiet mode */
    if (event->level >= logger.default_level)
      log_default(event);
  }

  /* process log for any user-defined output */
  for (size_t i = 0; i < vrna_array_size(logger.callbacks); i++) {
    logger_callback *cb = &(logger.callbacks[i]);

    if (event->level >= cb->level)
      cb->cb(event, cb->cb_data);
  }

  unlock();
}


static void
lock(void)
{
  if (logger.lock)
#if VRNA_WITH_PTHREADS
  {
#endif
    logger.lock(1, logger.lock_data);

#if VRNA_WITH_PTHREADS
} else {
  pthread_mutex_lock(&(logger.mtx));
}
#endif
}


static void
unlock(void)
{
  if (logger.lock)
#if VRNA_WITH_PTHREADS
  {
#endif
    logger.lock(0, logger.lock_data);

#if VRNA_WITH_PTHREADS
} else {
  pthread_mutex_unlock(&(logger.mtx));
}
#endif
}

/* uniform random number generator; vrna_urn() is in [0,1] */
/* uses a linear congruential library routine */
/* 48 bit arithmetic */
double vrna_urn(void){
#ifdef HAVE_ERAND48
  extern double erand48(unsigned short[]);


  return erand48(xsubi);
#else
  return ((double)rand()) / RAND_MAX;
#endif
}