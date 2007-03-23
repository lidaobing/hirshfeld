// $Id: port.h 14 2005-03-09 17:09:10Z nichloas $
#ifndef PORT_H
#define PORT_H

#ifdef __cplusplus
extern "C" {
#endif

#if HAVE_CONFIG_H
# include <config.h>
#endif /* HAVE_CONFIG_H */

#if HAVE_LIBGEN_H
#  include <libgen.h>
#else
  char *basename(char *path);
#endif /* HAVE_LIBGEN_H */

#ifdef __cplusplus
}
#endif

#endif

