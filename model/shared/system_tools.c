#include <stdio.h>
#include <sys/types.h>
#include <dirent.h>
#include <string.h>

int c_field_list(char *dir, char **list, int nl, int ls)
{
  DIR *d;
  struct dirent *entry;
  int n;

  printf("Got request: %s %d %d\n", dir, nl, ls);

  if ((d = opendir (dir)) != NULL) {
    n = 0;
    while ((entry = readdir (d)) != NULL) {
      if ( n >= nl ) return -3; /* list too long */
      /* printf ("%s\n", entry->d_name); */
      if ( strlen(entry->d_name) > ls-1 ) return -2; /* name too long */
      strncpy(list[n], entry->d_name, ls);
      n++;
    }
    closedir (d);
  } else {
    /* could not open directory */
    /* perror ("c_field_list: could not open directory"); */
    return -1;
  }

  return n;
}
