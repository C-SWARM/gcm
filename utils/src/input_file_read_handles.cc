#include "input_file_read_handles.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/// read a line from a file pointer
int goto_valid_line(FILE *in,
                    int line_size)
{
  int err = 0;
  static const char delim[] = " \t\n";
  char L[1024];
  int N = 1024;

  char *line = L;
  
  if(line_size>0)
  {  
    line = (char *) malloc(line_size);
    N = line_size;
  }

  fpos_t pos;
  char *tok = NULL;

  do
  {
    err += fgetpos(in,&pos);
    if( fgets(line,N,in) == NULL)
    {
      err++;
      break;
    }
    // make sure got whole line (check is '\n')
    if (line[strlen(line)-1] != '\n' && !feof(in))
    {
      printf("ERROR: line too long (>%d chars)! %s(%s)\n", N, __func__, __FILE__);
      err++;
      break;
    }
        
    tok = strtok(line,delim);

    if (tok == NULL) tok = line + strlen(line);      

  }
  while(tok[0] == '#' || tok[0] == '\0');
  
  if(err==0)
    err += fsetpos(in,&pos);
  
  // free if line was created
  if(line_size>0)
    free(line);
  
  return err;
} 

#define _XX_XX_1 1
#ifndef _XX_XX_1
#define _XX_XX_1
int scan_for_valid_line(FILE *in)
{
  static const size_t line_length = 1024;
  static const char delim[] = " \t\n";

  int err = 0;
  char *line = malloc(line_length);
  char *tok = NULL;
  fpos_t pos;

  /* scan for non-comment/blank line */
  do{
    /* get the starting file position for the line */
    err += fgetpos(in,&pos);

    /* get a line and exit if there is an error */
    if ( fgets(line,line_length,in) == NULL) {
      err++;
      goto exit_err;
    }

    /* make sure got whole line (last char is '\n') */
    if ( line[strlen(line) - 1] != '\n' && !feof(in)) {
      fprintf(stderr,"ERROR: line too long (>%zd chars)! %s(%s)\n",
              line_length, __func__, __FILE__);
      err++;
      goto exit_err;
    }

    /* get first token */
    tok = strtok(line,delim);
    if (tok == NULL) tok = line + strlen(line);
  } while ( tok[0] == '#' || tok[0] == '\0');

  /* return the file pointer to the beginning of the valid line */
  err += fsetpos(in,&pos);

 exit_err:
  free(line);
  return err;
}

#endif