#define H__H__GCM_GOTO_VALID_LINE__H__H

#include <stdio.h>

/// move file pointer position to valid line. Skip line start with #
///
/// \param[in,out] in file pointer to be moved to a valid line
/// \param[in]     line_size maximum line size in reading a line from file (in), default is 1024
/// \return        non-zero on interal error
int goto_valid_line(FILE *in,
                    int line_size = -1);