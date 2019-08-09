#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

int main (int argc, char *argv[])
{
  char Buffer[4096];
  char LineID[256];
  char LineAC[256];
  const char DirectoryTypeM[] = "MATRIX";
  const char DirectoryTypeP[] = "PATTERN";
  
	if (argc != 2) {
		printf("Usage: %s <concatenated profiles file>\n"
		       "\t This will generate folders MATRIX and PATTERN, and populate them with single file profiles.\n"
		       "\t Typical use is with prosite.dat file\n",
		       argv[0]); 
	}
	
  FILE* in = fopen(argv[1], "r");
  if (in == NULL) return 0;
 
  /* Generating directories if required */
  {
      struct stat st;
      
      if ( stat("MATRIX", &st) != 0) {
	  if (mkdir("MATRIX", S_IRWXU) != 0) {
	      perror("Creating drirectory MATRIX");
	      exit(1);
	  }
      }
      if ( stat("PATTERN", &st) != 0) {
	  if (mkdir("PATTERN", S_IRWXU) != 0) {
	      perror("Creating drirectory PATTERN");
	      exit(1);
	  }
      }
  }
   
  size_t iline = 0;
  _Bool OpenProfile = false;
  FILE * newFile = NULL;
  char newFileName[128];
  const char * DirType = NULL; 
  while (1) {

      const char * Line = fgets(Buffer, 4096, in);
      iline++;
      if (Line == NULL) {
				if (feof(in)) 
					break;
				else {
					fprintf(stderr, "Error reading line %lu in %s\n", iline, argv[1]);
					fclose(in);
					return 1;
				}
      }
      
      /* Remove heading spaces */
      while (*Line == ' ' || *Line == '\t') Line++;
      
      /* What line ID to we have ? */
      if (Line[0] == '/' && Line[1] == '/') {
				// New profile ?
				if (OpenProfile) {
						fputs(Buffer, newFile);
						newFileName[0] = '\0';
						fclose(newFile);
						newFile = NULL;
						DirType = NULL;
						OpenProfile = false;
				}
      }
      else if (Line[0] == 'I' && Line[1] == 'D') {
				// MATRIX or PATTERN
				if (strstr(Buffer, "PATTERN") == NULL) {
						DirType = &DirectoryTypeM[0];
				}
				else {
						DirType = &DirectoryTypeP[0];
				}
				strncpy(LineID, Line, 256);
      }
      else if (Line[0] == 'A' && Line[1] == 'C') {
				if (DirType == NULL) {
					fprintf(stderr,"No ID line found before AC line %lu > %sLast ID was > %s\n", iline, Line, LineID);
					exit(1);
				}
				strncpy(LineAC, Line, 256);
				Line += 2;
				while (*Line == ' ' || *Line == '\t') {Line++; }
				int i = 0;
				while (Line[i] != ';' && Line[i] != '\n') {
					newFileName[i] = Line[i];
					++i;
				}
				newFileName[i] = '\0';
				snprintf(Buffer, 4096, "%s/%s.prf", DirType, newFileName);
				newFile = fopen(Buffer, "w");
				if (newFile == NULL) {
					fprintf(stderr, "Line %lu, Cannot create file %s\n>%s", iline, Buffer, LineAC);
					exit(1);
				}
				fputs(LineID, newFile);
				fputs(LineAC, newFile);
				OpenProfile = true;
						}
						else {
				if (!OpenProfile ) {
					if (Line[0] != 'C' && Line[1] != 'C') {
						fprintf(stderr, "Line %lu, should not have no profile opened\n", iline);
						exit(1);
					}
				} else 
					fputs(Buffer, newFile);
      }  
  }
  
  fclose(in);
  return 0;
}
