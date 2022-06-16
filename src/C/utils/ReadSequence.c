/* TODO:
 *  Use the information from the filesystem to adjust buffer size reading.
 */
#define _GNU_SOURCE
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <string.h>
#include "../include/pfSequence.h"

/* Beware of alignment before changing these !!! */
#define BUFFER_SIZE         (1024*1024)
#define CHAIN_SEGMENT_SIZE      16384

/* NOTE: For mmap reading we extend the sequence count by one and store
 *       the file size in the extra data, that is the sequence count + 1
 * 	 structure should have
 * 		off_t Offset = file size in byte;
 *		size_t HeaderLength = 0;
 *		size_t SequenceLength = 0;
 */

int AnalyzeFASTAStructure(const char * const FileName, DBSequence_t * const Info)
{
	typedef struct chain_s {
		off_t Offset;
		size_t HeaderLength;
		off_t QualityOffset;
		size_t QualityHeaderLength;
	} chain_t;
  char Buffer[BUFFER_SIZE] __attribute__((aligned(16)));
  chain_t * * chain_ptr ;
  s_Data * restrict DataPtr = NULL;
	size_t MaxSequenceSize = 0Ul;
	size_t chain_size = 819UL;
	struct stat FileStat;

  strncpy(Info->FileName, FileName, 256);

  const int FileDescriptor = open(FileName, O_RDONLY);
  if ( FileDescriptor == -1) {
    fprintf(stderr, "Error occured while accessing file %s\n", FileName);
    perror("The error reads ");
    return 1;
  }

  if (fstat(FileDescriptor, &FileStat) == -1 ) {
    fprintf(stderr, "Error occured while accessing file %s\n", FileName);
    perror("The error reads ");
    return 1;
  }

  if (FileStat.st_size == 0 )
    return 1;
  else
    Info->FileSize = FileStat.st_size;

  /* Allocate initial CHAIN */
	{
		register chain_t * const data = malloc(CHAIN_SEGMENT_SIZE*sizeof(chain_t));
		chain_ptr = (chain_t**) malloc(chain_size*sizeof(chain_t*));
		if (data == NULL || chain_ptr == NULL) {
			fputs("Unable to allocate sufficient memory.\n", stderr);
			return 1;
		}
		chain_ptr[0] = data;
	}

  size_t ChainCount = 1UL;
  ssize_t SequenceCount  = -1L;
	chain_t * CurrentChain = chain_ptr[0] - 1 ;

  off_t BufferOffset  = lseek (FileDescriptor, 0, SEEK_CUR);
  size_t length       = read(FileDescriptor, Buffer, BUFFER_SIZE*sizeof(char));
	_Bool HeaderPending = false;
  if (length > 0UL) {
		_Bool First = true;
		do {
			for (size_t i=0UL; i<length; ++i) {
				const size_t CurrentOffset = BufferOffset + (off_t) i;
				if (Buffer[i] == '>') {
					if (!First)
						CurrentChain->QualityOffset = CurrentOffset;
					else
						First = false;
					if ( ++SequenceCount == CHAIN_SEGMENT_SIZE) {
						if (++ChainCount == chain_size) {
							chain_size += 64;
							chain_ptr = realloc(chain_ptr, chain_size*sizeof(chain_t));
							if (chain_ptr == NULL) {
								fputs("Unable to reallocate memory\n", stderr);
								goto FreeMemory;
							}
						}
						{
							register void * const data = malloc(CHAIN_SEGMENT_SIZE*sizeof(chain_t));
							if (data == NULL) {
								fputs("Unable to allocate sufficient memory.\n", stderr);
								goto FreeMemory;
							}
							chain_ptr[ChainCount-1] = data;
						}
						CurrentChain = chain_ptr[ChainCount-1];
						SequenceCount = 0UL;
					}
					else {
						CurrentChain++;
					}
					CurrentChain->Offset = CurrentOffset;
					CurrentChain->QualityOffset = 0UL;
					HeaderPending  = true;
				}
				else if (Buffer[i] == '\n') {
					if (HeaderPending) {
						CurrentChain->HeaderLength = (size_t) ( CurrentOffset - CurrentChain->Offset );
						HeaderPending = false;
					}
				}
			}

			BufferOffset = lseek (FileDescriptor, 0, SEEK_CUR);
			length = read(FileDescriptor, Buffer, BUFFER_SIZE*sizeof(char));
		} while (length > 0);

		if (/*CurrentChain >= chain_ptr[0] && */CurrentChain->QualityOffset == 0UL)
			CurrentChain->QualityOffset = FileStat.st_size;


    /* Package the overall data */
    SequenceCount += 1UL + (ChainCount-1)*CHAIN_SEGMENT_SIZE;

    DataPtr = (s_Data*) malloc((1+SequenceCount)*sizeof(s_Data));
    if ( DataPtr == NULL) {
      fputs("Unable to allocate sufficient memory\n", stderr);
      goto FreeMemory;
    }


    /* Fill up the data */
    CurrentChain = chain_ptr[0];
    ChainCount = 0;
    for (size_t i=0; i<SequenceCount-1; ++i) {
      const size_t index = i % CHAIN_SEGMENT_SIZE;
      if ( index != (CHAIN_SEGMENT_SIZE-1) ) {
        DataPtr[i].Sequence.Offset                  = CurrentChain[index].Offset;
        DataPtr[i].Sequence.HeaderLength            = (unsigned int) CurrentChain[index].HeaderLength;
        register const size_t SequenceLength        =  CurrentChain[index].QualityOffset -  CurrentChain[index].Offset;
        DataPtr[i].Sequence.SequenceLength          = (unsigned int) (SequenceLength - CurrentChain[index].HeaderLength - 1);

				DataPtr[i].Quality.Offset                   = CurrentChain[index].QualityOffset;
				DataPtr[i].Quality.HeaderLength             = (unsigned int) CurrentChain[index].QualityHeaderLength;
				register const size_t QualitySequenceLength = CurrentChain[index+1].Offset -  CurrentChain[index].QualityOffset;
				DataPtr[i].Quality.SequenceLength           = (unsigned int) (QualitySequenceLength - CurrentChain[index].QualityHeaderLength);
				if (DataPtr[i].Quality.SequenceLength	> 0) DataPtr[i].Quality.SequenceLength -= 1;

				const size_t FullItemLength = CurrentChain[index].QualityOffset - CurrentChain[index].Offset;
				MaxSequenceSize             = FullItemLength > MaxSequenceSize? FullItemLength : MaxSequenceSize;

      } else {
        DataPtr[i].Sequence.Offset           = CurrentChain[index].Offset;
        DataPtr[i].Sequence.HeaderLength     = (unsigned int) CurrentChain[index].HeaderLength;
				register const size_t SequenceLength = CurrentChain[index].QualityOffset -  CurrentChain[index].Offset;
				DataPtr[i].Sequence.SequenceLength   = (unsigned int) (SequenceLength - DataPtr[i].Sequence.HeaderLength - 1);

				DataPtr[i].Quality.Offset       = CurrentChain[index].QualityOffset;
				DataPtr[i].Quality.HeaderLength = (unsigned int) CurrentChain[index].QualityHeaderLength;

				const size_t FullItemLength = CurrentChain[index].QualityOffset - CurrentChain[index].Offset;
				MaxSequenceSize             = FullItemLength > MaxSequenceSize? FullItemLength : MaxSequenceSize;

        free(chain_ptr[ChainCount++]);
        CurrentChain = chain_ptr[ChainCount];

				register const size_t QualitySequenceLength = CurrentChain[0].Offset -  DataPtr[i].Quality.Offset;
				DataPtr[i].Quality.SequenceLength           = (unsigned int) (QualitySequenceLength - DataPtr[i].Quality.HeaderLength);
				if (DataPtr[i].Quality.SequenceLength	> 0) DataPtr[i].Quality.SequenceLength -= 1;
      }
    }

    const size_t index = (SequenceCount-1) % CHAIN_SEGMENT_SIZE;

    DataPtr[SequenceCount-1].Sequence.Offset         = CurrentChain[index].Offset;
    DataPtr[SequenceCount-1].Sequence.HeaderLength   = (unsigned int) CurrentChain[index].HeaderLength;
    register const size_t SequenceLength             = CurrentChain[index].QualityOffset - CurrentChain[index].Offset;
		DataPtr[SequenceCount-1].Sequence.SequenceLength = (unsigned int) (SequenceLength - CurrentChain[index].HeaderLength - 1);
		register const size_t QualitySequenceLength      = FileStat.st_size -  CurrentChain[index].QualityOffset;
		DataPtr[index].Quality.SequenceLength = (unsigned int) (QualitySequenceLength - CurrentChain[index].QualityHeaderLength);
		if ( DataPtr[index].Quality.SequenceLength	> 0) DataPtr[index].Quality.SequenceLength -= 1;
		DataPtr[index].Quality.HeaderLength = CurrentChain[index].QualityHeaderLength;
		DataPtr[index].Quality.Offset = CurrentChain[index].QualityOffset;

		register const size_t FullItemLength = CurrentChain[index].QualityOffset - CurrentChain[index].Offset;
		MaxSequenceSize = FullItemLength > MaxSequenceSize? FullItemLength : MaxSequenceSize;

    free(chain_ptr[ChainCount]);
  }
  close(FileDescriptor);

  /* Add extra file size at last DataPtr */
  DataPtr[SequenceCount].Sequence.Offset = FileStat.st_size;
  DataPtr[SequenceCount].Sequence.HeaderLength = 0;
  DataPtr[SequenceCount].Sequence.SequenceLength = 0;
	DataPtr[SequenceCount].Quality.Offset = FileStat.st_size;
	DataPtr[SequenceCount].Quality.HeaderLength = 0;
	DataPtr[SequenceCount].Quality.SequenceLength = 0;

  Info->DataPtr          = DataPtr;
  Info->SequenceCount    = SequenceCount;
  Info->MaxSequenceSize  = MaxSequenceSize+1;
	Info->LastModification = FileStat.st_mtim;

  return 0;

  FreeMemory:
    close(FileDescriptor);
    for (size_t i=0; i<ChainCount; ++i) {
      free(chain_ptr[i]);
    }
  return 1;
}

int AnalyzeFASTQStructure(const char * const FileName, DBSequence_t * const Info)
{
	typedef struct chain_s {
		off_t Offset;
		size_t HeaderLength;
		off_t QualityOffset;
		size_t QualityHeaderLength;
	} chain_t;
  char Buffer[BUFFER_SIZE] __attribute__((aligned(16)));
  chain_t * * chain_ptr ;
  s_Data * restrict DataPtr = NULL;
	size_t MaxSequenceSize = 0Ul;
	size_t chain_size = 819UL;
	struct stat FileStat;

  strncpy(Info->FileName, FileName, 256);

  const int FileDescriptor = open(FileName, O_RDONLY);
  if ( FileDescriptor == -1) {
    fprintf(stderr, "Error occured while accessing file %s\n", FileName);
    perror("The error reads ");
    return 1;
  }

  if (fstat(FileDescriptor, &FileStat) == -1 ) {
    fprintf(stderr, "Error occured while accessing file %s\n", FileName);
    perror("The error reads ");
    return 1;
  }

  if (FileStat.st_size == 0 )
    return 1;
  else
    Info->FileSize = FileStat.st_size;

  /* Allocate initial CHAIN */
	{
		register chain_t * const data = malloc(CHAIN_SEGMENT_SIZE*sizeof(chain_t));
		chain_ptr = (chain_t**) malloc(chain_size*sizeof(chain_t*));
		if (data == NULL || chain_ptr == NULL) {
			fputs("Unable to allocate sufficient memory.\n", stderr);
			return 1;
		}
		chain_ptr[0] = data;
	}

  size_t ChainCount = 1UL;
  ssize_t SequenceCount  = -1L;
	chain_t * CurrentChain = chain_ptr[0] - 1 ;

  off_t BufferOffset  = lseek (FileDescriptor, 0, SEEK_CUR);
  size_t length       = read(FileDescriptor, Buffer, BUFFER_SIZE*sizeof(char));
	_Bool HeaderPending = false;
  if (length > 0UL) {
		_Bool QualityScorePending = false;
		_Bool QualityHeaderPending = false;
		do {
			for (size_t i=0UL; i<length; ++i) {
// 					fputc(Buffer[i], stdout);
				const size_t CurrentOffset = BufferOffset + (off_t) i;
				if (Buffer[i] == '\n') {
					if (HeaderPending) {
						CurrentChain->HeaderLength = (size_t) ( CurrentOffset - CurrentChain->Offset );
						HeaderPending = false;
// 							printf("HEADER STOP\n");
					}
					else if (QualityHeaderPending) {
						CurrentChain->QualityHeaderLength = (size_t) ( CurrentOffset - CurrentChain->QualityOffset );
						QualityHeaderPending = false;
// 							printf("QUALITY HEADER STOP\nQUALITY_SCORE START\n");
						QualityScorePending = true;
					}
					else if (QualityScorePending) {
						QualityScorePending = false;
// 							printf("QUALITY STOP\n");
					}
				}
				else if ( Buffer[i] == '@' && !QualityScorePending && !HeaderPending && !QualityHeaderPending) {
					if ( ++SequenceCount == CHAIN_SEGMENT_SIZE) {
						if (++ChainCount == chain_size) {
							chain_size += 64;
							chain_ptr = realloc(chain_ptr, chain_size*sizeof(chain_t));
							if (chain_ptr == NULL) {
								fputs("Unable to reallocate memory\n", stderr);
								goto FreeMemory;
							}
						}
						{
							register void * const data = malloc(CHAIN_SEGMENT_SIZE*sizeof(chain_t));
							if (data == NULL) {
								fputs("Unable to allocate sufficient memory.\n", stderr);
								goto FreeMemory;
							}
							chain_ptr[ChainCount-1] = data;
						}
						CurrentChain = chain_ptr[ChainCount-1];
						SequenceCount = 0UL;
					}
					else {
						++CurrentChain;
					}
					CurrentChain->Offset = CurrentOffset;
					HeaderPending  = true;
// 						printf("HEADER START\n");
				}
				else if (Buffer[i] == '+' && !QualityHeaderPending && !HeaderPending && !QualityScorePending) {
					CurrentChain->QualityOffset = CurrentOffset;
					QualityHeaderPending = true;
// 						printf("QUALITY HEADER START\n");
				}
			}

			BufferOffset  = lseek (FileDescriptor, 0, SEEK_CUR);
			length = read(FileDescriptor, Buffer, BUFFER_SIZE*sizeof(char));
		} while (length > 0);


    /* Package the overall data */
    SequenceCount += 1UL + (ChainCount-1)*CHAIN_SEGMENT_SIZE;

    DataPtr = (s_Data*) malloc((1+SequenceCount)*sizeof(s_Data));
    if ( DataPtr == NULL) {
      fputs("Unable to allocate sufficient memory\n", stderr);
      goto FreeMemory;
    }


    /* Fill up the data */
    CurrentChain = chain_ptr[0];
    ChainCount = 0;
    for (size_t i=0; i<SequenceCount-1; ++i) {
      const size_t index = i % CHAIN_SEGMENT_SIZE;
      if ( index != (CHAIN_SEGMENT_SIZE-1) ) {
        DataPtr[i].Sequence.Offset                  = CurrentChain[index].Offset;
        DataPtr[i].Sequence.HeaderLength            = (unsigned int) CurrentChain[index].HeaderLength;
        register const size_t SequenceLength        =  CurrentChain[index].QualityOffset -  CurrentChain[index].Offset;
        DataPtr[i].Sequence.SequenceLength          = (unsigned int) (SequenceLength - CurrentChain[index].HeaderLength - 1);

				DataPtr[i].Quality.Offset                   = CurrentChain[index].QualityOffset;
				DataPtr[i].Quality.HeaderLength             = (unsigned int) CurrentChain[index].QualityHeaderLength;
				register const size_t QualitySequenceLength = CurrentChain[index+1].Offset -  CurrentChain[index].QualityOffset;
				DataPtr[i].Quality.SequenceLength           = (unsigned int) (QualitySequenceLength - CurrentChain[index].QualityHeaderLength);
				if (DataPtr[i].Quality.SequenceLength	> 0) DataPtr[i].Quality.SequenceLength -= 2;

				const size_t FullItemLength = CurrentChain[index].QualityOffset - CurrentChain[index].Offset;
				MaxSequenceSize             = FullItemLength > MaxSequenceSize? FullItemLength : MaxSequenceSize;

      } else {
        DataPtr[i].Sequence.Offset           = CurrentChain[index].Offset;
        DataPtr[i].Sequence.HeaderLength     = (unsigned int) CurrentChain[index].HeaderLength;
				register const size_t SequenceLength = CurrentChain[index].QualityOffset -  CurrentChain[index].Offset;
				DataPtr[i].Sequence.SequenceLength   = (unsigned int) (SequenceLength - DataPtr[i].Sequence.HeaderLength - 1);

				DataPtr[i].Quality.Offset       = CurrentChain[index].QualityOffset;
				DataPtr[i].Quality.HeaderLength = (unsigned int) CurrentChain[index].QualityHeaderLength;

				const size_t FullItemLength = CurrentChain[index].QualityOffset - CurrentChain[index].Offset;
				MaxSequenceSize             = FullItemLength > MaxSequenceSize? FullItemLength : MaxSequenceSize;

        free(chain_ptr[ChainCount++]);
        CurrentChain = chain_ptr[ChainCount];

				register const size_t QualitySequenceLength = CurrentChain[0].Offset -  DataPtr[i].Quality.Offset;
				DataPtr[i].Quality.SequenceLength           = (unsigned int) (QualitySequenceLength - DataPtr[i].Quality.HeaderLength);
				if (DataPtr[i].Quality.SequenceLength	> 0) DataPtr[i].Quality.SequenceLength -= 2;
      }
    }

    const size_t index = (SequenceCount-1) % CHAIN_SEGMENT_SIZE;

    DataPtr[SequenceCount-1].Sequence.Offset         = CurrentChain[index].Offset;
    DataPtr[SequenceCount-1].Sequence.HeaderLength   = (unsigned int) CurrentChain[index].HeaderLength;
    register const size_t SequenceLength             = CurrentChain[index].QualityOffset - CurrentChain[index].Offset;
		DataPtr[SequenceCount-1].Sequence.SequenceLength = (unsigned int) (SequenceLength - CurrentChain[index].HeaderLength - 1);
		register const size_t QualitySequenceLength      = FileStat.st_size -  CurrentChain[index].QualityOffset;
		DataPtr[index].Quality.SequenceLength = (unsigned int) (QualitySequenceLength - CurrentChain[index].QualityHeaderLength);
		if ( DataPtr[index].Quality.SequenceLength	> 0) DataPtr[index].Quality.SequenceLength -= 1;
		DataPtr[index].Quality.HeaderLength = CurrentChain[index].QualityHeaderLength;
		DataPtr[index].Quality.Offset = CurrentChain[index].QualityOffset;

		register const size_t FullItemLength = CurrentChain[index].QualityOffset - CurrentChain[index].Offset;
		MaxSequenceSize = FullItemLength > MaxSequenceSize? FullItemLength : MaxSequenceSize;

    free(chain_ptr[ChainCount]);
  }
  close(FileDescriptor);

  /* Add extra file size at last DataPtr */
  DataPtr[SequenceCount].Sequence.Offset = FileStat.st_size;
  DataPtr[SequenceCount].Sequence.HeaderLength = 0;
  DataPtr[SequenceCount].Sequence.SequenceLength = 0;
	DataPtr[SequenceCount].Quality.Offset = FileStat.st_size;
	DataPtr[SequenceCount].Quality.HeaderLength = 0;
	DataPtr[SequenceCount].Quality.SequenceLength = 0;

  Info->DataPtr          = DataPtr;
  Info->SequenceCount    = SequenceCount;
  Info->MaxSequenceSize  = MaxSequenceSize+1;
	Info->LastModification = FileStat.st_mtim;

  return 0;

  FreeMemory:
    close(FileDescriptor);
    for (size_t i=0; i<ChainCount; ++i) {
      free(chain_ptr[i]);
    }
  return 1;
}

int AnalyzeEMBLStructure(const char * const FileName, DBSequence_t * const Info)
{
	typedef struct chain_s {
		off_t Offset;
		size_t HeaderLength;
		size_t SequenceEnd;
	} chain_t;

  chain_t * * chain_ptr ;
  s_Data * restrict DataPtr = NULL;
	size_t MaxSequenceSize = 0Ul;
	size_t chain_size = 819UL;
	struct stat FileStat;

  strncpy(Info->FileName, FileName, 256);

	FILE* const FileDescriptor = fopen(FileName, "r" );
  if ( FileDescriptor == NULL) {
    fprintf(stderr, "Error occured while accessing file %s\n", FileName);
    perror("The error reads ");
    return 1;
  }

  if (stat(FileName, &FileStat) == -1 ) {
    fprintf(stderr, "Error occured while accessing file %s\n", FileName);
    perror("The error reads ");
    return 1;
  }

  if (FileStat.st_size == 0 )
    return 1;
  else
    Info->FileSize = FileStat.st_size;

  /* Allocate initial CHAIN */
	{
		register chain_t * const data = malloc(CHAIN_SEGMENT_SIZE*sizeof(chain_t));
		chain_ptr = (chain_t**) malloc(chain_size*sizeof(chain_t*));
		if (data == NULL || chain_ptr == NULL) {
			fputs("Unable to allocate sufficient memory.\n", stderr);
			return 1;
		}
		chain_ptr[0] = data;
	}

  size_t ChainCount = 1UL;
  ssize_t SequenceCount  = -1L;
	chain_t * CurrentChain = chain_ptr[0] - 1 ;

  off_t BufferOffset = 0UL;
  char * Line = NULL;
	size_t LineSize = 0UL;
	while(!feof(FileDescriptor)) {
			ssize_t length = getline(&Line, &LineSize, FileDescriptor);
			if (length >= 2) {
				if (Line[0] == 'S' && Line[1] == 'Q') {
					if ( ++SequenceCount == CHAIN_SEGMENT_SIZE) {
						if (++ChainCount == chain_size) {
							chain_size += 64;
							chain_ptr = realloc(chain_ptr, chain_size*sizeof(chain_t));
							if (chain_ptr == NULL) {
								fputs("Unable to reallocate memory\n", stderr);
								goto FreeMemory;
							}
						}
						{
							register void * const data = malloc(CHAIN_SEGMENT_SIZE*sizeof(chain_t));
							if (data == NULL) {
								fputs("Unable to allocate sufficient memory.\n", stderr);
								goto FreeMemory;
							}
							chain_ptr[ChainCount-1] = data;
						}
						CurrentChain = chain_ptr[ChainCount-1];
						SequenceCount = 0UL;
					}
					else {
						++CurrentChain;
					}
					CurrentChain->Offset = BufferOffset + 2;
					CurrentChain->HeaderLength = length - 2;
				}
				else if (Line[0] == '/' && Line[1] == '/') {
					CurrentChain->SequenceEnd = BufferOffset;
				}
			}
			BufferOffset += (off_t) length;
	}


	/* Package the overall data */
	SequenceCount += 1UL + (ChainCount-1)*CHAIN_SEGMENT_SIZE;

	DataPtr = (s_Data*) malloc((1+SequenceCount)*sizeof(s_Data));
	if ( DataPtr == NULL) {
		fputs("Unable to allocate sufficient memory\n", stderr);
		goto FreeMemory;
	}


	/* Fill up the data */
	CurrentChain = chain_ptr[0];
	ChainCount = 0;
	for (size_t i=0; i<SequenceCount; ++i) {
		const size_t index = i % CHAIN_SEGMENT_SIZE;
		DataPtr[i].Sequence.Offset         = CurrentChain[index].Offset;
		DataPtr[i].Sequence.HeaderLength   = (unsigned int) CurrentChain[index].HeaderLength - 1;
		const size_t FullItemLength        = (unsigned int) (CurrentChain[index].SequenceEnd - CurrentChain[index].Offset) + 1;
		DataPtr[i].Sequence.SequenceLength = (unsigned int) (FullItemLength - CurrentChain[index].HeaderLength);

		DataPtr[i].Quality.Offset          = (off_t) 0;
		DataPtr[i].Quality.HeaderLength    = 0U;
		DataPtr[i].Quality.SequenceLength  = 0U;

		MaxSequenceSize             = FullItemLength > MaxSequenceSize? FullItemLength : MaxSequenceSize;

		if (index == (CHAIN_SEGMENT_SIZE-1)) {
			free(chain_ptr[ChainCount++]);
			CurrentChain = chain_ptr[ChainCount];
		}
	}

	free(chain_ptr[ChainCount]);
  fclose(FileDescriptor);

  /* Add extra file size at last DataPtr */
  DataPtr[SequenceCount].Sequence.Offset = FileStat.st_size;
  DataPtr[SequenceCount].Sequence.HeaderLength = 0;
  DataPtr[SequenceCount].Sequence.SequenceLength = 0;
	DataPtr[SequenceCount].Quality.Offset = FileStat.st_size;
	DataPtr[SequenceCount].Quality.HeaderLength = 0;
	DataPtr[SequenceCount].Quality.SequenceLength = 0;

  Info->DataPtr          = DataPtr;
  Info->SequenceCount    = SequenceCount;
  Info->MaxSequenceSize  = MaxSequenceSize+1;
	Info->LastModification = FileStat.st_mtim;
  free(Line);

  return 0;

  FreeMemory:;
    if (Line) free(Line);
    fclose(FileDescriptor);
    for (size_t i=0; i<ChainCount; ++i) {
      free(chain_ptr[i]);
    }
  return 1;
}

int ExportDBStructure(FILE* const stream, const DBSequence_t * const Info)
{

  unsigned short FileNameLength = (unsigned short) strlen(Info->FileName);
  if (fwrite(&(FileNameLength), sizeof(unsigned short), 1, stream) != 1) return 1;
  if (fwrite(Info->FileName, sizeof(char), (size_t) FileNameLength, stream) != (size_t) FileNameLength) return 1;
  if (fwrite(&(Info->FileSize), sizeof(off_t), 1, stream) != 1) return 1;
  unsigned long tmp = (unsigned long) Info->SequenceCount;
  if (fwrite(&tmp, sizeof(unsigned long), 1, stream) != 1) return 1;
  tmp = (unsigned long) Info->MaxSequenceSize;
  if (fwrite(&tmp, sizeof(unsigned long), 1, stream) != 1) return 1;
	if (fwrite(&(Info->LastModification), sizeof(struct timespec), 1, stream) != 1) return 1;
  if (fwrite(Info->DataPtr, sizeof(s_Data), 1+Info->SequenceCount, stream) != 1+Info->SequenceCount) return 1;

 return 0;
}

int ImportDBStructure(FILE* const stream, const char * const DBFileName, DBSequence_t * const Info)
{
  unsigned short FileNameLength;
  unsigned long SequenceCount, MaxSequenceSize;
  memset(Info->FileName, 0, 256*sizeof(char));


  if (fread(&FileNameLength, sizeof(unsigned short), 1, stream) != 1) return 1;
  if (fread(Info->FileName, sizeof(char), (size_t) FileNameLength, stream) != FileNameLength) return 1;
  if (fread(&(Info->FileSize), sizeof(off_t), 1, stream) != 1) return 1;
  if (fread(&SequenceCount, sizeof(unsigned long), 1, stream) != 1) return 1;
  Info->SequenceCount = (size_t) SequenceCount;
  if (fread(&MaxSequenceSize, sizeof(unsigned long), 1, stream) != 1) return 1;
	if (fread(&(Info->LastModification), sizeof(struct timespec), 1, stream) != 1) return 1;
  Info->MaxSequenceSize = (size_t) MaxSequenceSize;
  Info->DataPtr = (s_Data*) malloc((1+Info->SequenceCount)*sizeof(s_Data));
  if ( Info->DataPtr == NULL) {
    fputs("Unable to allocate sufficient memory to hold FASTA structure\n", stderr);
    return 1;
  }
  if (fread(Info->DataPtr, sizeof(s_Data), (1+Info->SequenceCount), stream) != 1+Info->SequenceCount) return 1;

	/* Check size and time of last modification */
	struct stat st;
	if (stat(DBFileName, &st) != 0) {
		perror("ImportDBStructure");
		return 1;
	}

	if (st.st_size != Info->FileSize || st.st_mtim.tv_nsec != Info->LastModification.tv_nsec ||
		st.st_mtim.tv_sec != Info->LastModification.tv_sec) {
		fprintf(stderr, "Index file does not match %s\n", DBFileName);
		return 1;
	}

 return 0;
}

#ifdef _TEST
int main(int argc, char *argv[])
{
  char Buffer[2048] __attribute__((aligned(16)));
  DBSequence_t FASTA;

  if (argc < 2 || argc > 3) { fputs("Provide just a FASTA file\n", stderr); return 1; }
  printf("Starting analysis of %s...\n", argv[1]);

  const int res = AnalyzeFASTAStructure(argv[1], &FASTA);
  if (res != 0) {
    fputs("Error found.\n", stderr);
  } else {
    printf("FASTA file %s analyzed\n\tFound %lu sequences within %lu bytes\n\tBiggest sequence entry is %lu bytes\n",
           argv[1],
           FASTA.SequenceCount,
           FASTA.FileSize,
           FASTA.MaxSequenceSize);
    if (argc > 2) {
      const int index = (size_t) atoi(argv[2]);
      if ( index < 0) {
        sprintf(Buffer, "%s.copy", argv[1]);
        printf("Recreating global file for diff command in %s\n", Buffer);
        FILE* const out = fopen(Buffer, "w");
        FILE* const in  = fopen(argv[1], "r");
        for (size_t i=0; i<FASTA.SequenceCount; ++i) {
          fseek(in, (long) FASTA.DataPtr[i].Offset, SEEK_SET);
          fread(Buffer, FASTA.DataPtr[i].HeaderLength, sizeof(char), in);
          Buffer[FASTA.DataPtr[i].HeaderLength] = '\0';
          fprintf(out,"%s\n", Buffer);
          if (FASTA.DataPtr[i].SequenceLength >= 2048) {
            fprintf(stderr, "Sorry but allocated buffer cannot hold sequence length of %lu\n",FASTA.DataPtr[i].SequenceLength );
            return 1;
          }
          fread(Buffer, FASTA.DataPtr[i].SequenceLength, sizeof(char), in);
          Buffer[FASTA.DataPtr[i].SequenceLength] = '\0';
          fprintf(out, "%s\n", &Buffer[1]);
        }
        fclose(in);
        fclose(out);
      } else if (index < FASTA.SequenceCount) {
        printf("Sequence %i is at offset %li with header size of %lu and sequence size of %lu\n",
               index,
               FASTA.DataPtr[index].Offset,
               FASTA.DataPtr[index].HeaderLength,
               FASTA.DataPtr[index].SequenceLength);

        FILE* const in = fopen(argv[1], "r");
        fseek(in, (long) FASTA.DataPtr[index].Offset, SEEK_SET);
        fread(Buffer, FASTA.DataPtr[index].HeaderLength, sizeof(char), in);
        Buffer[FASTA.DataPtr[index].HeaderLength] = '\0';
        printf("Header : %s\n", Buffer);
        if (FASTA.DataPtr[index].SequenceLength < 2048) {
          fread(Buffer, FASTA.DataPtr[index].SequenceLength, sizeof(char), in);
          Buffer[FASTA.DataPtr[index].SequenceLength] = '\0';
          printf("Sequence : %s\n", Buffer);
        } else {
          fputs("Sorry but allocated buffer cannot hold sequence length\n", stderr);
        }
        fclose(in);
      }
    }

  }

  return 0;
}
#endif
