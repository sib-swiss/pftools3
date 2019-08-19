/***************************************************************************************************
                        PFTOOLS
 ***************************************************************************************************
  Oct 3, 2011 pfSequenceInline.h
 ***************************************************************************************************
 (C) 2011 SIB Swiss Institute of Bioinformatics
     Thierry Schuepbach (thierry.schuepbach@sib.swiss)
 ***************************************************************************************************/
extern inline char * __ALWAYS_INLINE CleanSequence(PFSequence * const Sequence)
{
  register size_t counter = 0;
  unsigned char * restrict const CharPtr = Sequence->ProfileIndex;

  for (size_t i=0; i<Sequence->Length; ++i) {
    register const unsigned char c = ((unsigned char) CharPtr[i] >= (unsigned char) 'a' ) 
                                    ? (unsigned char) CharPtr[i] - ((unsigned char) 'a' - (unsigned char) 'A')
				    : (unsigned char) CharPtr[i];
    if ( c >= (unsigned char) 'A' && c <= (unsigned char) 'Z' ) {
      CharPtr[counter++] = c;
    } 
  }
  
  Sequence->Length = counter;
  return Sequence->ProfileIndex;
}

static inline PFSequence * ReadSequenceIndex(Sequence * const Seq, FILE * const stream, const ss_Data * const DataPtr)
{
  /* Position into file */
  fseek(stream, (long) DataPtr->Offset, SEEK_SET);

  /* Read into memory */
  register const size_t Size = DataPtr->HeaderLength + DataPtr->SequenceLength + 1; // + 1 accounts for \n in between
  if (fread(Seq->Data.Header, sizeof(char), Size, stream) != Size) return NULL;

  /* Bound text */
  Seq->Data.Header[DataPtr->HeaderLength] = '\0';
  Seq->Data.Header[Size] = '\0';

  /* Set Sequence data start */
  Seq->ProfileData.ProfileIndex = (unsigned char*) &(Seq->Data.Header[DataPtr->HeaderLength+1]);

  /* Set Sequence length */
  Seq->ProfileData.Length = DataPtr->SequenceLength;

  return &(Seq->ProfileData);
}

#ifdef __USE_MMAP__
static inline PFSequence * MMAP_ReadSequenceIndex(Sequence * const Seq, const char * const restrict Array,
                                                  const ss_Data * const DataPtr, const off_t InitialArrayOffset
#ifdef MMAP_DEBUG
                                                  ,const size_t ThreadId, const size_t NodeId, const size_t length
#endif						
)
{
  /* Position into Array */
  register const size_t ArrayOffset = (size_t) ( DataPtr->Offset - InitialArrayOffset);
//   fseek(stream, (long) DataPtr->Offset, SEEK_SET);

  /* Read into memory */
  register const size_t Size = DataPtr->HeaderLength + DataPtr->SequenceLength + 1; // + 1 accounts for \n in between
#ifdef MMAP_DEBUG
  if ( Size + ArrayOffset > length) {
    fprintf(stderr,"Thread %lu from Node %lu will read beyond mmap %lu > %lu\n", ThreadId, NodeId, Size + ArrayOffset, length);
  }
#endif
  memcpy(Seq->Data.Header, &Array[ArrayOffset], sizeof(char)*Size);
  
//   if (fread(Seq->Data.Header, sizeof(char), Size, stream) != Size) return NULL;

  /* Bound text */
  Seq->Data.Header[DataPtr->HeaderLength] = '\0';
  Seq->Data.Header[Size] = '\0';

  /* Set Sequence data start */
  Seq->ProfileData.ProfileIndex = (unsigned char*) &(Seq->Data.Header[DataPtr->HeaderLength+1]);

  /* Set Sequence length */
  Seq->ProfileData.Length = DataPtr->SequenceLength;

  return &(Seq->ProfileData);
}
#endif

static inline void ReadSequenceNameIndex(char * const Name, FILE * const stream, const ss_Data * const DataPtr)
{
  /* Position into file */
  fseek(stream, (long int) DataPtr->Offset, SEEK_SET);

  /* Read into memory */
  if (fscanf(stream, ">%s",Name) != 1) {
    fprintf(stderr, "Read error for sequence name @ offset %lu\n", DataPtr->Offset);
  }
}

extern inline unsigned char __ALWAYS_INLINE TranslateCharToIndex(const char letter, const unsigned char * restrict const Alphabet)
{
  const unsigned char lletter = (unsigned char) letter;
  register size_t index = (size_t) ( ( lletter >= (unsigned char) 'a' ) ? lletter - ((unsigned char) 'a' - (unsigned char) 'A') : lletter );
  if ( index >= (size_t) 'A' && index <= (size_t) 'Z' ) {
    return Alphabet[index - (size_t) 'A'];
  } else {
    return 0;
  }
}

/* WARNING: NEED TO OPTIMIZE THIS PIECE OF JUNK */
extern inline PFSequence * __ALWAYS_INLINE TranslateSequenceToIndex(PFSequence * const Sequence, const unsigned char * restrict const Alphabet )
{
  register size_t counter = 0;
  unsigned char * restrict const CharPtr = Sequence->ProfileIndex;
  Sequence->OriginalSequence = malloc(Sequence->Length * sizeof(unsigned char)); // FIXME - not sure where to free this

  for (size_t i=0; i<Sequence->Length; ++i) {
    Sequence->OriginalSequence[counter] = ( CharPtr[i] >= (unsigned char) 'a' ) ? CharPtr[i] - ((unsigned char) 'a' - (unsigned char) 'A') : CharPtr[i];
    register size_t index = (size_t) ( ( CharPtr[i] >= (unsigned char) 'a' ) ? CharPtr[i] - ((unsigned char) 'a' - (unsigned char) 'A') : CharPtr[i] );
    if ( index >= (size_t) 'A' && index <= (size_t) 'Z' ) {
#ifdef XALIT_DEBUG
      fprintf(stderr,"CharPtr[%zu++] = Alphabet[%zu (%c) - (size_t) 'A'] = %u\n", counter, index, CharPtr[i], Alphabet[index - (size_t) 'A']);
#endif
      CharPtr[counter++] = Alphabet[index - (size_t) 'A'];
    }
  }

  Sequence->Length = counter;
  return Sequence;
}

extern inline void __ALWAYS_INLINE ReverseTranslatedSequence(PFSequence * const Sequence)
{
  unsigned char * restrict const CharPtr = Sequence->ProfileIndex;
  const size_t SeqLength = Sequence->Length;
  unsigned char * BackPtr = &CharPtr[SeqLength-1];

  for (size_t i=0; i<SeqLength/2; ++i) {
      const unsigned char c = CharPtr[i];
      CharPtr[i] = *BackPtr;
      *BackPtr-- = c;
  }
  if (SeqLength & 0x1) {
      const unsigned char c = CharPtr[SeqLength/2];
      CharPtr[Sequence->Length/2] = *BackPtr;
      *BackPtr = c;
  }
}
