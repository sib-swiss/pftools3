/*******************************************************
                        PFTOOLS
 *******************************************************
  Jan 18, 2011 output.c
 *******************************************************
 (C) 2011 SIB Swiss Institute of Bioinformatics
     Thierry Schuepbach (thierry.schuepbach@sib.swiss)
 *******************************************************/
#include "config.h"
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include "../include/pfProfile.h"
#include "../include/pfSequence.h"

/*
 * WARNING: alignement does not start at 0 but 1 !!!
 *          NEEDS TO BE FIXED SOME DAY
 */
static const char Dflt[] = "unknown"; 
	
void PrintSimple(const struct Profile * const prf, const char * * const AlignedSequence,
                 const struct Alignment * const alignment, char * const Header,
                 const size_t SequenceLength, const float RAVE, const int N, const PrintInput_t * const extra)
{
    char * cptr = Header;
    while ( *cptr != ' ' && *cptr != '\0') ++cptr;
    *cptr = '\0';
    RawToNormalizedFunctionPtr RawToNormalizedFunction = prf->RawToNormalized;
    const float * const restrict NormCoefs = prf->NormalizationCoefs;
    for (int i=0; i<N; ++i) {
        const float normtest = (RawToNormalizedFunction == NULL)
                             ? 0.0f
                             : RawToNormalizedFunction(alignment[i].Score, NormCoefs, RAVE, SequenceLength);
        fprintf(stdout, "%s %i %f\n%s\n",
                Header,
                alignment[i].Score,
                normtest,
                &AlignedSequence[i][1]);
    }
}


// profile_ac start strop raw_score aln_string
void PrintTSV(const struct Profile * const prf, const char * * const AlignedSequence,
              const struct Alignment * const alignment, char * const Header,
              const size_t SequenceLength, const float RAVE, const int N, const PrintInput_t * const extra)
{
    char * cptr;
    // Extract Header
    const char * hptr = Header;
    while ( *hptr != ' ' && *hptr != '\t' && *hptr != '\0') ++hptr;
    int HeaderLength = (int) ((uintptr_t) hptr - (uintptr_t) Header);
    hptr = Header;
    if (*hptr == '>' || *hptr == '@') {
       hptr++;
       HeaderLength--;
    }   

    // Extract ID (id)
    char * id = calloc(1+strlen(prf->Identification),sizeof(char));
    strcpy(id,prf->Identification);
    cptr = id;
    while ( *cptr != ';' && *cptr != '\0') ++cptr;
    *cptr = '\0';

    char * ac  = calloc(strlen(prf->AC_Number),sizeof(char));
    strcpy(ac,prf->AC_Number);
    cptr = ac;
    while ( *cptr != ';' && *cptr != '\0' ) ++cptr;
    *cptr = '\0';
    RawToNormalizedFunctionPtr RawToNormalizedFunction = prf->RawToNormalized;
    const float * const restrict NormCoefs = prf->NormalizationCoefs;
		
		char strand;
		switch(extra->SeqId & 0xC0000000) {
			case 0x00000000: strand = '+'; break;
			case 0x40000000: strand = 'r'; break;
			case 0x80000000: strand = 'c'; break;
			case 0xC0000000: strand = '-'; break;
		}
// 		if ( extra->SeqId & 0xC0000000 ) { // FIXME: this does not cover correctly the 'c' cases ??? (MP)
    if( extra->SeqId == 0x40000000 || extra->SeqId == 0xC0000000 ) {
       for (int i=0; i<N; ++i) {
					const float normtest = (RawToNormalizedFunction == NULL)
															? 0.0f
															: RawToNormalizedFunction(alignment[i].Score, NormCoefs, RAVE, SequenceLength);
					fprintf(stdout, "%s|%s\t%i\t%i\t%.*s\t%i\t%i\t%i\t%f\t%c\t%s\n",
									ac,
									id,
                                    alignment[i].IPMB,
                                    alignment[i].IPME,
									HeaderLength,
									hptr,
									(int) SequenceLength - alignment[i].Region.Sequence.Begin + 1,
									(int) SequenceLength - alignment[i].Region.Sequence.End + 1,
									alignment[i].Score,
									normtest,
									strand,
									&AlignedSequence[i][1]);
			}
		}
		else {
			for (int i=0; i<N; ++i) {
					const float normtest = (RawToNormalizedFunction == NULL)
															? 0.0f
															: RawToNormalizedFunction(alignment[i].Score, NormCoefs, RAVE, SequenceLength);
					fprintf(stdout, "%s|%s\t%i\t%i\t%.*s\t%i\t%i\t%i\t%f\t%c\t%s\n",
									ac,
									id,
                                    alignment[i].IPMB,
                                    alignment[i].IPME,
									HeaderLength,
									hptr,
									alignment[i].Region.Sequence.Begin,
									alignment[i].Region.Sequence.End,
									alignment[i].Score,
									normtest,
									strand,
									&AlignedSequence[i][1]);
			}
		}
    
    
    free(ac);
    free(id);
}


void PrintDefault(const struct Profile * const prf, const char * * const AlignedSequence,
                  const struct Alignment * const alignment, char * const Header,
                  const size_t SequenceLength, const float RAVE, const int N, const PrintInput_t * const extra)
{
    char * des;
    char * cptr = Header;

    // Extract ID (id), short ID (cptr) and description (des)
    while ( *cptr != ' ' && *cptr != '\0') ++cptr;
    des = ( *cptr == '\0' ) ? cptr : cptr + 1;
    *cptr = '\0';
    while ( *cptr != '|' && *cptr != '>' && *cptr != '@') --cptr;
    ++cptr;

		// Extract Header
		const char * hptr = Header;
		if (*hptr == '>' || *hptr == '@') hptr++;
		const char * startHdr = hptr;
		while ( *hptr != ' ' && *hptr != '\t' && *hptr != '\0') hptr++;
		const int HeaderLength = (int) ((uintptr_t) hptr - (uintptr_t) startHdr);
		
		if (startHdr[0] == '\0') startHdr = Dflt;
		
		RawToNormalizedFunctionPtr RawToNormalizedFunction = prf->RawToNormalized;
    const float * const restrict NormCoefs = prf->NormalizationCoefs;
		
		if (extra->SeqId & 0xC0000000) {
			for (unsigned int i=0; i<N; ++i) {
					const float normtest = (RawToNormalizedFunction == NULL)
															? 0.0f
															: RawToNormalizedFunction(alignment[i].Score, NormCoefs, RAVE, SequenceLength);
					fprintf(stdout, "%8.3f %7i pos. %7i - %7i %.*s(-) %s\n",
									normtest,
									alignment[i].Score,
									(int) SequenceLength - alignment[i].Region.Sequence.Begin + 1,
									(int) SequenceLength - alignment[i].Region.Sequence.End + 1,
									HeaderLength, startHdr,
									des
								);
			}
		}
		else {
			for (unsigned int i=0; i<N; ++i) {
					const float normtest = (RawToNormalizedFunction == NULL)
															? 0.0f
															: RawToNormalizedFunction(alignment[i].Score, NormCoefs, RAVE, SequenceLength);
					fprintf(stdout, "%8.3f %7i pos. %7i - %7i %.*s %s\n",
									normtest,
									alignment[i].Score,
									alignment[i].Region.Sequence.Begin,
									alignment[i].Region.Sequence.End,
									HeaderLength, startHdr,
									des
								);
			}
		}
}


void PrintInterpro(const struct Profile * const prf, const char * * const AlignedSequence,
                   const struct Alignment * const alignment, char * const Header,
                   const size_t SequenceLength, const float RAVE, const int N, const PrintInput_t * const extra)
{ // replicates old pfsearch -lxz output
    char * id;
    char * des;
    char * buffer = calloc(OutputPrintWidth+1,sizeof(char));
    char * cptr = Header;
		
    // Extract ID (id), short ID (cptr) and description (des) of matched entry
    while ( *cptr != ' ' && *cptr != '\0') ++cptr;
    des = ( *cptr == '\0' ) ? cptr : cptr + 1;
    *cptr = '\0';
    while ( *cptr != '|' && *cptr != '>' && *cptr != '@') --cptr;
    ++cptr;
    id = Header;
    if (*Header == '>' || *cptr == '@') ++id;
		
    RawToNormalizedFunctionPtr RawToNormalizedFunction = prf->RawToNormalized;
    const float * const restrict NormCoefs = prf->NormalizationCoefs;
		const int level = prf->CutOffData.Values[ (int) prf->LevelIndex ].MCLE;
		
    for (unsigned int i=0; i<N; ++i) {
        const float normtest = (RawToNormalizedFunction == NULL)
                             ? 0.0f
                             : RawToNormalizedFunction(alignment[i].Score, NormCoefs, RAVE, SequenceLength);
        
        if (N > 1)
        {
            fprintf(stdout, ">%s_%i L=%i %.3f %6i pos. %8i -%8i [%5i, %5i] %s %s\n",
                cptr,
                i+1,
                level,
                normtest,
                alignment[i].Score,
                alignment[i].Region.Sequence.Begin,
                alignment[i].Region.Sequence.End,
                alignment[i].IPMB,
                alignment[i].IPME,
                id,
                des
               );
        }
        else
        {
            fprintf(stdout, ">%s L=%i %.3f %6i pos. %8i -%8i [%5i, %5i] %s %s\n",
                cptr,
                level,
                normtest,
                alignment[i].Score,
                alignment[i].Region.Sequence.Begin,
                alignment[i].Region.Sequence.End,
                alignment[i].IPMB,
                alignment[i].IPME,
                id,
                des
               );
        }

        *buffer = '\0';
        unsigned int offset = 0;
        unsigned int len    = strlen(&AlignedSequence[i][1]);

        while ( offset <= len ) {
            strncpy(buffer,&AlignedSequence[i][1]+offset,OutputPrintWidth);
            if ( '\0' != *buffer ) {
                fputs(buffer,stdout);
                fputc('\n',stdout);
            }
            offset += OutputPrintWidth;
        }
    }
    free(buffer);
}


void PrintPfscan(const struct Profile * const prf, const char * * const AlignedSequence,
                 const struct Alignment * const alignment, char * const Header,
                 const size_t SequenceLength, const float RAVE, const int N, const PrintInput_t * const extra)
{ // = ~old pfscan with -lxz option output (e.g. >DNA_primase_DnaG_arc L=0   32.064   9802 pos.        1 -     416 [    1,    -1] MF_00007|DNA primase DnaG [dnaG].\nMKYLIRAR... )
    char * buffer = calloc(OutputPrintWidth+1,sizeof(char));
    char * cptr;

    // Extract ID (id) of matching profile
    char * id = calloc(1+strlen(prf->Identification),sizeof(char));
    strcpy(id,prf->Identification);
    cptr = id;
    while ( *cptr != ';' && *cptr != '\0') ++cptr;
    *cptr = '\0';

    // Extract AC of matching profile
    char * ac  = calloc(1+strlen(prf->AC_Number),sizeof(char));
    strcpy(ac,prf->AC_Number);
    cptr = ac;
    while ( *cptr != ';' && *cptr != '\0' ) ++cptr;
    *cptr = '\0';
    RawToNormalizedFunctionPtr RawToNormalizedFunction = prf->RawToNormalized;
    const float * const restrict NormCoefs = prf->NormalizationCoefs;
    const char * des = prf->Description;

    for (unsigned int i=0; i<N; ++i) {
        const float normtest = (RawToNormalizedFunction == NULL)
                             ? 0.0f
                             : RawToNormalizedFunction(alignment[i].Score, NormCoefs, RAVE, SequenceLength);
        int level = -1;
        for ( unsigned int j = 0; j < MAXC; j++ ) { // find match level
            int icut = prf->CutOffData.Values[ j ].ICUT;
            int mcle = prf->CutOffData.Values[ j ].MCLE;
            if ( mcle == 0 && icut == 0 ) break;
            if ( alignment[i].Score >= icut ) { level = mcle; break; }
            level = mcle -1; // p.s. with -c option a match could have a score lower than the lowest defined level...
        } // p.s. prf->CutOffData.Values follows CUT_OFF line order in profile src; highest level should come first...

        if (N > 1)
        {
            fprintf(stdout, ">%s_%i L=%i %.3f %6i pos. %8i -%8i [%5i, %5i] %s|%s\n",
                id,
                i+1,
                level,
                normtest,
                alignment[i].Score,
                alignment[i].Region.Sequence.Begin,
                alignment[i].Region.Sequence.End,
                alignment[i].IPMB,
                alignment[i].IPME,
                ac,
                des
               );
        }
        else
        {
            fprintf(stdout, ">%s L=%i %.3f %6i pos. %8i -%8i [%5i, %5i] %s|%s\n",
                id,
                level,
                normtest,
                alignment[i].Score,
                alignment[i].Region.Sequence.Begin,
                alignment[i].Region.Sequence.End,
                alignment[i].IPMB,
                alignment[i].IPME,
                ac,
                des
               );
        } // p.s. as old pfscan could only scan one sequence at a time: here there is no reported seq id! so if used against
          //      multiple sequence output would be useless!

        *buffer = '\0';
        unsigned int offset = 0;
        unsigned int len    = strlen(&AlignedSequence[i][1]);

        while (offset <= len)
        {
            strncpy(buffer,&AlignedSequence[i][1]+offset,OutputPrintWidth);
            if ( '\0' != *buffer )
            {
                fputs(buffer,stdout);
                fputc('\n',stdout);
            }
            offset += OutputPrintWidth;
        }
    }
    free(buffer);
    free(id);
    free(ac);
}


void PrintClassification(const struct Profile * const prf, char * * const AlignedSequence,
                         const struct Alignment * const alignment, char * const Header,
                         const size_t SequenceLength, const float RAVE, const int N, const PrintInput_t * const extra)
{
    RawToNormalizedFunctionPtr RawToNormalizedFunction = prf->RawToNormalized;
    const float * const restrict NormCoefs = prf->NormalizationCoefs;
    const float normtest = (RawToNormalizedFunction == NULL)
                         ? 0.0f
                         : RawToNormalizedFunction(alignment[0].Score, NormCoefs, RAVE, SequenceLength);
    fprintf(stdout, "%s %i %f\n%s\n",
            Header,
            alignment[0].Score,
            normtest,
            AlignedSequence[0]);
    struct Profile * BestPrf = NULL;
    register struct Profile * tmpPrf = (struct Profile *) prf;
    int BestScore = NLOW;

    while (tmpPrf->next) {
        tmpPrf = tmpPrf->next;
        const int Score = RescoreAlignment(tmpPrf, &AlignedSequence[0][1], alignment, SequenceLength);
        fprintf(stdout,"%s : score %i\n", tmpPrf->Description, Score);
        if (BestScore < Score) {
            BestScore = Score;
            BestPrf = tmpPrf;
        }
    }

    fprintf(stdout, "Best : %s : score %i\n", BestPrf->Description, BestScore);
}


void PrintPsScan(const struct Profile * const prf, const char * * const AlignedSequence,
                 const struct Alignment * const alignment, char * const Header,
                 const size_t SequenceLength, const float RAVE, const int N, const PrintInput_t * const extra)
{
    char * id;
    char * buffer = calloc(OutputPrintWidth+1,sizeof(char));
    char * cptr = Header;

    RawToNormalizedFunctionPtr RawToNormalizedFunction = prf->RawToNormalized;
    const float * const restrict NormCoefs = prf->NormalizationCoefs;

    for (unsigned int i=0; i<N; ++i) {
        const float normtest = (RawToNormalizedFunction == NULL)
                             ? 0.0f
                             : RawToNormalizedFunction(alignment[i].Score, NormCoefs, RAVE, SequenceLength);
        int level = -1;
        for ( unsigned int j = 0; j < MAXC; j++ ) { // find match level
            int icut = prf->CutOffData.Values[ j ].ICUT;
            int mcle = prf->CutOffData.Values[ j ].MCLE;
            if ( mcle == 0 && icut == 0 ) break;
            if ( alignment[i].Score > icut ) { level = mcle; break; }
            level = mcle -1; // p.s. with -c option a match could have a score lower than the lowest defined level...
        } // p.s. prf->CutOffData.Values follows CUT_OFF line order in profile src; highest level should come first...

        {
            fprintf(stdout, ">%s L=%i %.3f %6i pos. %8i -%8i [%5i, %5i] %s|%s %s\n",
                    prf->Identification,
                    level,
                    normtest,
                    alignment[i].Score,
                    alignment[i].Region.Sequence.Begin,
                    alignment[i].Region.Sequence.End,
                    alignment[i].IPMB,
                    alignment[i].IPME,
                    prf->AC_Number,
                    prf->Identification,
                    prf->Description
                );
        }

        *buffer = '\0';
        unsigned int offset = 0;
        unsigned int len    = strlen(&AlignedSequence[i][1]);

        while (offset <= len)
        {
            strncpy(buffer,&AlignedSequence[i][1]+offset,OutputPrintWidth);
            if ( '\0' != *buffer )
            {
                fputs(buffer,stdout);
                fputc('\n',stdout);
            }
            offset += OutputPrintWidth;
        }
    }
    free(buffer);
}


void PrintIncmatch( const struct Profile * const prf, const char * * const AlignedSequence,
                    const struct Alignment * const alignment, char * const Header,
                    const size_t SequenceLength, const float RAVE, const int N, const PrintInput_t * const extra)
{ // replicates old pfsearch -zkx output
		const static char DefaultID[] = "unknown";
    const char * id;
    char * buffer = calloc(OutputPrintWidth+1,sizeof(char));
    char * cptr = Header;

    while ( *cptr != ' ' && *cptr != '\0' ) ++cptr;
    *cptr = '\0';
    id = Header;
		if (*Header == '>' || *cptr == '@') ++id;
		if (*id == '\0') id = DefaultID;

    RawToNormalizedFunctionPtr RawToNormalizedFunction = prf->RawToNormalized;
    const float * const restrict NormCoefs = prf->NormalizationCoefs;
    for ( unsigned int i=0; i<N; ++i ) {
        const float norm_score = (RawToNormalizedFunction == NULL)
                           ? 0.0f
                           : RawToNormalizedFunction(alignment[i].Score, NormCoefs, RAVE, SequenceLength);
													 
						// For compliance with fortran we add a extra end space
            fprintf(stdout, ">%s/%i-%i motif=%s|%s norm_score=%.3f raw_score=%i seq_end=%i motif_start=%i motif_end=%i \n",
                id,
                alignment[i].Region.Sequence.Begin,
                alignment[i].Region.Sequence.End,
                prf->AC_Number,
                prf->Identification,
                norm_score,
                alignment[i].Score,
								alignment[i].Region.Sequence.End - (int) SequenceLength - 1,
                alignment[i].IPMB,
                alignment[i].IPME
               );

        *buffer = '\0';
        unsigned int offset = 0U;
        unsigned int len    = strlen(&AlignedSequence[i][1]);

        while ( offset <= len ) {
            strncpy(buffer,&AlignedSequence[i][1]+offset,OutputPrintWidth);
            if ( '\0' != *buffer ) {
                fputs(buffer,stdout);
                fputc('\n',stdout);
            }
            offset += OutputPrintWidth;
        }
    }
    free(buffer);
}


void PrintPSMaker( const struct Profile * const prf, const char * * const AlignedSequence,
                   const struct Alignment * const alignment, char * const Header,
                   const size_t SequenceLength, const float RAVE, const int N, const PrintInput_t * const extra)
{ // replicates old pfsearch -x output
    char * buffer = calloc( OutputPrintWidth+1, sizeof( char ) );

    char * header = calloc( strlen( Header ) + 1, sizeof( char ) );
    strcpy( header, Header );

    char * cptr = Header;
    while ( *cptr != ' ' && *cptr != '\0') ++cptr;
    *cptr = '\0';
    while ( *cptr != '|' && *cptr != '>' && *cptr != '@') --cptr;
    ++cptr;

    RawToNormalizedFunctionPtr RawToNormalizedFunction = prf->RawToNormalized;
    const float * const restrict NormCoefs = prf->NormalizationCoefs;
    for ( unsigned int i=0; i<N; ++i ) {
        const float norm_score = (RawToNormalizedFunction == NULL)
                           ? 0.0f
                           : RawToNormalizedFunction(alignment[i].Score, NormCoefs, RAVE, SequenceLength);
            fprintf( stdout, ">%s %.3f %6i pos. %8i -%8i %s\n",
                cptr,
                norm_score,
                alignment[i].Score,
                alignment[i].Region.Sequence.Begin,
                alignment[i].Region.Sequence.End,
                header+1
            );

        *buffer = '\0';
        unsigned int offset = 0;
        unsigned int len    = strlen(&AlignedSequence[i][1]);

        while ( offset <= len ) {
            strncpy(buffer,&AlignedSequence[i][1]+offset,OutputPrintWidth);
            if ( '\0' != *buffer ) {
                fputs(buffer,stdout);
                fputc('\n',stdout);
            }
            offset += OutputPrintWidth;
        }
    }
    free(buffer);
    free(header);
}

void PrintxPSA( const struct Profile * const prf, const char * * const AlignedSequence,
								const struct Alignment * const alignment, char * const Header,
								const size_t SequenceLength, const float RAVE, const int N, const PrintInput_t * const extra)
{
	//PFSCANV3 NEWOUT   >fig|83333.1.peg.4317/36-219 motif=MF_00223|FolE norm_score=33.703 raw_score=4046 level_tag=! seq_end=-4 motif_start=1 motif_end=-1
	char * buffer = calloc(OutputPrintWidth+1,sizeof(char));
	char * cptr;
	
	// Extract Header
	const char * hptr = Header;
	if (*hptr == '>' || *hptr == '@') hptr++;
	const char * startHdr = hptr;
	while ( *hptr != ' ' && *hptr != '\t' && *hptr != '\0') hptr++;
	const int HeaderLength = (int) ((uintptr_t) hptr - (uintptr_t) startHdr);
	
	if (startHdr[0] == '\0') startHdr = Dflt;
	
	// Extract ID (id)
	char * id = calloc(1+strlen(prf->Identification),sizeof(char));
	strcpy(id,prf->Identification);
	cptr = id;
	while ( *cptr != ';' && *cptr != '\0') ++cptr;
	*cptr = '\0';
		    
  char * ac  = calloc(1+strlen(prf->AC_Number),sizeof(char));
	strcpy(ac,prf->AC_Number);
	cptr = ac;
	while ( *cptr != ';' && *cptr != '\0' ) ++cptr;
	*cptr = '\0';
  RawToNormalizedFunctionPtr RawToNormalizedFunction = prf->RawToNormalized;
	const float * const restrict NormCoefs = prf->NormalizationCoefs;
	const char * des = prf->Description;
	
	char strand;
	switch(extra->SeqId & 0xC0000000) {
		case 0x00000000: strand = '+'; break;
		case 0x40000000: strand = 'r'; break;
		case 0x80000000: strand = 'c'; break;
		case 0xC0000000: strand = '-'; break;
	}
	
	const int CurrentMode = prf->NormalizationData.Values[(int) prf->ModeIndex].NNOR;
#ifdef XALIT_DEBUG
	for (unsigned int i=0; i<N; ++i) {
	  fprintf(stderr,"%u/%d %u %s\n",i,N,AlignedSequence[i][0],AlignedSequence[i] + 1);
	}
#endif
	
	if (!(extra->SeqId & 0xC0000000)) {
		for (unsigned int i=0; i<N; ++i) {   
			const float normtest = (RawToNormalizedFunction == NULL) ? 0.0f : RawToNormalizedFunction(alignment[i].Score, NormCoefs, RAVE, SequenceLength);
			const static char DashTag[] = "#";
			const char * tag = DashTag;
			int BestFound = NLOW;
			for ( unsigned int j = 0; j < prf->CutOffData.JCUT; j++ ) {
				// Is this in the correct mode 
				_Bool Found = false;
				for (int k=0; k<prf->CutOffData.Values[j].JCNM; k++) {
					if (prf->CutOffData.Values[j].MCUT[k] == CurrentMode) {
						Found = true;
						break;
					}
				}
				if (Found) {
					const int icut = prf->CutOffData.Values[ j ].ICUT;
					if (icut > BestFound) {
						if ( alignment[i].Score >= icut) {
							BestFound = icut;
							if (prf->CutOffData.Values[j].CCUT[0] != '\0') tag = prf->CutOffData.Values[j].CCUT;
						}
					}
				}
			} 
		
			fprintf(stdout,	">%.*s/%i-%i motif=%s|%s norm_score=%.3f raw_score=%i level_tag=%s seq_end=%i motif_start=%i motif_end=%i motif_rev=%c strand=%c\n",
							HeaderLength,
							startHdr,
							alignment[i].Region.Sequence.Begin,
							alignment[i].Region.Sequence.End,
							ac,
							id,
							normtest,
							alignment[i].Score,
							tag,
							alignment[i].Region.Sequence.End - (int) SequenceLength, 
							alignment[i].IPMB,
							alignment[i].IPME,
							prf->isReversed ? (int) 'T' : (int) 'F', 
							(int) strand
			);

			*buffer = '\0';
			unsigned int offset = 0;
			const unsigned int len = strlen(&AlignedSequence[i][1]);

			while (offset <= len)
			{
				strncpy(buffer,&AlignedSequence[i][1]+offset,OutputPrintWidth);
				if ( '\0' != *buffer ) {
					fputs(buffer,stdout);
					fputc('\n',stdout);
				}
				offset += OutputPrintWidth;
			}
		}
	}
	else {
		for (unsigned int i=0; i<N; ++i) {   
			const float normtest = (RawToNormalizedFunction == NULL) ? 0.0f : RawToNormalizedFunction(alignment[i].Score, NormCoefs, RAVE, SequenceLength);
			const static char DashTag[] = "#";
			const char * tag = DashTag;
			int BestFound = NLOW;
			for ( unsigned int j = 0; j < prf->CutOffData.JCUT; j++ ) {
				// Is this in the correct mode 
				_Bool Found = false;
				for (int k=0; k<prf->CutOffData.Values[j].JCNM; k++) {
					if (prf->CutOffData.Values[j].MCUT[k] == CurrentMode) {
						Found = true;
						break;
					}
				}
				if (Found) {
					const int icut = prf->CutOffData.Values[ j ].ICUT;
					if (icut > BestFound) {
						if ( alignment[i].Score >= icut) {
							BestFound = icut;
							if (prf->CutOffData.Values[j].CCUT[0] != '\0') tag = prf->CutOffData.Values[j].CCUT;
						}
					}
				}
			} 
		
			fprintf(stdout,	">%.*s/%i-%i motif=%s|%s norm_score=%.3f raw_score=%i level_tag=%s seq_end=%i motif_start=%i motif_end=%i motif_rev=%c strand=%c\n",
							HeaderLength,
							startHdr,
							(int) SequenceLength - alignment[i].Region.Sequence.Begin + 1,
							(int) SequenceLength - alignment[i].Region.Sequence.End + 1,
							ac,
							id,
							normtest,
							alignment[i].Score,
							tag,
							alignment[i].Region.Sequence.End - (int) SequenceLength, 
							alignment[i].IPMB,
							alignment[i].IPME,
							prf->isReversed ? (int) 'T' : (int) 'F', 
							(int) strand
			);

			*buffer = '\0';
			unsigned int offset = 0;
			const unsigned int len = strlen(&AlignedSequence[i][1]);

			while (offset <= len)
			{
				strncpy(buffer,&AlignedSequence[i][1]+offset,OutputPrintWidth);
				if ( '\0' != *buffer ) {
					fputs(buffer,stdout);
					fputc('\n',stdout);
				}
				offset += OutputPrintWidth;
			}
		}
	}
	free(buffer);
	free(id);
	free(ac);
}

void PrintPfscanLOpt(const struct Profile * const prf, const char * * const AlignedSequence,
										 const struct Alignment * const alignment, char * const Header,
                     const size_t SequenceLength, const float RAVE, const int N, const PrintInput_t * const extra)
{ // replicates old pfscan with -l option output (e.g. L=0  32.064    9802 pos.       1 -     416 MF_00007|DNA primase DnaG [dnaG].)
  // could be used by ps_scan.pl (>=1.87)
    RawToNormalizedFunctionPtr RawToNormalizedFunction = prf->RawToNormalized;
    const float * const restrict NormCoefs = prf->NormalizationCoefs;

    for (unsigned int i=0; i<N; ++i) {
        const float normtest = (RawToNormalizedFunction == NULL)
                             ? 0.0f
                             : RawToNormalizedFunction(alignment[i].Score, NormCoefs, RAVE, SequenceLength);
        int level = -1;
        for ( unsigned int j = 0; j < MAXC; j++ ) { // find match level
            int icut = prf->CutOffData.Values[ j ].ICUT;
            int mcle = prf->CutOffData.Values[ j ].MCLE;
            if ( mcle == 0 && icut == 0 ) break;
            if ( alignment[i].Score >= icut ) { level = mcle; break; }
            level = mcle -1; // p.s. with -c option a match could have a score lower than the lowest defined level...
        } // p.s. prf->CutOffData.Values follows CUT_OFF line order in profile src; highest level should come first...

        fprintf(stdout, "L=%i %.3f %6i pos. %8i - %7i %s|%s %s\n",
                                level,
                                normtest,
                                alignment[i].Score,
                                alignment[i].Region.Sequence.Begin,
                                alignment[i].Region.Sequence.End,
                                prf->AC_Number,
                                prf->Identification,
                                prf->Description
                               );
    }
}

void PrintTurtle(const struct Profile * const prf, const char * * const AlignedSequence,
										 const struct Alignment * const alignment, char * const Header,
                     const size_t SequenceLength, const float RAVE, const int N, const PrintInput_t * const extra)
{ // HAMAP as SPARQL compliant output.
    RawToNormalizedFunctionPtr RawToNormalizedFunction = prf->RawToNormalized;
    const float * const restrict NormCoefs = prf->NormalizationCoefs;

    for (unsigned int i=0; i<N; ++i) {
        const float normtest = (RawToNormalizedFunction == NULL)
                             ? 0.0f
                             : RawToNormalizedFunction(alignment[i].Score, NormCoefs, RAVE, SequenceLength);
        int level = -1;
        for ( unsigned int j = 0; j < MAXC; j++ ) { // find match level
            int icut = prf->CutOffData.Values[ j ].ICUT;
            int mcle = prf->CutOffData.Values[ j ].MCLE;
            if ( mcle == 0 && icut == 0 ) break;
            if ( alignment[i].Score >= icut ) { level = mcle; break; }
            level = mcle -1; // p.s. with -c option a match could have a score lower than the lowest defined level...
        } // p.s. prf->CutOffData.Values follows CUT_OFF line order in profile src; highest level should come first...
        const char *firstpipe = strchr(Header+1, '|');
        const char * seqid;
        int _length;
        if (firstpipe != NULL) {
            const char *secondpipe = strchr(firstpipe+1, '|');
            if (secondpipe != NULL) {
                _length = secondpipe - (firstpipe+1);
                seqid=firstpipe+1;
            } else {
                seqid=firstpipe+1;
                _length = strlen(seqid);
            }
        } else {
            seqid=Header+1;
            _length = strlen(seqid);
        }
        fprintf(stdout, "yr:%.*s up:sequence ys:%.*s ;\n", _length, seqid, _length, seqid);
        fprintf(stdout, "  rdfs:seeAlso profile:%s .\n[ edam:is_output_of [\n  a edam:operation_0300 ;\n  edam:has_input profile:%s \n  ] ;\n",
                                prf->AC_Number,
                                prf->AC_Number
                               );
        fprintf(stdout, "  faldo:region [\n    faldo:begin [\n      faldo:position %d ;\n      faldo:reference ys:%.*s . ] ;\n    faldo:end [\n      faldo:position %d ;\n      faldo:reference ys:%.*s .] \n  ];\n  rdf:value",
                                alignment[i].Region.Sequence.Begin,
                                _length,
                                seqid,
                                alignment[i].Region.Sequence.End,
                                _length,
                                seqid);
        
        fprintf(stdout," '%s' . \n]\n", &AlignedSequence[i][1]);
    }
}

void PrintSAM(const struct Profile * const prf, const char * * const AlignedSequence,
              const struct Alignment * const alignment, char * const Header,
              const size_t SequenceLength, const float RAVE, const int N, const PrintInput_t * const extra)
{
	char buffer[512];
	Sequence SeqData = { .Data.Memory=NULL} ;
	PFSequence * PFSeq = &(SeqData.ProfileData);
	
// 	printf("\nFASTA file %s, sequence %lu, sequence len %u, sequence header len %u, seq offset %lu, quality len %u, quality header len %u, quality offset %lu\n",
// 				FASTQ->FASTA->FileName,
// 				FASTQ->SeqId,
// 				FASTQ->FASTA->DataPtr[FASTQ->SeqId].Sequence.SequenceLength,
// 				FASTQ->FASTA->DataPtr[FASTQ->SeqId].Sequence.HeaderLength,
// 				FASTQ->FASTA->DataPtr[FASTQ->SeqId].Sequence.Offset,
// 				FASTQ->FASTA->DataPtr[FASTQ->SeqId].Quality.SequenceLength,
// 				FASTQ->FASTA->DataPtr[FASTQ->SeqId].Quality.HeaderLength,
// 				FASTQ->FASTA->DataPtr[FASTQ->SeqId].Quality.Offset);
	
	for (int iAlign=0; iAlign<N; iAlign++) {
//  		printf("ALN  : %s\n", &AlignedSequence[iAlign][1]);
		char * restrict bufferptr = buffer;
		size_t len = 512;
		/* QNAME */
		{
			char * cptr = Header;
			if (*cptr == '>' || *cptr == '@') ++cptr;
			while ( *cptr == ' ' && *cptr != '\0') ++cptr;
			while ( *cptr != ' ' && *cptr != '\0') ++cptr;
			*cptr = '\0';
			while ( *cptr != '|' && *cptr != '>' && *cptr != '@') --cptr;
			++cptr;
			
			register const size_t written = snprintf(bufferptr, len, "%s\t", cptr);
			if (len <= written) goto bail2;
			len -= written;
			bufferptr += written;
		}
		
		/* FLAG */
		{
			const unsigned int value = (extra->SeqId & 0xC0000000) ? 18U : 2U;
			register const size_t written = snprintf(bufferptr, len, "%i\t", value);
			if (len <= written) goto bail2;
			len -= written;
			bufferptr += written;
		}
		
		/* RNAME */
		{
			const char * cptr = prf->AC_Number;
			while(*cptr != ';' && *cptr != '\0') { cptr++; }
			const int AClen = (int) ((uintptr_t) cptr - (uintptr_t) prf->AC_Number);
			register const size_t written = snprintf(bufferptr, len, "%.*s|%s\t", AClen, prf->AC_Number, prf->Identification);
			if (len <= written) goto bail2;
			len -= written;
			bufferptr += written;
		}
		
		/* Compute CIGAR */
		int InitialDeletion = 0;
		struct Cigar { unsigned int cnt; char state; } cigar[128];
		int From, To;
		{
			size_t index = 0UL;
			const char * restrict cptr = &AlignedSequence[iAlign][1];
					
			if (*cptr >= 'A' && *cptr <= 'Z')
				cigar[0].state = 'M';
			else if (*cptr == '-')
				cigar[0].state = 'D';
			else if (*cptr >= 'a' && *cptr <= 'z')
				cigar[0].state = 'I';
			else 
				cigar[0].state = '?';
			cigar[0].cnt = 1;
			cptr++;
			while (*cptr != '\0') {
				if (*cptr >= 'A' && *cptr <= 'Z') {
					if (cigar[index].state == 'M')
						cigar[index].cnt++;
					else {
						if (++index == 128) goto bail;
						cigar[index].state = 'M';
						cigar[index].cnt = 1;
					}
				}
				else if (*cptr == '-') {
					if (cigar[index].state == 'D')
						cigar[index].cnt++;
					else {
						if (++index == 128) goto bail;
						cigar[index].state = 'D';
						cigar[index].cnt = 1;
					}
				}
				else if (*cptr >= 'a' && *cptr <= 'z') {
					if (cigar[index].state == 'I')
						cigar[index].cnt++;
					else {
						if (++index == 128) goto bail;
						cigar[index].state = 'I';
						cigar[index].cnt = 1;
					}
				}
				else {
					if (cigar[index].state == '?')
						cigar[index].cnt++;
					else {
						if (++index == 128) goto bail;
						cigar[index].state = '?';
						cigar[index].cnt = 1;
					}
				}
				cptr++;
			}
			
			From = 0; 
			To = index;
			if (cigar[0].state == 'D') {
				InitialDeletion = cigar[0].cnt;
				From = 1;
			}
			
			if(cigar[index].state == 'D') {
				To = index - 1;
			}
		}
		
		/* POS */
		{
			register const size_t written = snprintf(bufferptr, len, "%i\t", alignment[iAlign].IPMB);
			if (len <= written) goto bail2;
			len -= written;
			bufferptr += written;
		}
				
		/* MAPQ */
		{
			register const size_t written = snprintf(bufferptr, len, "255\t");
			if (len <= written) goto bail2;
			len -= written;
			bufferptr += written;
		}
		
		/* CIGAR */
		{
			for (int i=From; i<=To; i++) {
				register const size_t written = snprintf(bufferptr, len, "%u%c", cigar[i].cnt, cigar[i].state);
				if (len <= written) goto bail2;
				len -= written;
				bufferptr += written;
			}
			if (len <= 1) goto bail2;
			*bufferptr++ = '\t';
			len -= 1;
		}
		
		/* RNEXT */
		{
			register const size_t written = snprintf(bufferptr, len, "*\t");
			if (len <= written) goto bail2;
			len -= written;
			bufferptr += written;
		}
		
		/* PNEXT */
		{
			register const size_t written = snprintf(bufferptr, len, "0\t");
			if (len <= written) goto bail2;
			len -= written;
			bufferptr += written;
		}
		
		/* TLEN */
		{
			register const size_t written = snprintf(bufferptr, len, "%lu\t", prf->Length - alignment[iAlign].IPMB + 1 + alignment[iAlign].IPME + 1);
			if (len <= written) goto bail2;
			len -= written;
			bufferptr += written;
		}
		
		/* SEQ */
		{
			const char * in = &AlignedSequence[iAlign][1];
			char * out = (char*) &AlignedSequence[iAlign][1];
			
			while (*in != '\0') {
				if (*in != '-') {
					if (*in >= 'a' && *in <= 'z')
						*out = *in - 'a' + 'A';
					else 
						*out = *in;
					out++;
				}
				in++;
			}
			*out = '\0';
		}
		
		/* QUAL */
		const size_t index = extra->SeqId & 0x3FFFFFFF;
		static unsigned char Star = '*';
		int QualLength;
		if (extra->DB->DataPtr[index].Quality.SequenceLength) {
			SeqData.Data.Memory = malloc(extra->DB->MaxSequenceSize*sizeof(unsigned char));
			if (SeqData.Data.Memory == NULL) goto bail1;
			SETUP_DATABASE_ACCESS(extra->SequenceFile);
			PFSeq = GET_DATABASE_SEQUENCE(&SeqData, &(extra->DB->DataPtr[index].Quality));
			if (extra->SeqId & 0x40000000) ReverseTranslatedSequence(PFSeq);
			UNSET_DATABASE_ACCESS();
			PFSeq->ProfileIndex += alignment[iAlign].Region.Sequence.Begin - 1;
			QualLength = (int) strlen(&AlignedSequence[iAlign][1]);	
		}
		/* Only treat sequence */
		else {
			QualLength = 1;
			PFSeq->ProfileIndex = &Star;
		}
		
		
		/* PRINT */
		{
			RawToNormalizedFunctionPtr RawToNormalizedFunction = prf->RawToNormalized;
			const float * const restrict NormCoefs = prf->NormalizationCoefs;
			const float NScore = (RawToNormalizedFunction == NULL)
                             ? 0.0f
                             : RawToNormalizedFunction(alignment[iAlign].Score, NormCoefs, RAVE, SequenceLength);
			
			printf("%s%s\t%.*s\tAS:i:%i\tNS:f:%.2f\n", buffer, &AlignedSequence[iAlign][1], QualLength,
						 PFSeq->ProfileIndex, alignment[iAlign].Score, NScore);
				
			
		}
	}
	
	if (SeqData.Data.Memory) free(SeqData.Data.Memory);
	
	return;
	bail2:
		fputs("Internal buffer of SAM format not big enough.\n", stderr);
		exit(1);
	bail1:
		fputs("Cannot allocate memory for sequence in SAM output.\n", stderr);
		exit(1);
	bail:
		fputs("The 128 space dedicated to cigar is not sufficient!\n", stderr);
		exit(1);
}
