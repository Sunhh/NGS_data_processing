#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <zlib.h>

#define LINESIZE 1024
#define COLSIZE 50
#define IDLEN 20
#define NUMLEN 10
// WM97_Chr04 //

char *curr_time_str (char *s) ; 

struct refChr {
	// The string length must be in fixed length if I want to use malloc amd p++/p-- to operate the data. 
	char chrID[IDLEN]; 
	long int position;
	char refBase[2]; 
	int changeChr; 
}; 

struct colSNP {
	char colBase[COLSIZE]; 
	int isDiff; 
}; 

char log_msg[LINESIZE * 4]; // Cannot understand why 'char *log_msg' doesn't work. 

// char *(argv[]) is equivalent to char *argv[] ; 
int main ( int argc, char *(argv[]) ) 
{

	FILE *ipRefChr; 

	if (argc <= 3) {
		printf("Not enough arguments. Should be \n\n"); 
		printf("%s max_line_num WM97_v6.chrom.fa.refChr GS122_C.1col ...\n", argv[0]); 
		exit(1); 
	} else if (argc >= 4) {
		// ipRefChr = fopen( argv[1], "r"); 
		sprintf(log_msg, "Begin to run [%s].\n", argv[0]); tsmsg(log_msg); 
	} else {
		printf("What happens!\n"); 
		exit(1); 
	}

	// Load refChr file. 
	// struct refChr *rr = (struct refChr *) calloc ( argv[1] , sizeof(struct refChr) ); // The initialization step needs some time. 
//	struct refChr *(rr[100]); 
	struct refChr *rr; 
	rr = (struct refChr *) malloc( atol(argv[1]) * sizeof(struct refChr) ); 
	char *rr_headline; 
	rr_headline = (char *) malloc( 100 * sizeof(char) ); 
	long int rr_lineN; 

	struct tm *curr_time; 
	sprintf(log_msg, "Processing RefChr file [%s]\n", argv[2]); 
	rr_lineN = load_RefChr( argv[2], rr, rr_headline); 

	int argv_i ; 
	for (argv_i=3; argv_i<argc; argv_i++) {
		sprintf(log_msg, "Processing .1col file [%s]\n", argv[argv_i]); tsmsg(log_msg); 
		// Load in .1col files. 
		struct colSNP *cc; 
		cc = (struct colSNP *) malloc( atol(argv[1]) * sizeof(struct colSNP) ); 
		char *cc_headline; 
		cc_headline = (char *) malloc( COLSIZE * sizeof(char) ); 
		long int cc_lineN; 
		cc_lineN = load_1col( argv[argv_i], cc, cc_headline ); 
		
		// Do sth. with cc. 
		int within_dist = 5; 
		update_colSNP( rr , cc, rr_lineN, within_dist ); 
		char outfile[LINESIZE]; 
		strcpy( outfile, argv[argv_i] ); 
		strcat( outfile, ".maskClose"); 
		output_colSNP( cc , outfile, cc_lineN, cc_headline ); 
	
		// Free cc's memory. 
		free(cc); 

		sprintf(log_msg, "Got new file [%s]\n", outfile); tsmsg(log_msg); 
	}

	sprintf(log_msg, "All doen.\n"); tsmsg(log_msg); 
	return 0; 
}
//chr           pos     base
//WM97_Chr01    1       A
//WM97_Chr01    2       A

int output_colSNP ( struct colSNP *cc, char *outFile, long int max_lineN, char *headline ) {
	FILE *opFH; 
	if ( (opFH = fopen( outFile, "w" )) == NULL ) {
		perror( "Failed to open file" ); 
		exit(1); 
	}
	
	fprintf(opFH, "%s\n", headline); 
	long int count = 0; 
	while ( count < max_lineN ) {
		fprintf(opFH, "%s\n", cc->colBase); 
		count++; 
		cc++; 
	}

	fclose(opFH); 

	return 0; 
}


int update_colSNP ( struct refChr *rr, struct colSNP *cc, long int max_lineN, int within ) {
	long int count = 0, cc_start=-1, cc_end=-1; 
	while ( count < max_lineN ) {
		// Set cc->isDiff ; 
		if ( strcmp(rr->refBase, "N") == 0 || strcmp(cc->colBase, "N") == 0 || strcmp(rr->refBase, cc->colBase) == 0 ) {
			cc->isDiff = 0; 
		} else {
			cc->isDiff = 1; 
		}

		// Mask the colBase. 
		mask_colSNP( rr, cc, within, count, max_lineN ); 

		// Go to next; 
		count++; 
		rr++; 
		cc++; 
	}
	return count; 
}
int mask_colSNP ( struct refChr *rr, struct colSNP *cc, int within, long int cur_lineN, long int max_lineN ) {
	if ( within > cur_lineN+1 ) {
		within = cur_lineN+1; // _lineN is 0-based, while within is 1-based. 
		// within can be also described as 'the bases left for checking'. 
	}
	if ( cur_lineN == 0 || cur_lineN > max_lineN) {
		return 0; 
	}
	int diff_in_chk = 0, diff_high_i=-1, diff_low_i=within+1; 
	// Check from bigger to smaller positions. 
	// I should skip the checking if the biggest position is not changed. 
	if ( cc->isDiff != 1 ) {
		return 0; 
	}
	while ( within > 0 ) {
		if (cc->isDiff == 1) {
			diff_in_chk++; 
			if ( diff_high_i < within ) {
				diff_high_i = within; 
			}
			if ( diff_low_i > within ) {
				diff_low_i = within; 
			}
		}

		if ( rr->changeChr == 1 ) {
			// This line is the last line to be checked. The lines before should be ignored, and checked before. 
			break; 
		}

		within--; 
		rr--; 
		cc--; 
	}

	if ( diff_in_chk >= 2 ) {
		// Change the bases to N between the two SNPs. 
		while ( within <= diff_high_i ) {
			if ( within >= diff_low_i ) {
				strcpy(cc->colBase, "N"); 
			}

			within++; 
			cc++; 
		}
	}

	return 0; 
}// mask_colSNP() 



int load_1col ( char *inFile, struct colSNP *store, char *headline ) {
	FILE *ipFH; 
	if ( (ipFH = fopen( inFile, "r" )) == NULL ) {
		perror( "Failed to open file" ); 
		exit(1); 
	}
	
	char line[COLSIZE]; 
	long int lineN = 0; 

	while ( (fgets(line, COLSIZE, ipFH)) ) {
		// sscanf() doesn't work well, so I turn to use strtok() instead. 
		lineN++; 

		rm_return(line);   // Remove the tailing '\n' / '\r'; 

		// Skip the header line. 
		if (lineN == 1) {
			strcpy(headline, line); 
			continue; 
		}

		strcpy( (*store).colBase, line ); 

		// Output some processing information for control. 
		if ( lineN % 10000000 == 1 ) {
			sprintf(log_msg, "[LOADING_1COL] [%s] lineN=%d [%s]\n", inFile, lineN, store->colBase); tsmsg(log_msg); 
		}

		store++; 
	}
	lineN--; 
	fclose(ipFH); 

	return lineN; 
}

int load_RefChr ( char *inFile, struct refChr *store, char *headline ) {
	FILE *ipFH; 
	if ( (ipFH = fopen( inFile, "r" )) == NULL ) {
		perror( "Failed to open file" ); 
		exit(1); 
	}

	char line[LINESIZE]; 
	long int lineN = 0; 

	struct refChr *prev_store; 

	while ( (fgets(line, LINESIZE, ipFH)) ) {
		// sscanf() doesn't work well, so I turn to use strtok() instead. 
		lineN++; 

		rm_return(line);   // Remove the tailing '\n' / '\r'; 

		// Skip the header line. 
		if (lineN == 1) {
			strcpy(headline, line); 
			continue; 
		}
		
		// Read in and split each cell into structure. 
		char *part; 
		part = strtok(line, "\t"); 
		int i = 0; 
		while ( part != NULL && i < 3 ) {
			switch (i) {
				case 0: strcpy( (*store).chrID    , part); break; 
				case 1: (*store).position = atol((char *)part); break; 
				case 2: strcpy( (*store).refBase  , part); break; 
				default: break; 
			}
			part = strtok(NULL, "\t"); 
			i++; 
		}
		if ( lineN == 2 ) {
			store->changeChr = 1; 
			prev_store = store; 
		} else if ( strcmp( prev_store->chrID, store->chrID ) != 0 ) {
			store->changeChr = 1; 
			prev_store=store; 
		}

		// Output some processing information for control. 
		if ( lineN % 1000000 == 1 ) {
			sprintf(log_msg, "[LOADING_RefChr] [%s] lineN=%d [%s %d %s]\n", inFile, lineN, store->chrID, store->position, store->refBase); tsmsg(log_msg); 
		}

		store++; 
	}
	lineN--; 
	fclose(ipFH); 

	return lineN; 
}// load_RefChr() 

char *curr_time_str (char *s) {
	time_t rawtime; 
	struct tm * timeinfo; 
	time( &rawtime ); 
	timeinfo = localtime( &rawtime ); 
	strcpy( s, asctime( timeinfo ) ); 
	rm_return(s); 
	return s; 
}

int rm_return ( char *s ) {
	while ( s[strlen(s)-1] == '\n' || s[strlen(s)-1] == '\r' ) {
		s[strlen(s)-1] = '\0'; 
	}
	return 0; 
}

int tsmsg ( char *msg ) {
	char *time_txt; 
	time_txt = (char *) malloc ( 40 * sizeof(char) ); 
	curr_time_str( time_txt ); 
	printf("[%s] %s", time_txt, msg); 
	return 0; 
}


