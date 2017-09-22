#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <zlib.h>
#include <unistd.h>

#define LINESIZE 102400
// WM97_Chr04 //

// The speed of perl is acceptable if the column number is not too big, so there is no need to write c. 
// But when the column number is big, this c script will run faster. (For 135-columns SNP table, to process 1M lines takes 127 secs with perl, but only 21 secs with c.)

char *curr_time_str (char *s) ; 

char log_msg[LINESIZE * 4]; // Cannot understand why 'char *log_msg' doesn't work. 

// char *(argv[]) is equivalent to char *argv[] ; 
int main ( int argc, char *(argv[]) ) 
{
	// Initialize parameters. 
	int start_colN = 2; 

	if ( argc <= 1 && need_keyboard() ) {
		fprintf(stderr, "Not enough arguments. Should be \n\n"); 
		fprintf(stderr, "%s combined_1col.genotype [startColN_0-index] > combined_1col_diff.snp\n", argv[0]); 
		exit(1); 
	} else {
		sprintf(log_msg, "Begin to run [%s].\n", argv[0]); tsmsg(log_msg); 
	}
	
	FILE *inFp; 
	if ( need_keyboard() ) {
		// openFH( argv[1], inFp, "r" ); 
		if ( ( inFp = fopen( argv[1], "r" ) ) == NULL ) {
			perror("Error: Failed to open file.\n"); 
			exit(1); 
		}
		if ( argc >= 3 ) {
			start_colN = atoi( argv[2] ); 
		}
	} else {
		inFp = stdin; 
		if ( argc >= 2 ) {
			start_colN = atoi( argv[1] ); 
		}
	}

	sprintf(log_msg, "Starting 0-indexed col_num=[%d]\n", start_colN); tsmsg(log_msg); 

	
	char line[LINESIZE]; 
	char ref_str[LINESIZE] = "N"; 
	char raw_line[LINESIZE]; 
	//long int col_num = -2; 
	int isDiff = 0; 
	long int lineN = 0; 
	while ( fgets(line, LINESIZE, inFp) != NULL ) {

		lineN++; 
		rm_return(line); 

		// Reset variables 
		isDiff = 0; 
		strcpy(ref_str, "N"); 
		strcpy(raw_line, line); 
		
		//if ( col_num == -2 ) {
		//	col_num = get_col_num( line ); 
		//}
		
		// Check if this line contains differences. In function has_diff_genotype(), the string 'line' will be destroyed. 
		isDiff = has_diff_genotype(line, start_colN, ref_str); 

		if (isDiff == 1) {
			// printf("%d:%d:[%s] ref=%s\n", col_num, isDiff, raw_line, ref_str); 
			printf("%s\n", raw_line); 
		}
		
		if ( lineN % 1000000 == 1 ) {
			sprintf(log_msg, "[LOADING_SNPtbl] lineN=%d [%s]\n", lineN, raw_line); tsmsg(log_msg);
		}
	}
	fclose(inFp); 
	
	sprintf(log_msg, "All doen.\n"); tsmsg(log_msg); 
	return 0; 
}
//chr           pos     base
//WM97_Chr01    1       A
//WM97_Chr01    2       A

// Dynamic length array! 
// http://stackoverflow.com/questions/11198604/c-split-string-into-an-array-of-strings
// http://stackoverflow.com/questions/9210528/split-string-with-delimiters-in-c

// Here start_colN is 0-indexed. 
int has_diff_genotype ( char *line, int start_colN, char *refBase ) {
	int isDiff = 0; 
	char *part; 
	part = strtok(line, "\t"); 
	int i=0; 
	while ( part != NULL ) {
		if ( i < start_colN || strcmp(part, "N") == 0) {
			; 
		} else if ( strcmp(refBase, "N") == 0 ) {
			strcpy(refBase, part); 
		} else if ( strcmp(refBase, part) != 0 ) {
			isDiff = 1; 
			break; 
		}
		i++; 
		part = strtok(NULL, "\t"); 
	}
	return isDiff; 
}// has_diff_genotype() 

int get_col_num ( char in_line[] ) {
	char line[LINESIZE]; 
	strcpy(line, in_line); 
	long int count = 0; 
	char *part; 
	part = strtok(line, "\t"); 
	while ( part != NULL ) {
		count ++; 
		part = strtok(NULL, "\t"); 
	}
	return count; 
}

int need_keyboard () {
	// http://rosettacode.org/wiki/Check_input_device_is_a_terminal#C
	if ( isatty(fileno(stdin)) ) {
		// Need keyboard input, which means lack pipe-stdin or redirected-stdin ; 
		return 1; 
	} else {
		return 0; 
	}
}// need_keyboard() 

int openFH ( char *inFile, FILE *fh, char *type ) {
	if ( (fh = fopen( inFile, type ) ) == NULL ) {
		perror("Error: Failed to open file.\n"); 
		exit(1); 
	}
	return 0; 
}

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
	fprintf(stderr, "[%s] %s", time_txt, msg); 
	return 0; 
}


