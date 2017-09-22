#include <stdio.h>
#include <string.h>
#include <stdlib.h>

char *str_chop ( char **s, char delim, char *prev_delim, char **prev_s ) ; 

int main () {

	printf("char *str_chop( char **s, char delim, char *prev_delim, char **prev_s )\n"); 
	printf("int rm_return ( char *s )\n"); 
	printf("int need_keyboard ()\n"); 
	printf("int openFH ( char *inFile, FILE **fh, char *type )\n"); 


	return 0; 
}

// Usage : int openFH ( char *inFile, FILE **fh, char *type ) 
//   *inFile  : input file name. 
//   **fh     : &file_handle to receive result. 
//   *type    : string as "w"/"r" or something for fopen(); 
// Desc. : Open a file and set file handle value to *fh ; 
int openFH ( char *inFile, FILE **fh, char *type ) {
	if ( (*fh = fopen( inFile, type ) ) == NULL ) {
		perror("Error: Failed to open file.\n");
		exit(1);
	}
	return 0;
}


// Something to know is that if assign a (char *) pointer to a string, it will be having a 'const' attribution, so we may need to use strdup() to assign a manual string to it. 
// Example: 
//   char *s; s = strdup("some_string"); ... 

// Leaned and modifed from : http://rosettacode.org/wiki/Tokenize_a_string#C 
// Usage: char *str_chop ( char **s, char delim, char *prev_delim, char **prev_s 
//   **s   is for string_array; 
//   delim is for delimitation, such as : char delim = '\t'; 
//   *prev_delim should be the same character as delim at first, which is used for stream control. 
//   **prev_s is the sequence chopped out each time. 
// Example: 
//   char array[] = "1,2,3,,5,", *a1, **a2; a1=&(array[0]); a2=&a1; 
//   char delim=',', pd=',', *pdp; pdp = &pd; 
//
//   char *ps1, **ps2; ps1=array; ps2=&ps1; 
//   while ( str_chop( a2, delim, pdp, ps2 ) != NULL ) {
//   	printf("Element=|%s|\n", *ps2); 
//   }
char *str_chop ( char **s, char delim, char *prev_delim, char **prev_s ) {
	*prev_s = *s; // Previous start address in *s(tring) ; 

	if ( *prev_delim == '\0' ) {
		prev_delim = NULL; 
		return NULL; 
	}
	if ( **s == '\0' ) {
		// Get to the end. And *prev_delim is not '\0' or NULL, 
		//  so there is one more empty element to return; 
		// Here *prev_delim should equal to delim. 
		*prev_delim = '\0'; 
		return *prev_s; 
	}

	*prev_delim = delim; 
	while ( **s && (delim != **s) ) (*s)++; 
	**s = 0; 
	(*s)++; 

	return *prev_s; 
} //str_chop() 

// Usage : int rm_return ( char *s ) 
// Desc. : Remove the tailing continuous '\n's or '\r's until clean. 
//         s(tring) will be changed. 
int rm_return ( char *s ) {
	while ( s[strlen(s)-1] == '\n' || s[strlen(s)-1] == '\r' ) {
		s[strlen(s)-1] = '\0';
	}
	return 0;
}

// Usage : int need_keyboard() 
// Desc. : Check the stdin file pointer, and return 
//          1 - There is no pipe-stdin or redirected-stdin ; This means there should be a file in argv[] as input. 
//          0 - There is pipe-stdin or redirected-stdin input. 
//         This is usually used for checking if there is any input from file or stdin in a pipeline. 
int need_keyboard () {
	// http://rosettacode.org/wiki/Check_input_device_is_a_terminal#C
	if ( isatty(fileno(stdin)) ) {
		// Need keyboard input, which means lack pipe-stdin or redirected-stdin ;
		return 1;
	} else {
		return 0;
	}
}// need_keyboard()


