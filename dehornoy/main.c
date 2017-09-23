/* From: "Valentina V. Shvachko" <valya@infos.botik.ru>
To: "Sergej Chmutov" <chmutov@univaq.it>
Subject: main.c
Date: Mon, 11 May 1998 11:23:12 +0300 */
/* Main for testing braids */
/* #include "..\Include\braid.h" */

/* #include "braid.h" */
#include "word.c"

int *Reduction();

void main() { 
	FILE *fp1, *fp2;
	int i = 0, j, Length1;

	fp1 = fopen( "wordin", "r" );
	if( fp1 == NULL ) {
		printf( "File doesn't exit.\n" );
		return;
	}
	while( fscanf( fp1, "%d", &Word[i]) != EOF )
		i++;
	if( i > 1000 )  {
		printf( "NO memory" );
		return;
	}
	Length1 = i - 1;

	i = 0;
	while( Word[i] !=0 ) 
		i++;
	LENGTH = i;
	Word[i] = (-1) * Word[Length1];
	for( i = LENGTH+1 , j = Length1 - 1; i <= j; i++, j--) {
		int c;
		c = Word[i];
		Word[i] = (-1) * Word[j];
		Word[j] = (-1) * c;
	}
	LENGTH = Length1 ;

	fp2 = fopen( "wordout", "w" );

	if( Reduction( ) == 0 ) {
		fprintf( fp2, "The reduced Word is empty\n" );
		printf( "The reduced Word is empty\n" );
	} else {
		fprintf( fp2, "The reduced Word is\n " );
		printf( "The reduced Word is\n " );

		for( j = 0; j < LENGTH; j++) {
			fprintf( fp2, "%d ", Word[j]);
			printf( "%d ", Word[j]);
		}
	}
	
	
	fclose( fp1);
	fclose( fp2);
	getchar();
	return;
}
	
/*	
#include "..\Include\braid.h"

int *Reduction();

void main() { 
	FILE *fp1, *fp2;
	int i = 0, j;

	fp1 = fopen( "test6.txt", "r" );
	if( fp1 == NULL ) {
		printf( "File doesn't exit.\n" );
		return;
	}
	while( fscanf( fp1, "%d", &Word[i]) != EOF )
		i++;
	if( i > 1000 )  {
		printf( "NO memory" );
		return;
	}
	LENGTH = i; 
	fp2 = fopen( "res.txt", "w" );

	if( Reduction( ) == 0 ) {
		fprintf( fp2, "The reduced Word is empty\n" );
		printf( "The reduced Word is empty\n" );
	} else {
		fprintf( fp2, "The reduced Word is\n " );
		printf( "The reduced Word is\n " );

		for( j = 0; j < LENGTH; j++) {
			fprintf( fp2, "%d ", Word[j]);
			printf( "%d ", Word[j]);
		}
	}
	
	
	fclose( fp1);
	fclose( fp2);
	getchar();
	return;
}
*/


















