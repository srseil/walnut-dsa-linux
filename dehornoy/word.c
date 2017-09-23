/* From: "Valentina V. Shvachko" <valya@infos.botik.ru>
To: "Sergej Chmutov" <chmutov@univaq.it>
Subject: Re: Kosy
Date: Thu, 14 May 1998 23:08:54 +0300 

Ispravlennyj fail word.c
On pochemu-to ne doshel do tebya.

Valya */

#include <string.h>
/*#include "..\Include\braid.h" */

#include "braid.h"

int Abs( int number){
	return (number < 0)? -number : number;
}

void PrintWord() {
	int i;
	for( i=0; i < LENGTH; i++ ) 
		printf( "%d ", Word[i] );
	printf("  ---\n");
}


int FindMainGenerator( ) {

	int i, j = Abs( Word[0] );
	for( i = 1; i < LENGTH; i++ ) 
		j = ( j < Abs(Word[i]) )? j : Abs( Word[i] );
	return j;
}




void FreeReduction( ) {
	int r, l;
	int Wl;

	for( l= -1, Wl=0, r=0; r < LENGTH; r++ ) {
		if( Wl == (-1)*Word[r] ) {
			Word[l] = Word[r] = 0;
			l--;
			Wl = ( l < 0 ? 0 : Word[l] );
		} else {
			l++;
			if( l != r ) {
				Word[l] = Word[r];
				Word[r] = 0;
			}
			Wl = Word[l];
		}
	}
	LENGTH = l+1;
}


int FindGenerator( int Gen, int N1, int N2 ) {
	int i;
	for( i = N1; i < N2; i++ ) {
		if( Abs( Word[i]) == Gen )
			return i;
	}
	return -1;
}




struct handle FindHandle( int Gen, int L1, int L2 ) {
	int h1, h2; /* this are positions the handle, that we find */
	struct handle H;
	h1 = FindGenerator( Gen, L1, L2 );
	if( h1 == -1 || h1 == L2 - 1) {
		H.generator = 0;
		return H;
	}

	while( ( h2 = FindGenerator( Gen, h1 + 1, L2 ) ) != -1) {
		if( Word[h1] == Word[h2])
			h1 = h2;
		else {
			H.generator = Gen;
			H.pos1 = h1;
			H.pos2 = h2;
			return H;
		}
	}
	H.generator = 0;
	return H;
}


int HandleTransform( struct handle H ) {
	int i,  e = Word[H.pos1]/Abs( Word[H.pos1] );
	int Length1, c = 0;
	int Temp[1000];
	

	for( i = H.pos1 + 1; i < H.pos2; i++ ) {
		if( i + c + 1 < 1001) {
			if( Abs( Word[i] ) != H.generator + 1 )
				Temp[i + 2*c - 1] = Word[i];
			else {
				Temp[i + 2*c - 1] = (-1) * e * (H.generator + 1);
				Temp[i + 2*c] = Word[i]/Abs(Word[i])*H.generator;
				Temp[i + 2*c + 1] = e*(H.generator + 1);
				c++;
			}
		} else
			return 1000;

	}

	memmove( Word + H.pos2 + 2*c-1,
			 Word + H.pos2 + 1,
			 (LENGTH - 1 - H.pos2) * sizeof(int) );

	memmove( Word + H.pos1,
			 Temp + H.pos1,
			 (H.pos2 + 2*c - 1 - H.pos1) * sizeof(int) );

	/*
	for( i = LENGTH - 1; i > H.pos2; i--) 
		Word[i + 2*c - 2] = Word[i];

	for( i = H.pos1;  i < H.pos2 + 2*c -1; i++ ) 
		Word[i] = Temp[i];
	*/

	LENGTH = LENGTH + 2*c - 2;
	if( LENGTH > 1000 ) {
		printf( "No memory" );
		return 1000;
	}
	Length1 = LENGTH;

/*	PrintWord();
	printf("---\n");*/

	FreeReduction();

	return 2 *c - 2 - ( Length1 - LENGTH);
}


/* Returns the the difference of old and new lenghts of the Word */
int HandleReduction( struct handle H ) {
	struct handle H1;
	int Pos2Change = 0, k, HTr;
	while( 1 ) {
		H1 = FindHandle( H.generator + 1, H.pos1 + 1, H.pos2 );
		if( H1.generator == 0 ) {
			if( ( HTr = HandleTransform( H ) ) == 1000 ) {
				printf( "No memory" );
				return 1000;
			}
			return Pos2Change + HTr;
		}
		k = HandleReduction( H1 );
		Pos2Change += k;
		H.pos2 += k;
		if ( H.pos2 < H.pos1 )
			return Pos2Change +HTr;
	}
	return 0;
}


int *Reduction( ) {
	int J;
	struct handle H;
	FreeReduction( );
	if( LENGTH == 0 ) 
		return  0;

	J = FindMainGenerator( );
	H = FindHandle( J, 0, LENGTH );
	if( H.generator == 0 )
		return Word;

	HandleReduction( H );
	Reduction( );
}










