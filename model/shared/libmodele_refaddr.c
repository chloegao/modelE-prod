#include <stdio.h>

// Produce a reference address that can be used to match symbols to stacktrace
void libmodele_refaddr(void)
{
	printf("REFERENCE_ADDRESS libmodele_refaddr %p\n", libmodele_refaddr);
}
