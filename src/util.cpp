#include "util.h"

double volumeAtIndex(int i)
{
	// Returns the volume (cm^3) of volume element i
	return length(i) * DIVIDING_AREA;
}

/* Returns the length (cm) of volume element i */
double length(int i)
{
	if (i == SURFACE_INDEX)
		return 0;
	if (i == SUBSURFACE_INDEX)
		return SUBSURFACE_THICKNESS * NM_TO_CM;
	if (i == FIRST_BULK_INDEX)
		return FIRST_BULK_THICKNESS * NM_TO_CM;
	if (i < FIRST_EXP_INDEX)
		return ELEMENT_THICKNESS * NM_TO_CM;
	return ELEMENT_THICKNESS * pow(EXP_LENGTH_MULT, i - FIRST_EXP_INDEX) * NM_TO_CM;
}

/* Returns the length (cm) between centers of indices i and i-1 */
double lengthf(int i)
{
	return (length(i) + length(i - 1))/2.;
}

/* Returns the length (cm) between centers of indices i and i+1 */
double lengthb(int i)
{
	return (length(i) + length(i + 1))/2.;
}