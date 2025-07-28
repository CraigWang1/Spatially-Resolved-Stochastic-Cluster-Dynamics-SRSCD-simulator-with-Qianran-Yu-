#include "util.h"

double volumeAtIndex(int i)
{
	// Returns the volume (cm^3) of volume element i
	if (i >= 0 && i < POINTS)
	{
		if (i == 0)
			return SUBSURFACE_VOLUME; // index 0 should have 0 volume, but make it a nonzero volume to not break anything
		return DIVIDING_AREA * thicknessAtIndex(i);
	}

	return DIVIDING_AREA * FIRST_BULK_THICKNESS * NM_TO_CM; // return something nonzero to not break things
}

double thicknessAtIndex(int i)
{
	/*
	 * Returns the thickness of the mesh element at index i.
	 */
	if (i == 0)
		return 0;
	return FIRST_BULK_THICKNESS * pow(THICKNESS_SCALING_FACTOR, i - 1) * NM_TO_CM;
}

double lengthFrontAtIndex(int i)
{
	/* 
	 * Returns the distance between this mesh element's centroid and that of
	 * the element directly in front, with an exponential meshing scheme.
	 */
	if (i == 0)
		return 0;
	if (i == 1)
		return FIRST_BULK_THICKNESS / 2. * NM_TO_CM;
	return (thicknessAtIndex(i - 1) + thicknessAtIndex(i)) / 2.;
}

double lengthBackAtIndex(int i)
{
	/* 
	 * Returns the distance between this mesh element's centroid and that of
	 * the element directly behind, with an exponential meshing scheme.
	 */
	if (i == 0)
		return FIRST_BULK_THICKNESS / 2. * NM_TO_CM;
	return (thicknessAtIndex(i) + thicknessAtIndex(i + 1)) / 2.;
}