#include "util.h"

double volumeAtIndex(int i)
{
	// Returns the volume (cm^3) of volume element i
	if (i >= 0 && i < POINTS)
	{
		if (i == 0 || i == 1)
			return SUBSURFACE_VOLUME; // index 0 should have 0 volume, but make it a nonzero volume to not break anything
		if (i == FIRST_BULK_INDEX)
			return FIRST_BULK_VOLUME;
		if (i < FIRST_ELONGATED_INDEX)
			return VOLUME;
		return ELONGATED_VOLUME;

	}

	return VOLUME; // shouldn't get here
}