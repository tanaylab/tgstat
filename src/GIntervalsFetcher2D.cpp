#include "GIntervalsFetcher2D.h"

bool GIntervalsFetcher2D::compare_for_sort(const GInterval2D &interv1, const GInterval2D &interv2)
{
	return interv1 < interv2;
}

