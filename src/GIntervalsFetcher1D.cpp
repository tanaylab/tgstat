#include "GIntervalsFetcher1D.h"

bool GIntervalsFetcher1D::compare_by_start_coord(const GInterval &interv1, const GInterval &interv2)
{
	return interv1 < interv2;
}

bool GIntervalsFetcher1D::compare_by_end_coord(const GInterval &interv1, const GInterval &interv2)
{
	return interv1.chromid < interv2.chromid || interv1.chromid == interv2.chromid && interv1.end < interv2.end;
}

