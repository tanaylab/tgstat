#include "port.h"
BASE_CC_FILE

#include "GInterval.h"

//------------------------------------- GInterval -----------------------------------------------

const char *GInterval::COL_NAMES[GInterval::NUM_COLS] = { "chrom", "start", "end" };

double GInterval::dist2coord(int64_t coord, double margin) const
{
	double res = 0;
	double left_dist = strand == 1 ? coord - start : start - coord;
	double right_dist = strand == 1 ? coord - end : end - coord;

	if (!margin) {
		// is coord inside interv?
		if (coord >= start && coord <= end)
			return 0;

		res = fabs(left_dist) <= fabs(right_dist) ? left_dist : right_dist;
	} else {
		// is coord inside interv?
		if (coord >= start && coord <= end)
			res = (margin * (left_dist + right_dist)) / (double)(end - start);
		else {
			double offset = strand == 1 ? margin : -margin;
			res = fabs(left_dist) <= fabs(right_dist) ? left_dist - offset : right_dist + offset;
		}
	}
	return strand ? res : fabs(res);
}

int64_t GInterval::dist2interv(const GInterval &interv, bool touch_is_at_dist_one) const
{
	// do interv1 and interv2 overlap?
	if (max(start, interv.start) < min(end, interv.end))
		return 0;

	int64_t left_dist = (interv.strand == 1 ? -1 : 1) * (interv.start - end + touch_is_at_dist_one);
	int64_t right_dist = (interv.strand == 1 ? -1 : 1) * (interv.end - start - touch_is_at_dist_one);
	int64_t res = llabs(left_dist) <= llabs(right_dist) ? left_dist : right_dist;
	return interv.strand ? res : llabs(res);
}

double GInterval::coverage_ratio(const GInterval &interv) const
{
	double intersect = min(end, interv.end) - max(start, interv.start);

	if (intersect <= 0)
		return 0;

	double uni = max(end, interv.end) - min(start, interv.start);

	return intersect / uni;
}

char GInterval::char2strand(char c)
{
	if (c == '+')
		return 1;
	if(c == '-')
		return -1;
	TGLError<GInterval>(BAD_STRAND, "Bad strand character %c", c);
	return 0;
}

