#ifndef HANDED_HPP
#define HANDED_HPP

#include "data.hpp"
#include <cmath>
#include <vector>
#include <algorithm>

#define MIN_DIST 1.0
#define MAX_DIST 1.0
#define MOMDIST 0		/* if zero skip momentum test */

void EventEta(event &, double *, unsigned *);

#endif	// HANDED_HPP
