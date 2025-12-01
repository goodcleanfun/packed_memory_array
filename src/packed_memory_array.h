#ifndef PMA_H
#define PMA_H

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdarg.h>
#include <stdint.h>
#include <float.h>
#include <math.h>

#include "bit_utils/bit_utils.h"

static const double PMA_DENSITY_HIGH_ROOT = 0.75;
static const double PMA_DENSITY_HIGH_LEAF = 1.0;

static const double PMA_DENSITY_LOW_ROOT = 0.5;
static const double PMA_DENSITY_LOW_LEAF = 0.25;

#define PMA_LARGEST_MAX_SPARSITY (uint8_t)(1 / PMA_DENSITY_LOW_LEAF)
#define PMA_LARGEST_EMPTY_SEGMENT (PMA_LARGEST_MAX_SPARSITY)

static const uint64_t PMA_MAX_SIZE = 1ULL << 56; // Maximum size of the array

#endif // PMA_H

#ifndef PMA_NAME
#error "Must define PMA_NAME"
#endif

#ifndef PMA_TYPE
#error "Must define PMA_TYPE"
#endif

#ifndef PMA_LESS_THAN
#define PMA_LESS_THAN_DEFINED
#define PMA_LESS_THAN(a, b) ((a) < (b))
#endif

#ifndef PMA_EQUALS
#define PMA_EQUALS_DEFINED
#define PMA_EQUALS(a, b) ((a) == (b))
#endif

#ifndef PMA_EMPTY_VALUE
#define PMA_EMPTY_VALUE_DEFINED
#define PMA_EMPTY_VALUE ((PMA_TYPE)0)
#endif

#define PMA_CONCAT_(a, b) a ## b
#define PMA_CONCAT(a, b) PMA_CONCAT_(a, b)
#define PMA_TYPED(name) PMA_CONCAT(PMA_NAME, _##name)
#define PMA_FUNC(func) PMA_CONCAT(PMA_NAME, _##func)

#ifndef PMA_MALLOC
#define PMA_MALLOC malloc
#define PMA_MALLOC_DEFINED
#endif

#ifndef PMA_FREE
#define PMA_FREE free
#define PMA_FREE_DEFINED
#endif

#ifndef PMA_ARRAY_NAME
#define PMA_ARRAY_NAME PMA_TYPED(store)
#define PMA_ARRAY_DEFINED
#define ARRAY_NAME PMA_ARRAY_NAME
#define ARRAY_TYPE PMA_TYPE
#ifndef PMA_CONCURRENT
#include "aligned_array/aligned_array.h"
#else
#include "concurrent_array/concurrent_array.h"
#endif
#undef ARRAY_NAME
#undef ARRAY_TYPE
#endif

#define PMA_ARRAY_FUNC(func) PMA_CONCAT(PMA_ARRAY_NAME, _##func)

typedef struct {
    PMA_ARRAY_NAME array;
    size_t segment_size;
    size_t num_segments;
    size_t count;
    uint8_t height; // Height of the tree;
    double density_high; // Density threshold for high density
    double density_low; // Density threshold for low density
} PMA_NAME;

static bool PMA_FUNC(init_size)(PMA_NAME *pma, size_t size) {
    if (pma == NULL) return false;
    pma->segment_size = PMA_LARGEST_EMPTY_SEGMENT;
    size_t cap = (1ULL << pma->segment_size);
    if (cap < size) {
        size = cap;
    }
    PMA_ARRAY_FUNC(init_size_fixed)(&pma->array, size);

    pma->count = 0;
    pma->num_segments = cap / pma->segment_size;
    pma->height = floor_log2(pma->num_segments) + 1;
    pma->density_high = (PMA_DENSITY_HIGH_LEAF - PMA_DENSITY_HIGH_ROOT) / pma->height;
    pma->density_low = (PMA_DENSITY_LOW_ROOT - PMA_DENSITY_LOW_LEAF) / pma->height;
    return true;
}

static PMA_NAME *PMA_FUNC(new_size)(size_t size) {
    PMA_NAME *pma = PMA_MALLOC(sizeof(PMA_NAME));
    if (pma == NULL) {
        return NULL;
    }
    if (!PMA_FUNC(init_size)(pma, size)) {
        return NULL;
    }
    return pma;
}

static PMA_NAME *PMA_FUNC(new)(void) {
    return PMA_FUNC(new_size)(1ULL << PMA_LARGEST_EMPTY_SEGMENT);
}

static inline bool PMA_FUNC(empty_at_unchecked)(PMA_NAME *pma, size_t index) {
    return PMA_EQUALS(PMA_ARRAY_FUNC(get_unchecked)(&pma->array, index), PMA_EMPTY_VALUE);
}

static inline bool PMA_FUNC(empty_at)(PMA_NAME *pma, size_t index) {
    PMA_TYPE value;
    return pma != NULL && PMA_ARRAY_FUNC(get)(&pma->array, index, &value) && PMA_EQUALS(value, PMA_EMPTY_VALUE);
}


static bool PMA_FUNC(find)(PMA_NAME *pma, PMA_TYPE key, int64_t *index) {
    if (pma == NULL || index == NULL) {
        return false;
    }
    PMA_ARRAY_NAME *array = &pma->array;
    int64_t lo = 0;
    int64_t hi = PMA_ARRAY_FUNC(size)(array) - 1;
    while (lo <= hi) {
        int64_t mid = (lo + hi) / 2;
        int64_t i = mid;
        while (i >= lo && PMA_FUNC(empty_at_unchecked)(pma,i)) {
            i--;
        }
        if (i < lo) {
            lo = mid + 1;
        } else {
            PMA_TYPE array_key = PMA_ARRAY_FUNC(get_unchecked)(array, i);
            if (PMA_LESS_THAN(array_key, key)) {
                lo = mid + 1;
            } else if (PMA_EQUALS(array_key, key)) {
                *index = i;
                return true;
            } else { // array_key > key
                hi = i - 1;
            }
        }
    }
    // Didn't find the key
    int64_t pred = hi;
    while (pred >= 0 && PMA_FUNC(empty_at_unchecked)(pma, (size_t)pred)) {
        pred--;
    }
    *index = pred;
    return false;
}

/*
* Walk up tree until reaching suitable interval (within threshold)
* Density: elements per chunk / slots in array (how much is occupied)
* Density thresholds: depends on which level of the tree we are. Strictest at top, looser at bottom
*
* For node at depth d in this conceptual binary tree, the density thresholds are:
* density >= 1/2 - 1/4 * d / h (interpolates between 1/4 when d == h and 1/2 when d == 0)
* density <= 3/4 + 1/4 * d / h (interpolates between 3/4 when d == 0 and 1 when d == h)
* For the root, density must be between 1/2 and 3/4
* For the leaf node, density can be between 1/4 and 1
*/

static size_t PMA_FUNC(size)(PMA_NAME *pma) {
    if (pma == NULL) return 0;
    return PMA_ARRAY_FUNC(size)(&pma->array);
}

static size_t PMA_FUNC(count)(PMA_NAME *pma) {
    if (pma == NULL) return 0;
    return pma->count;
}

static bool PMA_FUNC(pack)(PMA_NAME *pma, size_t from, size_t to, size_t n) {
    if (pma == NULL || from >= to || to > PMA_ARRAY_FUNC(size)(&pma->array)) {
        printf("PMA_FUNC(pack) failed: pma == NULL: %d || from >= to %d || to > PMA_ARRAY_FUNC(size)(&pma->array) %d\n", pma == NULL, from >= to, to > PMA_ARRAY_FUNC(size)(&pma->array));
        return false;
    }
    // In-place packing of elements in the range [from, to)
    // from is inclusive, to is exclusive
    size_t read_index = from;
    size_t write_index = from;
    while (read_index < to) {
        if (!PMA_FUNC(empty_at_unchecked)(pma, read_index)) {
            if (read_index > write_index) {
                PMA_ARRAY_FUNC(set_unchecked)(&pma->array, write_index, PMA_ARRAY_FUNC(get_unchecked)(&pma->array, read_index));
                PMA_ARRAY_FUNC(set_unchecked)(&pma->array, read_index, PMA_EMPTY_VALUE);
            }
            write_index++;
        }
        read_index++;
    }

    return (n == write_index - from);
}

static bool PMA_FUNC(spread)(PMA_NAME *pma, size_t from, size_t to, size_t n) {
    if (pma == NULL || from >= to || to > PMA_ARRAY_FUNC(size)(&pma->array)) {
        printf("PMA_FUNC(spread) failed: pma == NULL: %d || from >= to %d || to > PMA_ARRAY_FUNC(size)(&pma->array) %d\n", pma == NULL, from >= to, to > PMA_ARRAY_FUNC(size)(&pma->array));
        return false;
    }
    // In-place spreading of elements in the range [from, to) using 8-bit fixed-point arithmetic
    // from is inclusive, to is exclusive
    size_t capacity = to - from;
    size_t frequency = (capacity << 8) / n;
    size_t read_index = from + n - 1;
    size_t write_index = (to << 8) - frequency;
    while ((write_index >> 8) > read_index) {
        PMA_ARRAY_FUNC(set_unchecked)(&pma->array, write_index >> 8, PMA_ARRAY_FUNC(get_unchecked)(&pma->array, read_index));
        PMA_ARRAY_FUNC(set_unchecked)(&pma->array, read_index, PMA_EMPTY_VALUE);
        read_index--;
        write_index -= frequency;
    }
    return true;
}

static bool PMA_FUNC(resize)(PMA_NAME *pma) {
    if (pma == NULL) return false;
    if (!PMA_FUNC(pack)(pma, 0, PMA_ARRAY_FUNC(size)(&pma->array), pma->count)) {
        printf("PMA_FUNC(resize) failed: !PMA_FUNC(pack)(pma, 0, %zu, %zu)\n", PMA_ARRAY_FUNC(size)(&pma->array), pma->count);
        return false;
    }
    size_t segment_size = ceil_log2(pma->count);
    size_t num_segments = ceil_div(pma->count, segment_size);
    num_segments = hyper_ceil(num_segments);
    // Update the segment size accordingly
    segment_size = ceil_div(pma->count, num_segments);
    size_t new_capacity = PMA_LARGEST_MAX_SPARSITY * segment_size * num_segments;
    segment_size = PMA_LARGEST_MAX_SPARSITY * segment_size;

    if (new_capacity > PMA_MAX_SIZE) {
        printf("PMA_FUNC(resize) failed: new_capacity > PMA_MAX_SIZE %zu > %zu\n", new_capacity, PMA_MAX_SIZE);
        return false;
    }
    pma->segment_size = segment_size;
    pma->num_segments = num_segments;
    pma->height = floor_log2(num_segments) + 1;
    pma->density_high = (PMA_DENSITY_HIGH_LEAF - PMA_DENSITY_HIGH_ROOT) / pma->height;
    pma->density_low = (PMA_DENSITY_LOW_ROOT - PMA_DENSITY_LOW_LEAF) / pma->height;
    PMA_ARRAY_FUNC(resize_fixed)(&pma->array, new_capacity);
    for (size_t i = pma->count; i < new_capacity; i++) {
        PMA_ARRAY_FUNC(set_unchecked)(&pma->array, i, PMA_EMPTY_VALUE);
    }
    return PMA_FUNC(spread)(pma, 0, new_capacity, pma->count);
}


static bool PMA_FUNC(rebalance)(PMA_NAME *pma, size_t i) {
    int64_t window_start, window_end;
    uint8_t height = 0;

    uint64_t occupied = PMA_FUNC(empty_at)(pma, i) ? 0 : 1;
    int64_t left_index = ((int64_t)i) - 1;
    int64_t right_index = ((int64_t)i) + 1;
    double density, density_high, density_low;

    do {
        uint64_t window_size = pma->segment_size * (1 << height);
        uint64_t window = i / window_size;
        window_start = window * window_size;
        window_end = window_start + window_size;
        while (left_index >= window_start) {
            if (!PMA_FUNC(empty_at)(pma, (size_t)left_index)) {
                occupied++;
            }
            left_index--;
        }
        while (right_index < window_end) {
            if (!PMA_FUNC(empty_at)(pma, (size_t)right_index)) {
                occupied++;
            }
            right_index++;
        }
        density = (double)occupied / (double)window_size;
        density_high = PMA_DENSITY_HIGH_LEAF - (height * pma->density_high);
        density_low = PMA_DENSITY_LOW_LEAF + (height * pma->density_low);
        height++;
    } while (((density < density_low) ||
              (density > density_high) ||
              (fabs(density_high - density) < DBL_EPSILON)
            ) && (height < pma->height));

    // Found a window within threshold
    if ((density > density_low || fabs(density_high - density) < DBL_EPSILON) &&
        (density < density_high || fabs(density_low - density) < DBL_EPSILON)
    ) {
        if (!PMA_FUNC(pack)(pma, window_start, window_end, occupied)) {
            printf("PMA_FUNC(rebalance) failed: !PMA_FUNC(pack)(pma, %zu, %zu, %zu)\n", window_start, window_end, occupied);
            return false;
        }
        if (!PMA_FUNC(spread)(pma, window_start, window_end, occupied)) {
            printf("PMA_FUNC(rebalance) failed: !PMA_FUNC(spread)(pma, %zu, %zu, %zu)\n", window_start, window_end, occupied);
            return false;
        }
    } else {
        if (!PMA_FUNC(resize)(pma)) {
            printf("PMA_FUNC(rebalance) failed: !PMA_FUNC(resize)(pma)\n");
            return false;
        }
    }
    return true;
}



static bool PMA_FUNC(insert_after)(PMA_NAME *pma, int64_t i, PMA_TYPE key) {
    if (pma == NULL || i < -1 || (i >= 0 && PMA_FUNC(empty_at)(pma, (size_t)i))) {
        printf("PMA_FUNC(insert_after) failed: pma == NULL: %d || i < -1: %d || (i >= 0 && PMA_FUNC(empty_at)(pma, (size_t)i)) %d\n", pma == NULL, i < -1, i >= 0 && PMA_FUNC(empty_at)(pma, (size_t)i));
        return false;
    }
    // Find an empty space to the right of i, which should be close if density thresholds are met
    int64_t j = i + 1;

    size_t n = PMA_ARRAY_FUNC(size)(&pma->array);
    while(j < n && !PMA_FUNC(empty_at_unchecked)(pma, (size_t)j)) {
        j++;
    }
    if (j < n) {
        // Found an empty space to the right of i
        while (j > i + 1) {
           // Push elements to make space for the new value
            PMA_ARRAY_FUNC(set_unchecked)(&pma->array, (size_t)j, PMA_ARRAY_FUNC(get_unchecked)(&pma->array, (size_t)(j - 1)));
            j--;
        }
        PMA_ARRAY_FUNC(set_unchecked)(&pma->array, (size_t)i + 1, key);
        i++;
    } else {
        // No empty space to the right of i
        j = i - 1;
        while (j >= 0 && !PMA_FUNC(empty_at_unchecked)(pma, (size_t)j)) {
            j--;
        }
        if (j >= 0) {
            // Found a space
            while (j < i) {
                // Push elements over to make space for the new value
                PMA_TYPE next = PMA_ARRAY_FUNC(get_unchecked)(&pma->array, (size_t)(j + 1));
                PMA_ARRAY_FUNC(set_unchecked)(&pma->array, (size_t)j, next);
                j++;
            }
            PMA_ARRAY_FUNC(set_unchecked)(&pma->array, (size_t)i, key);
        }
    }

    pma->count++;
    return PMA_FUNC(rebalance)(pma, (size_t)i);
}

static bool PMA_FUNC(insert)(PMA_NAME *pma, PMA_TYPE key) {
    if (pma == NULL) return false;
    int64_t i;
    if (!PMA_FUNC(find)(pma, key, &i)) {
        // No duplicates
        return PMA_FUNC(insert_after)(pma, i, key);
    }
    return false;
}

static bool PMA_FUNC(delete_at)(PMA_NAME *pma, size_t i) {
    if (pma == NULL) {
        return false;
    }
    if (!PMA_ARRAY_FUNC(set)(&pma->array, i, PMA_EMPTY_VALUE)) {
        return false;
    }
    pma->count--;
    // Rebalance the tree if necessary
    return PMA_FUNC(rebalance)(pma, i);
}

static bool PMA_FUNC(delete)(PMA_NAME *pma, PMA_TYPE value) {
    if (pma == NULL) return false;

    int64_t i;
    if (PMA_FUNC(find)(pma, value, &i)) {
        if (!PMA_FUNC(delete_at)(pma, (size_t)i)) return false;
        return true;
    }
    return false;
}

void PMA_FUNC(destroy)(PMA_NAME *pma) {
    if (pma == NULL) return;
    PMA_ARRAY_FUNC(destroy_data)(&pma->array);
    PMA_FREE(pma);
}


#ifdef PMA_LESS_THAN_DEFINED
#undef PMA_LESS_THAN
#undef PMA_LESS_THAN_DEFINED
#endif

#ifdef PMA_EQUALS_DEFINED
#undef PMA_EQUALS
#undef PMA_EQUALS_DEFINED
#endif

#ifdef PMA_EMPTY_VALUE_DEFINED
#undef PMA_EMPTY_VALUE
#undef PMA_EMPTY_VALUE_DEFINED
#endif

#ifdef PMA_MALLOC_DEFINED
#undef PMA_MALLOC
#undef PMA_MALLOC_DEFINED
#endif

#ifdef PMA_FREE_DEFINED
#undef PMA_FREE
#undef PMA_FREE_DEFINED
#endif

#ifdef PMA_ARRAY_DEFINED
#undef PMA_ARRAY_DEFINED
#undef PMA_ARRAY_NAME
#endif

#undef PMA_ARRAY_FUNC
#undef PMA_ARRAY_TYPE